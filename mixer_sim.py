#!/usr/bin/python3

import skrf as rf
from skrf.calibration import TwelveTerm
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import copy
import remove_noise

def pTodB(p):
    return 10*np.log10(p)

def dBToP(db):
    return 10**(db/10)

def vTodB(v):
    return 20*np.log10(p)

def dBToV(db):
    return 10**(db/20)

class Spectrum:
    def __init__(self, f, s, name='signal spectrum', watts=False, volts=False, z0=50, fmt='', color=None):
        self.z0 = z0
        self.fmt=fmt
        self.color=color
        self.f = np.array(f)
        self.components = []
        if watts:
            self.s = np.array(s)
        if volts:
            self.s = (np.array(s)**2)/self.z0
        else:
            self.s=np.array(s)
        self.name = name


    def scale_power(self, sf):
        self.s*=sf

    def plot_y(self, y):
        #print("plotting", self.name, "color", self.color)
        plt.plot(self.f, y, self.fmt, label=self.name, color=self.color)

    def plot_components_db(self, show=True):
        for c in self.components:
            c.plot_db(show=False)
        if show:
            plt.show()

    def plot_components_dbm(self, show=True):
        for c in self.components:
            c.plot_dbm(show=False)
        if show:
            plt.show()

    def plot_db(self, show=True):
        self.plot_y(10*np.log10(self.s))
        plt.grid(True)
        plt.title('Normalized power spectrum')
        plt.ylabel('Relative power [dB]')
        plt.xlabel('Frequency [Hz]')
        if show:
            plt.show()

    def plot_dbm(self, show=True):
        self.plot_y(10*np.log10(self.s)+30)
        plt.grid(True)
        plt.title('Power spectrum')
        plt.ylabel('Relative power [dBm]')
        plt.xlabel('Frequency [Hz]')
        if show:
            plt.show()

    def plot_power(self, show=True):
        self.plot_y(self.s)
        plt.grid(True)
        plt.title('Power spectrum')
        plt.ylabel('Power [W]')
        plt.xlabel('Frequency [Hz]')
        if show:
            plt.show()

    def plot_voltage(self, show=True):
        self.plot_y(np.sqrt(self.s*self.z0))
        plt.grid(True)
        plt.title('Amplitude spectrum')
        plt.ylabel('Voltage [V]')
        plt.xlabel('Frequency [Hz]')
        if show:
            plt.show()

    
    def frequency_scaled(self, scale):
        new = copy.deepcopy(self)
        new.f = new.f*scale
        if scale != 1:
            new.name = "H%d(%s)"%(scale, new.name)
        new.components=[]
        for c in self.components:
            new.components.append(c.frequency_scaled(scale))
        return new

    def power_scaled(self, scale):
        new = copy.deepcopy(self)
        new.s = new.s*scale
        new.name = 'dB(%d, %s)'%(pTodB(scale), new.name)
        new.components=[]
        for c in self.components:
            new.components.append(c.power_scaled(scale))
        return new

    def with_harmonics(self, levels):
        fmt = list(Line2D.markers.keys())[3:]
        out = copy.deepcopy(self)
        out.s = np.zeros(out.s.shape)
        h = 1
        coponents = [];
        for l in levels:
            harmonic = self.frequency_scaled(h).power_scaled(l)
            harmonic.fmt = fmt[h-1]
            out += harmonic
            h+=1
        return out

    def frequency_shifted(self, f, name=None):
        fnew = self.f+f
        neg = np.where(fnew < 0)
        pos = np.where(fnew >= 0)
        fneg = np.flip(np.abs(fnew[neg]))
        fpos = fnew[pos]

        fout = np.concatenate((fpos, fneg))
        fout.sort(kind='mergesort')
        snew = np.zeros(fout.shape)
        new = copy.deepcopy(self)
        new.components=[]
        if len(fpos):
            specpos = copy.deepcopy(self)
            specpos.s = self.s[pos]
            specpos.f = fpos
            if name:
                specpos.name ='(%s)%s'%(self.name, name);
            else:
                specpos.name ='(%s)%+.2f[GHz]'%(self.name, f/1e9);
            new.components.append(specpos)
            snew += np.interp(fout, fpos, self.s[pos], left=0, right=0)
        if len(fneg):
            specneg = copy.deepcopy(self)
            specneg.s = np.flip(self.s[neg])
            specneg.f = fneg
            if name:
                specneg.name ='M((%s)%s)'%(self.name, name);
            else:
                specneg.name ='M((%s)%+.2f[GHz])'%(self.name, f/1e9);
            specneg.fmt+='|'
            new.components.append(specneg)
            snew += np.interp(fout, fneg, np.flip(self.s[neg]), left=0, right=0)
        new.s = snew
        new.f = fout
        if name:
            new.name ='(%s)%s'%(self.name, name);
        else:
            new.name ='(%s)%+.2f[GHz]'%(self.name, f/1e9);
        for c in self.components:
            new.components.append(c.frequency_shifted(f, name=name))
        return new

    def __add__(self, other):
        new = copy.deepcopy(self)
        less = np.where(other.f < np.min(self.f))
        more = np.where(other.f > np.max(self.f))
        headF = other.f[less]
        headS = other.s[less]
        tailF = other.f[more]
        tailS = other.s[more]
        new.s += np.interp(new.f, other.f, other.s, left=0, right=0)
        new.s = np.concatenate([headS, new.s, tailS])
        new.f = np.concatenate([headF, new.f, tailF])
        new.name += ' + '+other.name
        if len(other.components):
            new.components+=other.components
        else:
            new.components.append(other)
        return new

    def __mul__(self, other):
        if not isinstance(other, rf.Network):
            raise ValueError("Spectrum multiplication operator is in tended to be used for filtering so right hand side must be skrf.Network")
        new = copy.deepcopy(self)
        fo = other.f;
        amp = np.abs(other.s[:,0,1])
        new.s *= np.interp(new.f, fo, amp**2)
        new.name = "(%s)*%s"%(new.name, other.name)
        new.components=[]
        for c in self.components:
            new.components.append(c*other)
        return new

    def normalized(self):
        new = copy.deepcopy(self)
        sf = 1/np.max(np.abs(self.s))
        new.s*=sf
        new.components=[]
        for c in self.components:
            new.components.append(c)
            new.components[-1].s*=sf
        return new
        

def spectrum_from_network(n, ref_level_dbm=0, ref_impedance=50, **kwargs):
    dBmVoltage = np.sqrt(0.001*ref_impedance)
    refFactor = dBmVoltage*10**(ref_level_dbm/10)
    return Spectrum(n.f, np.abs(n.s[:,1,0])*refFactor, n.name, volts=True, **kwargs)

def cw_spectrum(f, watts, pfx='', **kwargs):
    return Spectrum([f-1, f, f+1, ], [watts*0.1, watts, watts*0.1], pfx+'cw %.2f[GHz]'%(f/1e9), **kwargs);


class MixerModel:
    def __init__(self, lo_harmonics, rf_harmonics):
        self.lo_loss = (np.array([0, 27, 12, 33, 22])+12)[:lo_harmonics]
        self.rf_loss = dBToP(np.array([0, -40, -60, -80, -80]))[:rf_harmonics]
        #self.lo_loss = np.array([0])+12
        self.loIfIsolation=25;
        self.loRfIsolation=25;
        self.rfIfIsolation=25;

    def computeIfSpectrum(self, rfSpectrum, fLo, loPower):
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']
        rfSpectrum.color=colors[-1]
        outS = copy.deepcopy(rfSpectrum)
        outS.s = np.zeros(outS.s.shape)
        outS+= rfSpectrum.power_scaled(dBToP(-self.rfIfIsolation))

        for l, n in zip(self.lo_loss, range(len(self.lo_loss))):
            loh = (n+1)*fLo
            loS =  cw_spectrum(fLo, loPower*dBToP(-self.rfIfIsolation), pfx='lo ', color=colors[n], fmt='-x')
            #outS += cw_spectrum(loh, loPower*dBToP(-self.rfIfIsolation), color=colors[n], fmt='-x')
            if n > 0:
                loS = loS.frequency_scaled(n+1)
            outS+=loS
                
            cs = copy.deepcopy(rfSpectrum)
            cs.color = colors[n]
            #rfWithHarmonics = cs.with_harmonics(dBToP(np.array([0, -40, -60, -80])))
            #rfWithHarmonics = cs.with_harmonics(dBToP(np.array([0, -40, -60])))
            #rfWithHarmonics = cs.with_harmonics(dBToP(np.array([0])))
            rfWithHarmonics = cs.with_harmonics(self.rf_loss)
            lsb = rfWithHarmonics.frequency_shifted(-loh, name="%dLO"%(-(n+1)))
            usb = rfWithHarmonics.frequency_shifted(loh, name="%dLO"%(n+1))
            lsb = lsb.power_scaled(dBToP(-l))
            usb = usb.power_scaled(dBToP(-l))
            for c in usb.components:
                c.fmt='--'+c.fmt
            for l in lsb.components:
                l.fmt='-'+l.fmt
            lsb.color=colors[n]
            usb.color=colors[n]
            outS+=lsb
            outS+=usb
            
        outS.name = "mixer output"
        return outS

if __name__ == "__main__":
    import argparse
    import filter_identification

    def loadNetworkList(fileList):
        allNetwork = []
        for file in fileList:
            allNetwork.append(rf.Network(file))
        return allNetwork;
        
    parser = argparse.ArgumentParser(description="Calculate harmonic transfer functions for mixers")
    parser.add_argument('--rf-filter', '-rffi', help='S2P file containing rf (preselector) filter. S21 will be used', nargs = '*')
    parser.add_argument('--rf-frequency', '-rffr', help='Single tone representing RF signal',type = float)
    parser.add_argument('--if-filter', '-iff', help='S2P file containing if (roofing) filter. S21 will be used', nargs='*')
    parser.add_argument('--lo-frequency', '-lof', help="Frequency of LO", type=float)
    parser.add_argument('--lo-auto-lsb', '-lal', help="Automatically determine LO frequency so rf lsb will end up in middle of IF-filter", action='store_true')
    parser.add_argument('--lo-auto-usb', '-lau', help="Automatically determine LO frequency so rf usb will end up in middle of IF-filter", action='store_true')
    parser.add_argument('--lo-harmonics', '-lh', help="Number of lo harmonics to consider, default = 5", type=int, default=5, choices=range(1,6))
    parser.add_argument('--rf-harmonics', '-rh', help="Number of rf harmonics to consider, default = 3", type=int, default=3, choices=range(1,6))
    parser.add_argument('--no-legend', '-nl', help="Do not display legend in plots containing spectral components", action='store_true')
    args = parser.parse_args()
    ifFilters = None
    if args.if_filter:
        ifFilters = loadNetworkList(args.if_filter)
        ifFilter = rf.cascade_list(ifFilters)
        #ifFilter = rf.Network(args.if_filter.name)
        ifFreq = filter_identification.getCfBw(filter_identification.identifyBands(ifFilter)[0], ifFilter.f)[0]
        ifFilter.name='if filter'

    if args.rf_filter and args.rf_frequency:
        print("Specify --rf-filter OR -rf-frequency")
        quit()

    if args.rf_frequency:
        rfSpectrum = cw_spectrum(args.rf_frequency, 1e-3, "rf ")
        rfFreq = args.rf_frequency

    rfFilter = None
    rfFilters = None
    if args.rf_filter:
        rfFilters = loadNetworkList(args.rf_filter)
        rfFilter = rf.cascade_list(rfFilters)
        rfFilter.name='rf filter'
        rfSpectrum = spectrum_from_network(rfFilter)
        rfFreq = filter_identification.getCfBw(filter_identification.identifyBands(rfFilter)[0], rfFilter.f)[0]

    if args.lo_frequency:
        lof = args.lo_frequency
    if args.lo_auto_usb:
        lof = abs(ifFreq-rfFreq)
    if args.lo_auto_lsb:
        lof = ifFreq+rfFreq
    m = MixerModel(args.lo_harmonics, args.rf_harmonics)

    plt.figure();
    cw_spectrum(lof, 0.001, 'lo ').plot_dbm(show=False)
    if rfFilter is None:
        rfSpectrum.plot_dbm(show=False)
    else:
        if len(rfFilters) > 1:
            for rff in rfFilters:
                rff.name = 'rf '+rff.name
                rff.plot_s_db(m=1, n=0)
        rfFilter.plot_s_db(m=1, n=0)
    if len(ifFilters) > 1:
        for iff in ifFilters:
            iff.name = 'if '+iff.name
            iff.plot_s_db(m=1, n=0)
    ifFilter.plot_s_db(m=1, n=0)
    plt.grid(True)
    plt.title('input parameters')
        
    ifSpectrum = m.computeIfSpectrum(rfSpectrum, lof, dBToP(13-30))
    plt.figure()
    ifSpectrum.plot_components_dbm(show=False)
    ifFilter.plot_s_db(m=1,n=0, marker='d')
    if not args.no_legend:
        plt.legend(loc='right');
    else:
        plt.legend('', frameon=False)
    plt.title('Components for lo %e'%(lof))

    ifFiltered = ifSpectrum*ifFilter
    plt.figure()
    ifFiltered.normalized().plot_components_db(show=False)
    ifFiltered.fmt='d'
    ifFiltered.normalized().plot_db(show=False)
    if not args.no_legend:
        plt.legend(loc='right');
    else:
        plt.legend('', frameon=False)
    plt.xlim([0, 25e9])
    plt.title('Normalized IF spectrun filtered by if filter')

    plt.show();
