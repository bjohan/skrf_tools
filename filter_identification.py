#!/usr/bin/python3
import skrf as rf
from matplotlib import pyplot as plt
import numpy as np
import copy
import remove_noise

def pTodB(p):
    return 10*np.log10(p)

def dBToP(db):
    return 10**(db/10)

def vTodB(v):
    return 20*np.log10(np.abs(v))

def dBToV(db):
    return 10**(db/20)


def splitBands(idx):
    #print(idx)
    bands = []
    bandEdge = np.where(np.squeeze(np.diff(idx))>1)[0]
    #print("bandEdge", bandEdge)
    if len(bandEdge) == 0:
        bands.append(idx)
    else:
        start = 0
        for be in bandEdge:
            #print("band from", start, "to", idx[be])
            bands.append(idx[start:be+1])
            start = be+1;
        #print("band from", start, "to", idx[-1])
        bands.append(idx[start:])
    bret = []
    for b in bands:
        if len(b):
            bret.append(b)
    return bret

def getFrequencies(b, fax):
    return fax[b][0], fax[b][-1]

def getCfBw(b, fax):
    f1, f2 = getFrequencies(b, fax)
    return 0.5*(f1+f2), f2-f1

def mergableBands(ba, bb, fax, dist):
    af1, af2 = getFrequencies(ba, fax)
    acf, abw = getCfBw(ba, fax)

    bf1, bf2 = getFrequencies(bb, fax)
    bcf, bbw = getCfBw(bb, fax)

    f = 1+dist

    if ((abw+bbw)*f)/2 > np.abs(bcf-acf): #Todo, more precise check for overlap 
        return True

    return False


def mergeAdjacentBands(bands, fax, dist):
    nb = len(bands)
    for i in range(nb):
        for j in range(nb):
            if not(i == j):
                #print("band %d %d"%(i,j), mergableBands(bands[i], bands[j], fax, 0.01))
                if mergableBands(bands[i], bands[j], fax, dist):
                    #print(np.concatenate((bands[i],bands[j])))
                    low = np.min(np.concatenate((bands[i],bands[j])))
                    high = np.max(np.concatenate((bands[i],bands[j])))
                    nb=np.arange(low, high+1)
                    bands[i]=nb
                    del bands[j]
                    return mergeAdjacentBands(bands, fax, dist)
    return bands
                    
def filterType(bands, f):
    if len(bands) == 0:
        return "isolation"
    if len(bands) > 1:
        return "multiband"
    if (bands[0][0] == 0) and (bands[0][-1] == len(f)-1):
        return "through"
    if bands[0][0] == 0:
        return "lowpass"
    if bands[0][-1] == len(f)-1:
        return "high pass"
    return "bandpass"

def bandType(band, f):
    if len(band)==0:
        return "empty"
    if (band[0] == 0) and (band[-1] == len(f)-1):
        return "through"
    if band[0] == 0:
        return "lowpass"
    if band[-1] == len(f)-1:
        return "high pass"
    return "bandpass"


def identifyBands(n, passbandThreshold=-3, isolationThreshold=-10, mergeWidth=0.05):
    f = n.f
    mag = vTodB(np.squeeze(n.s[:,1,0]))
    maxMag = max(mag)
    if maxMag < isolationThreshold:
        return []
    pb = np.where(mag > maxMag+passbandThreshold)[0]
    bands = splitBands(pb)
    return mergeAdjacentBands(bands, f, mergeWidth)

def describeBand(band, n):
    bt = bandType(band, n.f)
    if bt == "through":
        return "Through, loss %.2f[dB]"%(np.min(vTodB(np.squeeze(n.s[:,1,0]))))
    if bt == "bandpass":
        cf, bw = getCfBw(band, n.f)
        return "Bandpass fc %.2f[MHz] bw %.2f[MHz]"%(cf/1e6, bw/1e6)
    if bt == "lowpass":
        return "Lowpass cf %.2f[MHz]"%(getFrequencies(band, n.f)[1]/1e6)
    if bt == "high pass":
        return "Highpass cf %.2f[MHz]"%(getFrequencies(band, n.f)[0]/1e6)
    
def describeFilter(n, passbandThreshold=-3, isolationThreshold=-10, mergeWidth=0.05):
    bands = identifyBands(n, passbandThreshold, isolationThreshold, mergeWidth)
    ft = filterType(bands, n.f)
    if ft == "through":
        return "Through, loss %.2f[dB]"%(np.min(vTodB(np.squeeze(n.s[:,1,0]))))
    if ft == "isolation":
        return "Isolation %.2f[dB]"%(np.max(vTodB(np.squeeze(n.s[:,1,0]))))
    if ft == "bandpass":
        cf, bw = getCfBw(bands[0], n.f)
        return "Bandpass filter. Center frequency %.2f[MHz] bandwidth %.2f[MHz]"%(cf/1e6, bw/1e6)
    if ft == "lowpass":
        return "Lowpass filter. 3dB frequency %.2f[MHz]"%(getFrequencies(bands[0], n.f)[1]/1e6)
    if ft == "high pass":
        return "Highpass filter. 3dB frequency %.2f[MHz]"%(getFrequencies(bands[0], n.f)[0]/1e6)
    if ft == "multiband":
        s = "Multi band filter."
        for b in bands:
            cf, bw = getCfBw(b, n.f)
            s+='\n\t'+describeBand(b, n)
        return s;
    return ft

def identifyFilterType(n, plot = False, passbandThreshold=-3, isolationThreshold=-10, mergeWidth=0.05):
    print(n.name)
    print('\t', describeFilter(n, passbandThreshold, isolationThreshold, mergeWidth))
    if plot:
        f = n.f
        mag = vTodB(np.squeeze(n.s[:,1,0]))
        bands = identifyBands(n, passbandThreshold, isolationThreshold, mergeWidth)
        plt.plot(f, mag)
        for b in bands:
            if len(b):
                f1=f[b[0]]/1e6
                f2=f[b[-1]]/1e6
                plt.plot(f[b], mag[b], label=describeBand(b, n))
            else:
                print("Warning", b)
        plt.legend()
        plt.grid(True)
        plt.title(n.name+" "+filterType(bands, n.f))
        plt.show()

if __name__ == "__main__":
    import argparse

    def loadNetworkList(fileList):
        allNetwork = []
        for file in fileList:
            allNetwork.append(rf.Network(file))
        return allNetwork;
        
    parser = argparse.ArgumentParser(description="Automatic identification of filters in s-parameters")
    parser.add_argument('touchstone',  help='S2P files to identify', nargs = '*')
    parser.add_argument('--passband-threshold', '-pt', help="Threshold for determining passband relative to max value of S21. Default -3dB", type=float, default=-3)
    parser.add_argument('--isolation-threshold', '-it', help="Threshold for determining isolation relative to max value of s21. Default -10dB", type=float, default=-10)
    parser.add_argument('--merge-distance', '-md', help="For denoising. If passbands are closer (1+md)*bw they are merged in to one. Default 0.05", type=float, default=0.05)
    parser.add_argument('--plot', '-p', help="Plot filter data", action='store_true')
    args = parser.parse_args()

    for f in args.touchstone:
        n = rf.Network(f)
        identifyFilterType(n, plot=args.plot, passbandThreshold=args.passband_threshold, isolationThreshold=args.isolation_threshold, mergeWidth=args.merge_distance)
