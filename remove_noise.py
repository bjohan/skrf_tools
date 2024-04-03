#!/usr/bin/python3

import skrf as rf
from skrf.calibration import TwelveTerm
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import copy

def detect_signal(v, filter_length=50, threshold=4, dilation=0.1):
    d = np.abs(np.diff(np.diff(20*np.log10(np.abs(v)))))
    davg = np.convolve(d, np.ones(filter_length)/filter_length, mode='same')
    valid = np.where(davg<threshold)
    numValid = np.zeros(davg.shape)
    numValid[valid]=1.0
    validAvg = np.convolve(numValid, np.ones(filter_length)/filter_length, mode='same')
    valid = np.where(validAvg>dilation) #dilate and smooth
    return valid

def detect_noise(v, filter_length=50, threshold=4, dilation=0.1):
    d = np.abs(np.diff(np.diff(20*np.log10(np.abs(v)))))
    davg = np.convolve(d, np.ones(filter_length)/filter_length, mode='same')
    valid = np.where(davg<threshold)
    numValid = np.zeros(davg.shape)
    numValid[valid]=1.0
    validAvg = np.convolve(numValid, np.ones(filter_length)/filter_length, mode='same')
    valid = np.where(validAvg<dilation) #dilate and smooth
    return valid

def clean_network(net, eps=10**(-130/20), filter_length=50, threshold=4, dilation=0.1):
        for m in range(2):
            for n in range(2):
                noiseSamples = detect_noise(net.s[:, m,n], filter_length=filter_length, threshold=threshold, dilation=dilation)
                signalSamples = detect_signal(net.s[:, m,n], filter_length=filter_length, threshold=threshold, dilation=dilation)
                net.s[noiseSamples, m,n] = eps
        return net
            

if __name__ == "__main__":
    class ArgRange(object):
        def __init__(self, start, end):
            self.start = start
            self.end = end
        def __eq__(self, other):
            return self.start <= other <= self.end

        def __str__(self):
            return "[%f-%f]"%(self.start, self.end)

    import argparse
    parser = argparse.ArgumentParser(description="Condition s2p-files by removing trace noise du to low signal level")
    parser.add_argument('input', help='S2P file to condition',type = argparse.FileType("r"))
    parser.add_argument('--write', '-w', help='Write conditioned data to this s2p-file', type = argparse.FileType('wb'))
    parser.add_argument('--epsilon','-e',  help='Replace noise with this value in dB. Default -130 dB', type = float, nargs='?', default=-130)
    parser.add_argument('--filter-length', '-f', help='Number of taps in averaging filter. Default 50 taps', type = int, nargs='?', default=50)
    parser.add_argument('--threshold', '-t', help='Threshold in dB for noise. Applied to second derivative of s-parameter filtered using filter above. Default 4 dB ', type = float, nargs='?', default=4)
    parser.add_argument('--dilation', '-d', help='Dilation threshold to make contiguous noise regions. Default 0.1', type = float, nargs='?', default=0.1, choices=[ArgRange(0.0, 1.0)])
    parser.add_argument('--plot', '-p', help='Plot result', action='store_true')
    args = parser.parse_args()

    n = clean_network(rf.Network(args.input.name), eps=10**(args.epsilon/20), filter_length=args.filter_length, threshold = args.threshold, dilation=args.dilation)
    if args.plot:
        n.plot_s_db();
        plt.show();
    if args.write:
        print("Writing conditioned data to", args.write.name)
        n.write_touchstone(args.write.name)
    

