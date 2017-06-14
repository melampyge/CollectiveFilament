#!/usr/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import sys
import math
import os
from scipy import linalg

try:
    ampfile = sys.argv[1]
    evfile = sys.argv[2]
    acfile = sys.argv[3]
    ccfile = sys.argv[4]
except:
    print 'Usage: ' + sys.argv[0] + '           amplitude evolution file        eigenvalue file        acfile         ccfile'
    exit()


##################################################################

def read_amplitudes(ampfile):
    """ read in the amplitudes"""
    t = []
    amp1 = []
    amp2 = []
    amp3 = []
    amp4 = []

    ifile = open(ampfile)
    for i in range(3):
        ifile.readline()
    for line in ifile:
        line = line.split()
        t.append(float(line[0]))
        amp1.append(float(line[1]))
        amp2.append(float(line[2]))
        amp3.append(float(line[3]))
        amp4.append(float(line[4]))
    ifile.close()

    t = np.array(t, dtype = np.int32)
    amp1 = np.array(amp1, dtype = np.float64)
    amp2 = np.array(amp2, dtype = np.float64)
    amp3 = np.array(amp3, dtype = np.float64)
    amp4 = np.array(amp4, dtype = np.float64)
    
    return t, amp1, amp2, amp3, amp4

##################################################################

def read_eigenvalues(evfile):
    """ read in the eigenvalues"""
    ifile = open(evfile)
    for i in range(3):
        ifile.readline()
    line = ifile.readline()
    line = line.split()
    ev = []
    for i in range(1, len(line)):
        ev.append(float(line[i]))
    ifile.close()
    ev = np.array(ev, dtype = np.float64)
    return ev

##################################################################

def read_cc(ccfile):
    """ read in the autocorrelation functions"""
    ifile = open(ccfile)
    for i in range(3):
        ifile.readline()
    t = []
    c12 = []
    sc12 = []
    c13 = []
    sc13 = []
    c14 = []
    sc14 = []
    c23 = []
    sc23 = []
    c24 = []
    sc24 = []
    c34 = []
    sc34 = []
    for line in ifile:
        line = line.split()
        t.append(int(line[0]))
        c12.append(float(line[1]))
        sc12.append(float(line[2]))
        c13.append(float(line[3]))
        sc13.append(float(line[4]))
        c14.append(float(line[5]))
        sc14.append(float(line[6]))
        c23.append(float(line[7]))
        sc23.append(float(line[8]))
        c24.append(float(line[9]))
        sc24.append(float(line[10]))
        c34.append(float(line[11]))
        sc34.append(float(line[12]))
    ifile.close()
    t = np.array(t, dtype = np.int32)
    c12 = np.array(c12, dtype = np.float64)
    sc12 = np.array(sc12, dtype = np.float64)
    c13 = np.array(c13, dtype = np.float64)
    sc13 = np.array(sc13, dtype = np.float64)
    c14 = np.array(c14, dtype = np.float64)
    sc14 = np.array(sc14, dtype = np.float64)
    c23 = np.array(c23, dtype = np.float64)
    sc23 = np.array(sc23, dtype = np.float64)
    c24 = np.array(c24, dtype = np.float64)
    sc24 = np.array(sc24, dtype = np.float64)
    c34 = np.array(c34, dtype = np.float64)
    sc34 = np.array(sc34, dtype = np.float64)



    return t, c12, sc12, c13, sc13, c14, sc14, c23, sc23, c24, sc24, c34, sc34

##################################################################

def read_ac(acfile):
    """ read in the autocorrelation functions"""
    ifile = open(acfile)
    for i in range(3):
        ifile.readline()
    t = []
    a1 = []
    sa1 = []
    a2 = []
    sa2 = []
    a3 = []
    sa3 = []
    a4 = []
    sa4 = []
    for line in ifile:
        line = line.split()
        t.append(int(line[0]))
        a1.append(float(line[1]))
        sa1.append(float(line[2]))
        a2.append(float(line[3]))
        sa2.append(float(line[4]))
        a3.append(float(line[5]))
        sa3.append(float(line[6]))
        a4.append(float(line[7]))
        sa4.append(float(line[8]))
    ifile.close()
    t = np.array(t, dtype = np.int32)
    a1 = np.array(a1, dtype = np.float64)
    sa1 = np.array(sa1, dtype = np.float64)
    a2 = np.array(a2, dtype = np.float64)
    sa2 = np.array(sa2, dtype = np.float64)
    a3 = np.array(a3, dtype = np.float64)
    sa3 = np.array(sa3, dtype = np.float64)
    a4 = np.array(a4, dtype = np.float64)
    sa4 = np.array(sa4, dtype = np.float64)

    return t, a1, sa1, a2, sa2, a3, sa3, a4, sa4

##################################################################

def multivariate_histogram(x,y,fname):
    """ compute multivariate histograms"""
    ### generate a histogram of amp1 and amp1
    fig = plt.figure()

    nullfmt = NullFormatter()
    
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_hist2d = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    ax = plt.axes(rect_hist2d)
    ax1 = plt.axes(rect_histx)
    ax2 = plt.axes(rect_histy)

    ax1.xaxis.set_major_formatter(nullfmt)
    ax2.yaxis.set_major_formatter(nullfmt)

    nbins = 50
    ax.hist2d(x, y, bins = nbins, cmap = 'copper')
    ax1.hist(x, bins = nbins)
    ax2.hist(y, bins = nbins, orientation = 'horizontal')
    
    plt.savefig(fname)
    plt.close()

    return
    
##################################################################
    
def main():
    """ main function, called when the script is started"""
    ### read in the amplitutes
    t, amp1, amp2, amp3, amp4 = read_amplitudes(ampfile)

    ### read in the eigenvalues
    ev = read_eigenvalues(evfile)

    ### read in the autocorrelations
    tc, a1, sa1, a2, sa2, a3, sa3, a4, sa4 = read_ac(acfile)
    ### read in the crosscorrelations
    tc, c12, sc12, c13, sc13, c14, sc14, c23, sc23, c24, sc24, c34, sc34 = read_cc(ccfile)


    ### plot eigenvalues
    fig = plt.figure()
    ax = plt.subplot()
    ax.semilogy(ev, ls = ':', marker = 'o')
    ax.set_ylim([4e-3, 1.1e-1])
    plt.savefig('ev.png')
    plt.close()
    
    ### generate multivariate histograms
    multivariate_histogram(amp1, amp2, 'a1a2.png')
    multivariate_histogram(amp1, amp3, 'a1a3.png')
    multivariate_histogram(amp1, amp4, 'a1a4.png')
    multivariate_histogram(amp2, amp3, 'a2a3.png')
    multivariate_histogram(amp2, amp4, 'a2a4.png')
    multivariate_histogram(amp4, amp2, 'a4a2.png')
    multivariate_histogram(amp3, amp4, 'a3a4.png')

    ### plot autocorrelation functions
    fig = plt.figure()
    ax1 = plt.subplot(411)
    ax2 = plt.subplot(412)
    ax3 = plt.subplot(413)
    ax4 = plt.subplot(414)

    ax1.errorbar(tc, a1, yerr = 2*sa1, ls = ':', marker = 'o', markersize = 2)
    ax2.errorbar(tc, a2, yerr = 2*sa2, ls = ':', marker = 'o', markersize = 2)
    ax3.errorbar(tc, a3, yerr = 2*sa3, ls = ':', marker = 'o', markersize = 2)
    ax4.errorbar(tc, a4, yerr = 2*sa4, ls = ':', marker = 'o', markersize = 2)

    plt.tight_layout()
    
    plt.savefig('ac.png')
    plt.close()


    # plot crosscorrelation functions
    fig = plt.figure()
    ax1 = plt.subplot(611)
    ax2 = plt.subplot(612)
    ax3 = plt.subplot(613)
    ax4 = plt.subplot(614)
    ax5 = plt.subplot(615)
    ax6 = plt.subplot(616)

    ax1.errorbar(tc, c12, yerr = 2*sc12, ls = ':', marker = 'o', markersize = 2)
    ax2.errorbar(tc, c13, yerr = 2*sc13, ls = ':', marker = 'o', markersize = 2)
    ax3.errorbar(tc, c14, yerr = 2*sc14, ls = ':', marker = 'o', markersize = 2)
    ax4.errorbar(tc, c23, yerr = 2*sc23, ls = ':', marker = 'o', markersize = 2)
    ax5.errorbar(tc, c24, yerr = 2*sc24, ls = ':', marker = 'o', markersize = 2)
    ax6.errorbar(tc, c34, yerr = 2*sc34, ls = ':', marker = 'o', markersize = 2)

    plt.tight_layout()

    plt.savefig('cc.png')
    plt.close()    

     
    return
    
    
##################################################################

if __name__ == '__main__':
    main()
    
