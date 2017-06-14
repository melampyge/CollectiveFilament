#!/usr/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import sys
import math
import os
from scipy import linalg
from scipy import stats

try:
    evfile = sys.argv[1]
    acfile = sys.argv[2]
    ampfile = sys.argv[3]
except:
    print 'Usage: ' + sys.argv[0] + '       eigenvalue file         amplitude correlation file             amplitude evolution file'
    exit()

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

def read_correlations(acfile):
    """ read in the correlation functions"""
    t = []
    ac1 = []
    ac1_std = []
    ac2 = []
    ac2_std = []
    a1a2 = []
    a1a2_std = []
    a2a1 = []
    a2a1_std = []

    ifile = open(acfile)
    for i in range(3):
        ifile.readline()
    for line in ifile:
        line = line.split()
        t.append(int(line[0]))
        ac1.append(float(line[1]))
        ac1_std.append(float(line[2]))
        ac2.append(float(line[3]))
        ac2_std.append(float(line[4]))
        a1a2.append(float(line[5]))
        a1a2_std.append(float(line[6]))
        a2a1.append(float(line[7]))
        a2a1_std.append(float(line[8]))
    ifile.close()

    t = np.array(t, dtype = np.int32)
    ac1 = np.array(ac1, dtype = np.float64)
    ac1_std = np.array(ac1_std, dtype = np.float64)
    ac2 = np.array(ac2, dtype = np.float64)
    ac2_std = np.array(ac2_std, dtype = np.float64)
    a1a2 = np.array(a1a2, dtype = np.float64)
    a1a2_std = np.array(a1a2_std, dtype = np.float64)
    a2a1 = np.array(a2a1, dtype = np.float64)
    a2a1_std = np.array(a2a1_std, dtype = np.float64)

    return t, ac1, ac1_std, ac2, ac2_std, a1a2, a1a2_std, a2a1, a2a1_std

##################################################################

def read_amplitudes(ampfile):
    """ read in the amplitudes"""
    t = []
    amp = []
    ifile = open(ampfile)
    for i in range(3):
        ifile.readline()
    for line in ifile:
        line = line.split()
        t.append(float(line[0]))
        amp.append([])
        for i in range(1, len(line)):
            amp[-1].append(float(line[i]))
    ifile.close()

    t = np.array(t, dtype = np.int32)
    amp = np.array(amp, dtype = np.float64)

    return t, amp

##################################################################

def zero_passing_time(t,hyac):
    """ find where the autocorrelation function passes through zero
        for the first time"""
    n = len(hyac)
    idx = 0
    for i in range(n-1):
        if hyac[i]*hyac[i+1] < 0:
            idx = i
            break
    # find zero passing time from linear interpolation
    tm = t[idx]
    tp = t[idx+1]
    hm = hyac[idx]
    hp = hyac[idx+1]
    
    t0 = tm + (tp - tm) / (hp - hm) * (0.0 - hm)
    return t0

##################################################################

def find_extrema(t,hyac,std):
    """ find where the autocorrelation function passes through zero
        for the first time"""
    n = len(hyac)
    te = []
    he = []
    stde = []
    for i in range(1, n-1):
        if (hyac[i-1] - hyac[i]) * ( hyac[i] - hyac[i+1])  < 0:
            te.append(t[i])
            he.append(hyac[i])
            stde.append(std[i])
    te = np.array(te, dtype = np.int32)
    he = np.array(he, dtype = np.float64)
    stde = np.array(stde, dtype =  np.float64)
    # loop over all extrema and check whether they are valid
    nextr = len(te)
    valid = 0
    for i in range(nextr):
        hi = math.fabs(he[i])
        stdi = stde[i]
        if hi - 2*stdi < 0:
            break
        valid = valid + 1
    return te, he, stde, valid

##################################################################

def cart2pol(x,y):
    """ transform cartesian coordinates to polar coordinates"""
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y,x)
    return rho,  phi

##################################################################

def kurtosis2d(x,y):
    """ compute the 2D kurtosis of data x and y"""
    # combine to a single data set
    data = np.array([x,y])
    # average
    av = np.average(data, axis = 1)
    # compute the covariance matrix
    S = np.cov(data)
    # compute the inverse of S
    Sm = linalg.inv(S)


    # compute the kurtosis
    n = len(x)
    b2 = 0
    for i in range(n):
        d = data[:,i] - av
        b2 += np.dot(d, np.dot(Sm, d))**2
    b2 /= n
  
    return b2

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
    ### read in the eigenvalues
    ev = read_eigenvalues(evfile)
    ### read in the autocorrelation of the amplitudes
    tc, ac1, ac1_std, ac2, ac2_std, a1a2, a1a2_std, a2a1, a2a1_std = read_correlations(acfile)
    ### read in the amplitutes
    t, amp = read_amplitudes(ampfile)

    ### extract desired quantities:
    # importance of the first two eigenvalues
    evsum = (ev[0] + ev[1])
    evsum_normed = (ev[0] + ev[1]) / np.sum(ev)
    

    # correlation at t = 0
    ac1_0 = ac1[0]
    ac2_0 = ac2[0]
    a1a2_0 = a1a2[0]
    a2a1_0 = a2a1[0]
    # zero passing time
    tac1_0 = zero_passing_time(tc,ac1)
    tac2_0 = zero_passing_time(tc,ac2)
    ta1a2_0 = zero_passing_time(tc[1:],a1a2[1:]) # don't use first value
    ta2a1_0 = zero_passing_time(tc[1:],a2a1[1:])
    # gradient at t = 0
    dac1_0 = ac1[1] - ac1[0]
    dac2_0 = ac2[1] - ac2[0]
    da1a2_0 = a1a2[1] - a1a2[0]
    da2a1_0 = a2a1[1] - a2a1[0]
    # find significant extrema
    tac1e, ac1e, ac1e_std, valid_ac1 = find_extrema(tc, ac1, ac1_std)
    tac2e, ac2e, ac2e_std, valid_ac2 = find_extrema(tc, ac2, ac2_std)
    ta1a2e, a1a2e, a1a2e_std, valid_a1a2 = find_extrema(tc, a1a2, a1a2_std)
    ta2a1e, a2a1e, a2a1e_std, valid_a2a1 = find_extrema(tc, a2a1, a2a1_std)

 
    ### compute the multivariate curtosis ov amp1 and amp2
    b2 = kurtosis2d(amp[:,0], amp[:,1])

    ### cross-correlations with a1**2 + a2**2
    namp = len(amp[0])
    c = np.zeros((namp), dtype = np.float64)
    p = np.zeros((namp), dtype = np.float64)
    for i in range(namp):
        c[i], p[i] = stats.pearsonr(amp[:,0]**2 + amp[:,1]**2, amp[:,i]**2)

    fig = plt.figure()
    ax1 = plt.subplot(121)
    ax2 = plt.subplot(122)

    ax1.plot(c, ls = '', marker = 'o')
    ax2.plot(p, ls = '', marker = 'o')

    plt.savefig('amplitude_correlations.png')
    plt.close()
    
    ### show results
    fig = plt.figure(figsize = (16,12))
    ax1 = plt.subplot(411)
    ax2 = plt.subplot(412)
    ax3 = plt.subplot(413)
    ax4 = plt.subplot(414)
    xmin = 0
    xmax = tac1e[2]*2
    ymin = -ac1_0
    ymax = ac1_0

    # ac1
    ax1.plot(tc,0*tc, color = '0.5')
    ax1.errorbar(tc,ac1,2*ac1_std)
    ax1.plot([tac1_0, tac1_0], [ymin, ymax])
    ax1.plot(tac1e[:3], ac1e[:3], ls = '', marker = 's')
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])

    # ac2
    ax2.plot(tc,0*tc, color = '0.5')
    ax2.errorbar(tc,ac2,2*ac2_std)
    ax2.plot([tac2_0, tac2_0], [ymin, ymax])
    ax2.plot(tac2e[:3], ac2e[:3], ls = '', marker = 's')
    ax2.set_xlim([xmin, xmax])
    ax2.set_ylim([ymin, ymax])

    # a1a2
    ax3.plot(tc,0*tc, color = '0.5')
    ax3.errorbar(tc,a1a2,2*a1a2_std)
    ax3.plot([ta1a2_0, ta1a2_0], [ymin, ymax])
    ax3.plot(ta1a2e[:3], a1a2e[:3], ls = '', marker = 's')
    ax3.set_xlim([xmin, xmax])
    ax3.set_ylim([ymin, ymax])

    # a2a1
    ax4.plot(tc,0*tc, color = '0.5')
    ax4.errorbar(tc,a2a1,2*a2a1_std)
    ax4.plot([ta2a1_0, ta2a1_0], [ymin, ymax])
    ax4.plot(ta2a1e[:3], a2a1e[:3], ls = '', marker = 's')
    ax4.set_xlim([xmin, xmax])
    ax4.set_ylim([ymin, ymax])

    # safe plot
    plt.savefig('correlation_funcs.png')
    plt.close()



    ### generate a histogram of amp1 and amp2
    multivariate_histogram(amp[:,0], amp[:,1], 'hist_a1a2.png')

    ### generate a plot of evsum / ev
    fig = plt.figure()
    plt.plot(ev / np.sum(ev), ls = '', marker = 'o')
    plt.savefig('eigenvalues.png')
    plt.close()
    

    ### write results to file
    ofile = open('results2.data', 'w')
    ofile.write('Selected results from the pca analysis\n\n')
    # importance of the first two eigenvalues
    ofile.write('ev2 = ' + str(evsum_normed) + '\t' + str(evsum) + '\n')
    # number of valid extrema
    ofile.write('valid_extrema = ' + str(valid_ac1) + '\t' + str(valid_ac2) + '\t' + str(valid_a1a2) + '\t' + str(valid_a2a1) + '\n')
    # kurtosis values
    ofile.write('2dhist_kurtosis = ' + str(b2) + '\n')
    # results for ac1
    ofile.write('ac1(t=0) = ' + str(ac1_0) + '\t+/-' + str(ac1_std[0]) + '\n')
    ofile.write('t(ac1=0) = ' + str(tac1_0) + '\n')
    ofile.write('dac1(t=0) = ' + str(dac1_0) + '\n')
    for i in range (3):
        ofile.write('extremum ' + str(i) + ' (t,ac1,std) = ' + str(tac1e[i]) + ' ' + str(ac1e[i]) + ' ' + str(ac1e_std[i]) + '\n')
    ofile.write('t(last_valid_extremum) = ')
    if valid_ac1 > 0:
        ofile.write(str(tac1e[valid_ac1-1]) + '\n')
    else:
        ofile.write('0\n')
    ofile.write('\n')

    # results for ac2
    ofile.write('ac2(t=0) = ' + str(ac2_0) + '\t+/-' + str(ac2_std[0]) + '\n')
    ofile.write('t(ac2=0) = ' + str(tac2_0) + '\n')
    ofile.write('dac2(t=0) = ' + str(dac2_0) + '\n')
    for i in range (3):
        ofile.write('extremum ' + str(i) + ' (t,ac2,std) = ' + str(tac2e[i]) + ' ' + str(ac2e[i]) + ' ' + str(ac2e_std[i]) + '\n')
    ofile.write('t(last_valid_extremum) = ')
    if valid_ac2 > 0:
        ofile.write(str(tac2e[valid_ac2-1]) + '\n')
    else:
        ofile.write('0\n')
    ofile.write('\n')

    # results for a1a2
    ofile.write('a1a2(t=0) = ' + str(a1a2_0) + '\t+/-' + str(a1a2_std[0]) + '\n')
    ofile.write('t(a1a2=0) = ' + str(ta1a2_0) + '\n')
    ofile.write('da1a2(t=0) = ' + str(da1a2_0) + '\n')
    for i in range (3):
        ofile.write('extremum ' + str(i) + ' (t,a1a2,std) = ' + str(ta1a2e[i]) + ' ' + str(a1a2e[i]) + ' ' + str(a1a2e_std[i]) + '\n')
    ofile.write('t(last_valid_extremum) = ')
    if valid_a1a2 > 0:
        ofile.write(str(ta1a2e[valid_a1a2-1]) + '\n')
    else:
        ofile.write('0\n')
    ofile.write('\n')
    
    # results for a2a1
    ofile.write('a2a1(t=0) = ' + str(a2a1_0) + '\t+/-' + str(a2a1_std[0]) + '\n')
    ofile.write('t(a2a1=0) = ' + str(ta2a1_0) + '\n')
    ofile.write('da2a1(t=0) = ' + str(da2a1_0) + '\n')
    for i in range (3):
        ofile.write('extremum ' + str(i) + ' (t,a2a1,std) = ' + str(ta2a1e[i]) + ' ' + str(a2a1e[i]) + ' ' + str(a2a1e_std[i]) + '\n')
    ofile.write('t(last_valid_extremum) = ')
    if valid_a2a1 > 0:
        ofile.write(str(ta2a1e[valid_a2a1-1]) + '\n')
    else:
        ofile.write('0\n')
    ofile.write('\n')
    
    # close the file
    ofile.close()
    return
    
    
##################################################################

if __name__ == '__main__':
    main()
    
