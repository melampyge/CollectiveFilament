#!/usr/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import os

try:
    infilename = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '       input file'
    exit()

##################################################################

def read_data(ifile):
    t = []
    hyac = []
    std = []
    ifile = open(ifile)
    for i in range(3):
        ifile.readline()
    for line in ifile:
        line = line.split()
        t.append(int(line[0]))
        hyac.append(float(line[1]))
        std.append(float(line[2]))
    ifile.close()
    t = np.array(t, dtype = np.int32)
    hyac = np.array(hyac, dtype = np.float64)
    std = np.array(std, dtype = np.float64)
    return t, hyac, std

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
    return te, he, stde
    
##################################################################
    
def main():
    """ main function, called when the script is started"""
    ### read in the autocorrelation function
    t, hyac, std = read_data(infilename)
    ### extract desired quantities
    # hyac[0]
    h_0 = hyac[0]
    # zero passing time
    t_0 = zero_passing_time(t,hyac)
    # gradient at t = 0
    dh_0 = hyac[1] - hyac[0]
    # first three extrema
    te, he, stde = find_extrema(t, hyac, std)

    ### show results in one plot
    plt.plot(t, 0*t, color = '0.5')
    plt.errorbar(t, hyac, 2*std)
    plt.plot([t_0,t_0], [-h_0, h_0])
    plt.plot(te[:3], he[:3], ls = '', marker = 's')
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,2*te[2],y1,y2))
    plt.ion()
    plt.show()
    ### have values been determined correctly
    valid = int(raw_input('number of valid extrema? ( 0-3)  '))
    plt.close()

    ### write results to file
    ofile = open('results2.data', 'w')
    ofile.write('Selected results from the hy autocorrelation function\n\n')
    ofile.write('h(t=0) = ' + str(h_0) + '\t+/-' + str(std) + '\n')
    ofile.write('t(h=0) = ' + str(t_0) + '\n')
    ofile.write('dh(t=0) = ' + str(dh_0) + '\n')
    ofile.write('valid = ' + str(valid) + '\n')
    for i in range(3):
        ofile.write('extremum ' + str(i) + ' (t,h,std) = ' + str(te[i]) + ' ' + str(he[i]) + ' ' + str(stde[i]) + '\n')
    ofile.close()
    return
    
    
##################################################################

if __name__ == '__main__':
    main()
    
