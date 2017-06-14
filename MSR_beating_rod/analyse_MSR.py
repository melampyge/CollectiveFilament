#!/usr/local/bin/python2.7

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

try:
    msrfile = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '      msr file'
    exit()

#####################################################################

def read_data(msrfile):
    """ read in the data"""
    ifile = open(msrfile)
    for i in range(3):
        ifile.readline()
    t_msr = []
    msr = []
    for line in ifile:
        line = line.split()
        t_msr.append(float(line[0]))
        msr.append(float(line[1]))
    
    # transfrom to arrays
    t_msr = np.array(t_msr, dtype = np.float64)
    msr = np.array(msr, dtype = np.float64)
    return t_msr, msr

#####################################################################

def find_dr(t,MSR):
    """ compute the rotational diffusion coefficient"""
    # compute log versions of the data
    tlog = np.log10(t)
    Mlog = np.log10(MSR)
    kmin = MSR[0] / t[0]
    kmax = MSR[-1] / t[-1]
    n = 6
    k = np.logspace(np.log10(kmin), np.log10(kmax),n)
    s = []
    for i in range(n):
        s.append(np.log10(k[i]*t))
    # find the linear regime of the data
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(tlog,Mlog)
    for i in range(n):
        ax.plot(tlog,s[i],lw = 0.5, color = '0.5')

    ax.grid(which = 'major')
    ax.grid(which = 'minor')
    plt.ion()
    plt.show()
    while True:
        try:
            # read in the input and output values
            inpt = raw_input('lower and upper bounds:  ')
            inpt = inpt.split()
            tmin = float(inpt[0])
            tmax = float(inpt[1])
            # compute the slope of the selected regime
            res = stats.linregress(tlog[(tlog >= tmin) & (tlog <= tmax)],Mlog[(tlog >= tmin ) & ( tlog <= tmax) ])
            lslope = res[0]
            inpt = raw_input('slope = ' + str(lslope) + '. Accept? (0 / 1):  ')
            inpt = int(inpt)
            if inpt == 1:
                break
        except:
            pass
    plt.close()
    # perform linear fit
    res  = stats.linregress(t[(tlog >= tmin) & (tlog <= tmax)],MSR[(tlog >= tmin ) & ( tlog <= tmax) ])
    slope = res[0]
    dr = 0.5*slope
    print 'Dr = ' + str(dr)
    return dr, lslope

#####################################################################

def  write_results(t,MSR,dr,lslope):
     ofile = open('diff_coeff.data', 'w')
     ofile.write('Results from analyzing the MSR\n\n')
     ofile.write('dr = ' + str(dr) + '\n')
     ofile.write('lslope = ' + str(lslope) + '\n')
     ofile.close()
     return
#####################################################################

def main():
    """ main function"""
    # read in the MSD data
    t_msr, msr = read_data(msrfile)
    # fit curves
    dr, lslope = find_dr(t_msr[1:],msr[1:])
    # write results to file
    write_results(t_msr,msr,dr,lslope)
    return

#####################################################################

if __name__ == '__main__':
    main()
