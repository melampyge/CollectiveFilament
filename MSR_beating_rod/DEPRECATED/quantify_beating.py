#!/usr/local/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import os

try:
    infilename = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '       msr file'
    exit()

##################################################################

def read_msr(infilename):
    """ read in the msr data"""
    t = []
    msr = []
    ifile = open(infilename)
    for i in range(3):
        ifile.readline()
    for line in ifile:
        line = line.split()
        t.append(int(line[0]))
        msr.append(float(line[1]))
    t = np.array(t, dtype = np.int32)
    msr = np.array(msr, dtype = np.float64)

    return t, msr
    
##################################################################

def compute_derivative_1(x, omit = 0):
    """ compute the first derivative based on central limit"""
    n = len(x)
    dx = np.zeros((n), dtype = np.float64)
    for i in range(1+omit,n-1-omit):
        dx[i] = 0.5*(x[i+1] - x[i])
    return dx

##################################################################

def compute_derivative_2(x, omit = 0):
    """ compute the first derivative smoothened by 2 values in each direction"""
    n = len(x)
    dx = np.zeros((n), dtype = np.float64)
    for i in range(2+omit,n-2-omit):
        xm2 = x[i-2]
        xm1 = x[i-1]
        xp1 = x[i+1]
        xp2 = x[i+2]
        dx[i] = 2*(xp1-xm1) + xp2 - xm2
    dx /= 8
    return dx

#################################################################

def compute_derivative_3(x, omit = 0):
    """ compute the first derivative smoothend by 3 values in each direction"""
    n = len(x)
    dx = np.zeros((n), dtype = np.float64)
    for i in range(3+omit,n-3-omit):
        xm3 = x[i-3]
        xm2 = x[i-2]
        xm1 = x[i-1]
        xp1 = x[i+1]
        xp2 = x[i+2]
        xp3 = x[i+3]
        dx[i] = 5*(xp1-xm1) + 4*(xp2-xm2) + xp3-xm3
    dx /= 32
    return dx

##################################################################

def compute_2nd_derivative(x):
    """ compute the second derivative using the central limit theorem"""
    n = len(x)
    d2x = np.zeros((n), dtype = np.float64)
    for i in range(2,n-1):
        xm1 = x[i-1]
        xi = x[i]
        xp1 = x[i+1]
        d2x[i] = xp1 - 2*xi + xm1
    return d2x

##################################################################

def compute_2nd_derivative_SGfilter(x, order = 2):
    """ compute the second derivative using an SGfilter"""
    n = len(x)
    d2x = np.zeros((n), dtype = np.float64)
    if order == 2:
        for i in range(2,n-2):
            xm2 = x[i-2]
            xm1 = x[i-1]
            xi = x[i]
            xp1 = x[i+1]
            xp2 = x[i+2]
            d2x[i] = 2*xm2  - xm1 - 2*xi - xp1 + 2*xp2
        d2x /= 7
    if order == 3:
        for i in range(3,n-3):
            xm3 = x[i-3]
            xm2 = x[i-2]
            xm1 = x[i-1]
            xi = x[i]
            xp1 = x[i+1]
            xp2 = x[i+2]
            xp3 = x[i+3]
            d2x[i] = 5*xm3 -3*xm1 -4*xi -3*xp1 + 5*xp3
        d2x /= 42
    if order == 4:
        for i in range(4,n-4):
            xm4 = x[i-4]
            xm3 = x[i-3]
            xm2 = x[i-2]
            xm1 = x[i-1]
            xi = x[i]
            xp1 = x[i+1]
            xp2 = x[i+2]
            xp3 = x[i+3]
            xp4 = x[i+4]
            d2x[i] = 28*xm4 +7*xm3 -8*xm2 -17*xm1 -20*xi -17*xp1 -8*xp2 +7*xp3 + 28*xp4
        d2x /= 462
    return d2x

##################################################################

def find_extrema(dmsr, omit = 0):
    """ search for extrema by checking where the derivative passes
        throug zero """
    n = len(dmsr)
    idxmin = -n
    idxmax = -n
    for i in range(1+omit,n-2-omit):
        if dmsr[i]*dmsr[i+1] < 0:
            if idxmax < 0:
                idxmax = i
            else:
                idxmin = i
                break
    return idxmin, idxmax

##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    print 'reading MSR'
    t_msr, msr = read_msr(infilename)
    ### compute first and second derivative
    print 'computing derivatives'
    dmsr = compute_derivative_1(msr)
    d2msr = compute_derivative_2(dmsr, 2)
    d2msr = compute_2nd_derivative_SGfilter(msr, order = 4a)

    # compute the position of the first maximum and minimum
    idxmin, idxmax = find_extrema(dmsr) 
    tmin = t_msr[idxmin]
    tmax = t_msr[idxmax]
    msr_min = msr[idxmin]
    msr_max = msr[idxmax]

    # compute the position of the inflection points
    idxt2, idxt1 = find_extrema(d2msr, omit = 3)
    tt1 = t_msr[idxt1]
    tt2 = t_msr[idxt2]
    msr_t1 = msr[idxt1]
    msr_t2 = msr[idxt2]

    # line with slope 2 that touches the second MSD point
    s2_start = msr[1]
    s2_end = max(msr)
    t2_start = t_msr[1]
    t2_end = np.sqrt(s2_end/s2_start)*t2_start

    # line with slope 1 that touches the las MSD point
    t1_end = t_msr[-1]
    s1_start = msr[1]
    t1_start = t_msr[1]
    s1_end = max(msr)
    t1_end = s1_end/s1_start*t1_start

    # compute the maximum slope in log log representation
    log_msr = np.log(msr[1:])
    log_time = np.log(t_msr[1:])
    d_log_msr = (log_msr[2:] - log_msr[:-2])/(log_time[2:] - log_time[:-2])
    #dmax = max(d_log_msr[log_time[1:-1] < 16])
    
    # generate a plot, show tmin and tmax in legend
    fig = plt.figure()
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(223)
    ax4 = plt.subplot(224)

    ax1.loglog(t_msr,msr)#, label = 'tmin = ' + str(tmin) + '\ntmax = ' + str(tmax))
    ax1.loglog([tmin,tmin],[min(msr)+0.0001,max(msr)])
    ax1.loglog([tmax,tmax],[min(msr)+ 0.0001,max(msr)])
    ax1.loglog([t1_start, t1_end],[s1_start,s1_end], color = 'k')
    ax1.loglog([t2_start, t2_end],[s2_start,s2_end], color = 'k')
    #ax1.legend(loc = 'lower right')
    
    ax2.semilogx(t_msr,dmsr)
    ax2.semilogx([tt1,tt1],[min(dmsr),max(dmsr)])
    ax2.semilogx([tt2,tt2],[min(dmsr),max(dmsr)])

    ax3.plot(t_msr, 0*t_msr, color = '0.5')
    ax3.plot(t_msr, d2msr)
    ax3.set_xlim([0,5*tt2])
    ax3.plot([tt1,tt1],[min(d2msr),max(d2msr)])
    ax3.plot([tt2,tt2],[min(d2msr),max(d2msr)])

    ax4.plot(log_time[1:-1], d_log_msr)
    
    plt.ion()
    plt.show()
    extrema_correct = int(raw_input('extrema correct? (0 = no, 1 = yes) '))
    inflection_correct = int(raw_input('inflection correct? (0 = no, 1 = yes) '))
    superdiffusive =  int(raw_input('superdiffusive? (0 = no, 1 = yes) '))
    tcrit = float(raw_input('tmax to determine maximum slope? '))
    dmax = max(d_log_msr[log_time[1:-1] < tcrit])
    plt.savefig('curves.png')
    plt.close()

    ofile = open('msr_extrema.data', 'w')
    ofile.write('First Maximum and minimum\n\n')
    ofile.write('extrema_correct = ' + str(extrema_correct) + '\n')
    ofile.write('inflection_correct = ' + str(inflection_correct) + '\n')
    ofile.write('superdiffusive = ' + str(superdiffusive) + '\n')
    ofile.write('dmax = ' + str(dmax) + '\n')
    ofile.write('t_max, msr_max = ' + str(tmax) + ' ' + str(msr_max) + '\n')
    ofile.write('t_min, msr_min = ' + str(tmin) + ' ' + str(msr_min) + '\n')
    ofile.write('t_t1, msr_t1 = ' + str(tt1) + ' ' + str(msr_t1) + '\n')
    ofile.write('t_t2, msr_t2 = ' + str(tt2) + ' ' + str(msr_t2) + '\n')
    ofile.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
