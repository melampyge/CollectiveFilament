#!/usr/local/bin/python2.7

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import sys, math

try:
    fname = sys.argv[1]
    nsteps = int(sys.argv[2])
except:
    print 'Usage: ' + sys.argv[0] + '      infilename         #(steps)'
    exit()
    
##################################################################

def read_snapshot(ifile):
    """ read in the coordinates from the xyz file"""
    # get the number of atoms
    line = ifile.readline()
    line = line.split()
    natoms = int(line[0])
    # allocate array to store data
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    # skip one line
    ifile.readline()
    # read in the coordinates
    for i in range(natoms):
        line = ifile.readline()
        line = line.split()
        x[i] = float(line[-3])
        y[i] = float(line[-2])
    # return the data array
    return x,y

##################################################################

def compute_rotation(x,y):
    """ compute the tangent correlation"""
    natom = len(x)
    # compute angle with respect to some arbitrary orientation
    phi = []
    for i in range(1,natom):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        dphi = math.atan2(dy,dx)
        phi.append(dphi)
    phi = np.array(phi)
    # correct for 2*pi periodicity
    

so
    return phi2

##################################################################

def run_analysis():
    """ detect existence of coils"""
    # open the file for reading
    ifile = open(fname)
    #s = []
    #s2 = []
    # loop over all snapshots
    for j in range(nsteps):
        if j % 1000 == 0:
            print j
        # read data of a single snapshot
        x,y = read_snapshot(ifile)
        # compute angular correlation function
        phi = compute_rotation(x,y)
        # make a linear regression to the data
        bondnumber = np.linspace(0,len(phi)-1,len(phi))
        #slope, intercept, r_value, p_value, std_err = stats.linregress(bondnumber,phi)
        # compute the sum of the squares between the regression line and phi
        #RMS = np.sqrt(np.average((bondnumber*slope+intercept - phi)**2))
        ### compute the slope and the intercept using an explicit algorithm

        t1 = 0
        t2 = 0
        t3 = 0
        t4 = 0
        n = len(phi)
        for i in range(n):
            xi = float(i+1)
            yi = phi[i]
            t1 = t1 + xi*yi
            t2 = t2 + xi
            t3 = t3 + yi
            t4 = t4 + xi*xi
        slope = (t1 - t2*t3/n) / (t4 - t2*t2/n)
        intercept = t3/n - slope*t2/n
        slope2 = np.average(phi[1:] - phi[:-1])
        slope3 = (phi[-1] - phi[0])/98.
        fig = plt.figure()
        ax1 = plt.subplot(211)
        ax2 = plt.subplot(212)
        ax1.plot(x,y)
        ax1.set_aspect('equal')
        ax1.set_xlabel(r'$x$ [$\sigma$]')
        ax1.set_ylabel(r'$y$ [$\sigma$]')
        ax2.plot(phi)
        ax2.plot(bondnumber*slope + intercept, label = 'coilicity = ' + str(slope))
        legend = plt.legend()
        ax2.set_xlabel(r'contour length [$\sigma]$]')
        ax2.set_ylabel(r'$\theta$')
        plt.savefig('fig_' + str(j) + '.png')
        plt.show()
        plt.close()
        #s.append(math.fabs(slope))
    # close the file
    ifile.close()
    #s = np.array(s)
    #plt.hist(s,bins=int(np.sqrt(len(s))))
    #plt.show()
    #plt.close()
    return

##################################################################

def main():
    """ main function, called when script is started"""
    # compute the power spectrum
    run_analysis()
    return

##################################################################

if __name__ == '__main__':
    main()
