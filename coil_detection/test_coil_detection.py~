#!/usr/local/bin/python2.7

import numpy as np
from scipy import linalg
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

def check_angles(x,y):
    """ examine whether angles suggest coil formation"""
    nout = 3 # number of atoms at the beginning and end of the
             #  polymer that are not included in the analysis
    natoms = len(x)
    # set initial value for return value
    sangle = 0
    # loop over selected atoms
    for i in range(nout,natoms-nout-2):
        x1 = x[i]
        x2 = x[i+1]
        x3 = x[i+2]
        y1 = y[i]
        y2 = y[i+1]
        y3 = y[i+2]
        sangle = sangle + math.copysign(1,(x2-x1)*(y3-y1) - (y2-y1)*(x3-x1))
        
    return sangle

##################################################################

def run_analysis():
    """ detect existence of coils"""
    # open the file for reading
    ifile = open(fname)
    # loop over all snapshots
    for i in range(nsteps):
        # read data of a single snapshot
        x,y = read_snapshot(ifile)
        # check angles for coil type behavior, if sangle > natoms/2,
        #   the existence of a coil is likely!
        sangle  = check_angles(x,y)
    # close the file
    ifile.close()
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
