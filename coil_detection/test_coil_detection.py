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

def check_coil(x,y):
    """ examine whether angles suggest coil formation"""
    natoms = len(x)
    # set initial value for return value
    # loop over selected atoms
    phimin = 0
    phimax = 0
    phi = 0
    coil = 0
    for i in range(natoms-2):
        x1 = x[i]
        x2 = x[i+1]
        x3 = x[i+2]
        y1 = y[i]
        y2 = y[i+1]
        y3 = y[i+2]
        sangle = math.copysign(1,(x2-x1)*(y3-y1) - (y2-y1)*(x3-x1))
        delx1 = x2 - x1
        delx2 = x3 - x2
        dely1 = y2 - y1
        dely2 = y3 - y2
        dcosphi = (delx1*delx2 + dely1*dely2)/np.sqrt(delx1**2 + dely1**2)/np.sqrt(delx2**2 + dely2**2)
        dphi = np.arccos(dcosphi)
        phi += dphi*sangle
        #print dcosphi,dphi,phi
        if phi < phimin:
            phimin = phi
        if phi > phimax:
            phimax = phi
        if phimax - phimin >= 2*np.pi:
            coil = 1
    if coil == 1:
        phimin = 0
        phimax = 0
        phi = 0
        coil = 0
        for i in range(natoms-2):
            x1 = x[i]
            x2 = x[i+1]
            x3 = x[i+2]
            y1 = y[i]
            y2 = y[i+1]
            y3 = y[i+2]
            sangle = math.copysign(1,(x2-x1)*(y3-y1) - (y2-y1)*(x3-x1))
            delx1 = x2 - x1
            delx2 = x3 - x2
            dely1 = y2 - y1
            dely2 = y3 - y2
            dcosphi = (delx1*delx2 + dely1*dely2)/np.sqrt(delx1**2 + dely1**2)/np.sqrt(delx2**2 + dely2**2)
            dphi = np.arccos(dcosphi)
            phi += dphi*sangle
            #print dcosphi,dphi,phi
            if phi < phimin:
                phimin = phi
            if phi > phimax:
                phimax = phi
            print dcosphi,dphi,phi,phimin,phimax
        plt.plot(x,y)
        plt.axis('equal')
        plt.show()
        plt.close()
    return coil

##################################################################

def run_analysis():
    """ detect existence of coils"""
    # open the file for reading
    ifile = open(fname)
    coil = np.zeros((nsteps))
    # loop over all snapshots
    for i in range(nsteps):
        # describe Progress
        prog = float(i+1)/float(nsteps)*100
        sys.stdout.write("\r    Progress: %1.2f%%" % prog)
        sys.stdout.flush()
        if i == nsteps-1:
            sys.stdout.write('\n')
        # read data of a single snapshot
        x,y = read_snapshot(ifile)
        # check angles for coil type behavior, if sangle > natoms/2,
        #   the existence of a coil is likely!
        coil[i] = check_coil(x,y)
    # close the file
    ifile.close()
    # check for coils
    for i in range(nsteps):
        if coil[i] == 1:
            print i
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
