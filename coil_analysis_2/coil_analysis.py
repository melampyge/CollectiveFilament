#!/usr/local/bin/python2.7

import numpy as np
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
    # read in the timestep
    line = ifile.readline()
    line = line.split()
    t = int(line[-1])
    # read in the coordinates
    for i in range(natoms):
        line = ifile.readline()
        line = line.split()
        x[i] = float(line[-3])
        y[i] = float(line[-2])
    # return the data array
    return x,y,t

##################################################################

def compute_curvatures(x,y):
    """ compute the local curvature at every atom"""
    natoms = len(x)
    curvatures = []
    for i in range(5,natoms-1):
        # compute all local derivatives
        x1 = x[i-1]
        x2 = x[i]
        x3 = x[i+1]
        y1 = y[i-1]
        y2 = y[i]
        y3 = y[i+1]
        dx = (x3-x1)/2
        dy = (y3-y1)/2
        d2x = (x3-2*x2+x1)
        d2y = (y3-2*y2+y1)
        ci = (dx*d2y - dy*d2x)/(np.sqrt(dx**2 + dy**2)**3)
        curvatures.append(ci)
    curvatures = np.array(curvatures)
    return curvatures

##################################################################

def compute_sangle(x,y):
    """ compute the how much the next point is to the left or right"""
    sangles = []
    natoms = len(x)
    for i in range(natoms-2):
        x1 = x[i]
        x2 = x[i+1]
        x3 = x[i+2]
        y1 = y[i]
        y2 = y[i+1]
        y3 = y[i+2]
        sangles.append((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1))
    return sangles

##################################################################

def tangential_correlation(x,y):
    """ compute the tangential correlatin of the beads"""
    tcorr = []
    natoms = len(x)
    n0 = 5
    delxc = x[n0] - x[n0+1]
    delyc = y[n0] - y[n0+1]
    for i in range(n0,natoms-1):
        delx = x[i] - x[i+1]
        dely = y[i] - y[i+1]
        tcorr.append(delx*delxc + dely*delyc)
    tcorr = np.array(tcorr)
    return tcorr

##################################################################

def std_on_the_fly(i,x,m,s):
    """ compute standard deviation and average on the fly"""
    d = x - m
    m += d/(i+1)
    s += d*(x-m)
    return

##################################################################

def run_analysis():
    """ detect whether coils are existent and how compact they are"""
    # open the file for reading
    ifile = open(fname)
    cav = np.zeros((nsteps))
    sav = np.zeros((nsteps))
    # loop over all snapshots
    for i in range(nsteps):
        # describe Progress
        prog = float(i+1)/float(nsteps)*100
        sys.stdout.write("\r    Progress: %1.2f%%" % prog)
        sys.stdout.flush()
        if i == nsteps-1:
            sys.stdout.write('\n')
        # read data of a single snapshot
        x,y,ti = read_snapshot(ifile)
        # check the distance of the first atom to all other atoms
        #curvatures = compute_curvatures(x,y)
        #sangles = compute_sangle(x,y)
        # compute the average and store it
        #sav[i] = np.average(sangles)
        # use correlation between tangential vectors
        tcorri = tangential_correlation(x,y)
        if i == 0:
            tcorr = np.zeros((len(tcorri)))
            tcorr_std = np.zeros((len(tcorri)))
        std_on_the_fly(i,tcorri,tcorr,tcorr_std)
    tcorr_std /= (nsteps-1)
    tcorr_std = np.sqrt(tcorr_std/nsteps)

    x = np.linspace(1,len(tcorr),len(tcorr))
    plt.errorbar(x,tcorr,tcorr_std)
    plt.show()
    plt.close()
    # close the input file
    ifile.close()
    return

##################################################################

def main():
    """ main function, called when script is started"""
    # run the analysis
    run_analysis()
    return

##################################################################

if __name__ == '__main__':
    main()
