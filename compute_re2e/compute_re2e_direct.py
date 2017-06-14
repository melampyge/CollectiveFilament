#!/usr/local/bin/python2.7
### read in xyz atom coordinates
### compute end-to-end distance

import numpy as np
import sys, os, math
import matplotlib.pyplot as plt
from scipy import signal

try:
    fname = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '       filename'
    exit()


#####################################################################
### define / read in some global variables
gamma = 2.0  # viscosity
kT = 1.0     # thermal energy

ifile = open('scale_params.data')
line = ifile.readline()
line = line.split()
L = float(line[-1])    # polymer length
line = ifile.readline()
line = line.split()
dt = float(line[-1])     # simulation timestep


############################################################################

def read_data():
    """ read in the dx and dy differences"""
    ifile = open(fname, 'r')
    ifile.readline()
    ifile.readline()
    dx = []
    dy = []
    for line in ifile:
        line = line.split()
        dx.append(float(line[1]))
        dy.append(float(line[2]))
    ifile.close()
    dx = np.array(dx)
    dy = np.array(dy)
    dx /= L
    dy /= L
    return dx, dy

############################################################################

def main():
    """ compute the average value of the end-to-end distance"""
    # read in the data
    dx,dy = read_data()
    # compute average end-to-end distance
    e2esq = dx**2 + dy**2
    rav = np.sqrt(np.average(e2esq))
    rstd = np.sqrt(np.std(e2esq)/np.sqrt(len(e2esq)))
    # write data to file
    ofile = open('re2e_av.data', 'w')
    ofile.write('Average end-to-end distance and standard deviation of the mean\n\n')
    ofile.write('r_av = ' + str(rav) + '\n')
    ofile.write('r_std = ' + str(rstd) + '\n')
    ofile.close()
    return

############################################################################

if __name__ == '__main__':
    main()
