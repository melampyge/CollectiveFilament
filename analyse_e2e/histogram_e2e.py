#!/usr/local/bin/python2.7

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


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
ttrans = gamma*L**3/6./kT

#####################################################################

def read_data():
    """ read in the data"""
    ifile = open('e2e.data')
    ifile.readline()
    ifile.readline()
    r = []
    for line in ifile:
        line = line.split()
        dx = float(line[1])
        dy = float(line[2])
        d = dx**2 + dy**2
        r.append(d)
    ifile.close()
    r = np.array(r)
    r = np.sqrt(r)
    r /= L
    return r

#####################################################################

def main():
    """ main function"""
    # read in the e2e vector
    e2e = read_data()
    # generate a histogram
    nbins = 500
    hist, edges = np.histogram(e2e, bins = nbins, range = (0,1), density = True)
    # create plot of the histogram
    plt.hist(e2e, bins = nbins, range = (0,1))
    plt.savefig('hist_e2e.png')
    plt.close()
    # write results to file
    ofile = open('hist_e2e.data', 'w')
    ofile.write('Histogram of the end-to-end distance\n\n')
    ofile.write('lower_edge\tupper_edge\tp\n')
    for i in range(len(hist)):
        ofile.write(str(edges[i]) + '\t' + str(edges[i+1]) + '\t' + str(hist[i]) + '\n')
    ofile.close()
    return

#####################################################################

if __name__ == '__main__':
    main()
