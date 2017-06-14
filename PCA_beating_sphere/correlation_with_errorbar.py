#!/usr/local/bin/python2.7

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import os
import performance_toolsWrapper

try:
    infilename = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '       parameter file'
    exit()

##################################################################

def read_settings():
    """ read in the settings from the parameter file"""
    # open file for reading
    ifile = open(infilename, 'r')
    # skip comment line
    ifile.readline()
    # filename
    line = ifile.readline()
    line = line.split()
    inputfile = line[-1]
    # compute correlation for every this many snapshots
    line = ifile.readline()
    line = line.split()
    dt = int(line[-1])
    # timestep up to which to compute the correlation (as fraction of all steps)
    line = ifile.readline()
    line = line.split()
    limit = float(line[-1])
    # close file
    ifile.close()
    # return input values
    return inputfile, dt, limit

##################################################################

def read_amplitudes(inputfile):
    """ read the y-position of the last bead"""
    ### open files for reading
    ifile = open(inputfile, 'r')
    ### generate lists for hy and time
    a1 = []
    a2 = []
    # time
    time = []
    ### skip header
    for i in range(3):
        ifile.readline()
    ### read the rest of the file
    for line in ifile:
        line = line.split()
        time.append(int(line[0]))
        a1.append(float(line[1]))
        a2.append(float(line[2]))
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    a1 = np.array(a1, dtype = np.float64)
    a2 = np.array(a2, dtype = np.float64)
    ### return data
    return time, a1, a2

##################################################################

def cross_correlate(time, x1, x2, dt, limit):
    """ compute the autocorrelation function of x"""
    ### compute required values
    nsteps = len(x1)
    limit = int(limit*nsteps)
    ncc = int(limit/dt)
    ### allocate output array to store the results
    t_cc = np.zeros((ncc), dtype = np.int32)
    cc = np.zeros((ncc), dtype = np.float64)
    linval = np.zeros((ncc), dtype = np.int32)
    std = np.zeros((ncc), dtype = np.float64)
    blockdata = np.zeros((nsteps), dtype = np.float64)
    # size of the array for the blocking method
    nblock = int(math.log(nsteps, 2)) + 2
    block_std = np.zeros((nblock), dtype = np.float64)
    block_uncert = np.zeros((nblock), dtype = np.float64)
    ### loop over all molecules, submit MSD computation, read and average
    ###   the results
    performance_tools = performance_toolsWrapper.Performance_tools()
    # submit MSR computation
    performance_tools.compute_crosscorrelation_with_errorbars(t_cc, time, x1, x2, cc, std, blockdata, block_std, block_uncert, linval, nsteps, ncc, limit, dt)
    ### return time and msr
    return t_cc, cc, std

##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    print 'reading settings'
    inputfile, dt, limit = read_settings()
    ### loop over the entire trajectory and read hy
    print 'reading y-position of the last bead'
    time, amp_pca1, amp_pca2 = read_amplitudes(inputfile)

    
    ### compute the crosscorrelation functions
    print 'computing the autocorrlation function'
    t_a1ac, a1ac, std_a1ac = cross_correlate(time, amp_pca1, amp_pca1, dt, limit)
    t_a2ac, a2ac, std_a2ac = cross_correlate(time, amp_pca2, amp_pca2, dt, limit)
    t_a1a2cc, a1a2cc, std_a1a2cc = cross_correlate(time, amp_pca1, amp_pca2, dt, limit)
    t_a2a1cc, a2a1cc, std_a2a1cc = cross_correlate(time, amp_pca2, amp_pca1, dt, limit)
         
    #### write tables
    # correlation functions
    ofile = open('amplitude_correlations.data', 'w')
    ofile.write('Auto- and cross-correlation functions of the amplitudes\n\n')
    ofile.write('lag time\ta1ac\tstd\ta2ac\tstd\ta1a2cc\tstd\ta2a1cc\tstd\n')
    nsteps = len(t_a1ac)
    for i in range(nsteps):
        ofile.write(str(t_a1ac[i]) + '\t' + str(a1ac[i]) + '\t' + str(std_a1ac[i]) + '\t' + str(a2ac[i]) + '\t' + str(std_a2ac[i]) + '\t' + str(a1a2cc[i]) + '\t' + str(std_a1a2cc[i]) + '\t' + str(a2a1cc[i]) + '\t' + str(std_a2a1cc[i]) + '\n')
    ofile.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
