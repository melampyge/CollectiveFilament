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
    # timestep up to which to compute the msr (as fraction of all steps)
    line = ifile.readline()
    line = line.split()
    limit = float(line[-1])
    # output folder
    line = ifile.readline()
    line = line.split()
    ofname = line[-1]
    # close file
    ifile.close()
    # return input values
    return inputfile, dt, limit, ofname

##################################################################

def read_hy(inputfile):
    """ read the y-position of the last bead"""
    ### open files for reading
    ifile = open(inputfile, 'r')
    ### generate lists for hy and time
    hy = []
    # time
    time = []
    ### skip header
    for i in range(3):
        ifile.readline()
    ### read the rest of the file
    for line in ifile:
        line = line.split()
        time.append(int(line[0]))
        hy.append(float(line[1]))
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    hy = np.array(hy, dtype = np.float64)
    ### return data
    return time, hy

##################################################################

def compute_autocorrelation(time, x, dt, limit):
    """ compute the autocorrelation function of x"""
    ### compute required values
    nsteps = len(x)
    limit = int(limit*nsteps)
    nac = int(limit/dt)
    ### allocate output array to store the results
    t_ac = np.zeros((nac), dtype = np.int32)
    ac = np.zeros((nac), dtype = np.float64)
    linval = np.zeros((nac), dtype = np.int32)
    std = np.zeros((nac), dtype = np.float64)
    blockdata = np.zeros((nsteps), dtype = np.float64)
    # size of the array for the blocking method
    nblock = int(math.log(nsteps, 2)) + 2
    block_std = np.zeros((nblock), dtype = np.float64)
    block_uncert = np.zeros((nblock), dtype = np.float64)
    ### loop over all molecules, submit MSD computation, read and average
    ###   the results
    performance_tools = performance_toolsWrapper.Performance_tools()
    # submit MSR computation
    performance_tools.compute_autocorrelation_with_errorbars(t_ac, time, x, ac, std, blockdata, block_std, block_uncert, linval, nsteps, nac, limit, dt)
    ### return time and msr
    return t_ac, ac, std

##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    print 'reading settings'
    inputfile, dt, limit, ofname = read_settings()
    ### loop over the entire trajectory and read hy
    print 'reading y-position of the last bead'
    time, hy = read_hy(inputfile)

    ### compute the autocorrelation function
    print 'computing the autocorrlation function'
    t_hyac, hyac, std = compute_autocorrelation(time, hy, dt, limit)
    
    ### generate output folder
    os.system('mkdir ' + str(ofname))
       
    #### write tables
    # correlation functions
    ofile = open(ofname + '/correlation_functions_with_std.data', 'w')
    ofile.write('Correlation functions\n\n')
    ofile.write('lag time\tHYAC\tstd\n')
    for i in range(len(t_hyac)):
        ofile.write(str(t_hyac[i]) + '\t' + str(hyac[i]) + '\t' + str(std[i]) + '\n')
    ofile.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
