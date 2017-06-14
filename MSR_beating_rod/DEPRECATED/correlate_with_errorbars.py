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
    # compute msr for every this many snapshots
    line = ifile.readline()
    line = line.split()
    dt = int(line[-1])
    # timestep up to which to compute the msd
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

def get_angle(inputfile):
    """ compute phi"""
    ifile = open(inputfile)
    for i in range(3):
        ifile.readline()
    time = []
    phi = []
    for line in ifile:
        line = line.split()
        time.append(int(line[0]))
        phi.append(float(line[1]))
    time = np.array(time, dtype = np.int32)
    phi = np.array(phi, dtype = np.float64)
    return time, phi

##################################################################

def compute_msr(time, phi, dt, limit):
    """ compute and average the MSR"""
    ### compute required values
    nsteps = len(phi)
    limit = int(limit*nsteps)
    nmsr = int(limit/dt)
    ### allocate output array to store the results
    t_msr = np.zeros((nmsr), dtype = np.int32)
    msr = np.zeros((nmsr), dtype = np.float64)
    std = np.zeros((nmsr), dtype = np.float64)
    linval = np.zeros((nmsr), dtype = np.int32)
    ### values for the blocking method
    blockdata = np.zeros((nsteps), dtype = np.float64)
    nblock = int(math.log(nsteps, 2)) + 2
    block_std = np.zeros((nblock), dtype = np.float64)
    block_uncert = np.zeros((nblock), dtype = np.float64)
    # submit MSR computation
    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.compute_MSR_with_errorbars(t_msr, time, phi, msr, std, blockdata, block_std, block_uncert, linval, nsteps, nmsr, limit, dt)
    ### return time and msr and std
    return t_msr, msr, std


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

def compute_derivatives(x):
    """ compute the first and second derivative using the central limit theorem"""
    n = len(x)
    dx = np.zeros((n), dtype = np.float64)
    d2x = np.zeros((n), dtype = np.float64)
    for i in range(1,n-1):
        dx[i] = x[i+1] - x[i-1]
        d2x[i] = x[i+1] + x[i-1] - 2*x[i]
    dx /= 2
    d2x /= 2
    # border values
    d2x[0] = d2x[1]
    d2x[-1] = d2x[-2]
    dx[0] = x[1] - x[0]
    dx[-1] = x[-1] - x[-2]
    return dx, d2x


##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    print 'reading settings'
    inputfile, dt, limit, ofname = read_settings()
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    print 'reading data'
    time, phi = get_angle(inputfile)
    # compute derivative of phi
    dphi, d2phi = compute_derivatives(phi)
    
    ### compute msr
    print 'computing mean square rotation'
    t_msr, msr, std_msr = compute_msr(time, phi, dt, limit)
    print 'computing the rotational velocity autocorrelation function'
    t_vac, vac, std_vac = compute_autocorrelation(time, dphi, dt, limit)
    
    ### generate output folder
    os.system('mkdir ' + str(ofname))
       
    #### write tables
    # correlation functions
    ofile = open(ofname + '/correlation_functions_with_errorbars.data', 'w')
    ofile.write('Correlation functions\n\n')
    ofile.write('lag time\tMSR\tstd\tVAC\tstd\n')
    for i in range(len(t_msr)):
        ofile.write(str(t_msr[i]) + '\t' + str(msr[i]) + '\t' + str(std_msr[i]) + '\t' + str(vac[i]) + '\t' + str(std_vac[i]) + '\n')
    ofile.close()
    
        
##################################################################

if __name__ == '__main__':
    main()
    
