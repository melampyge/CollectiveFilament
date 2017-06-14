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

def get_com(inputfile):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = open(inputfile, 'r')
    for i in range(3):
        ifile.readline()
    ### generate lists for ex and ey
    comx = []
    comy = []
    # time
    time = []
    ### read rest of the file
    for line in ifile:
        line = line.split()
        time.append(int(line[0]))
        comx.append(float(line[1]))
        comy.append(float(line[2]))
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    comx = np.array(comx, dtype = np.float64)
    comy = np.array(comy, dtype = np.float64)
    ### return data
    return time, comx, comy

##################################################################

def compute_msd_with_errorbars(time, comx, comy, dt, limit):
    """ compute and average the MSR"""
    ### compute required values
    nsteps = len(comx)
    limit = int(limit*nsteps)
    nmsd = int(limit/dt)
    ### allocate output array to store the results
    t_msd = np.zeros((nmsd), dtype = np.int32)
    msd = np.zeros((nmsd), dtype = np.float64)
    std = np.zeros((nmsd), dtype = np.float64)
    linval = np.zeros((nmsd), dtype = np.int32)
    ### values for the blocking method
    blockdata = np.zeros((nsteps), dtype = np.float64)
    nblock = int(math.log(nsteps, 2)) + 2
    block_std = np.zeros((nblock), dtype = np.float64)
    block_uncert = np.zeros((nblock), dtype = np.float64)
    ### loop over all molecules, submit MSD computation, read and average
    ###   the results
    performance_tools = performance_toolsWrapper.Performance_tools()
    # submit MSR computation
    performance_tools.compute_MSD_with_errorbars(t_msd, time, comx, comy, msd, std, blockdata, block_std, block_uncert, linval, nsteps, nmsd, limit, dt)
    ### return time and msr
    return t_msd, msd, std

##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    print 'reading settings'
    inputfile, dt, limit, ofname = read_settings()
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    print 'reading end-to-end vector'
    time, comx, comy = get_com(inputfile)
    
    ### compute msd
    print 'computing mean square displacement'
    t_msd, msd, std = compute_msd_with_errorbars(time, comx, comy, dt, limit)
    
    ### generate output folder
    os.system('mkdir ' + str(ofname))
       
    #### write tables
    # correlation functions
    ofile = open(ofname + '/correlation_functions_with_errorbars.data', 'w')
    ofile.write('Correlation functions\n\n')
    ofile.write('lag time\tMSD\tstd\n')
    for i in range(len(t_msd)):
        ofile.write(str(t_msd[i]) + '\t' + str(msd[i]) + '\t' + str(std[i]) + '\n')
    ofile.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
