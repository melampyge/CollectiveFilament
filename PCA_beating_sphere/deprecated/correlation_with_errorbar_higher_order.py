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

def read_amplitudes(inputfile):
    """ read the y-position of the last bead"""
    ### open files for reading
    ifile = open(inputfile, 'r')
    ### generate lists for hy and time
    a1 = []
    a2 = []
    a3 = []
    a4 = []
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
        a3.append(float(line[3]))
        a4.append(float(line[4]))
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    a1 = np.array(a1, dtype = np.float64)
    a2 = np.array(a2, dtype = np.float64)
    a3 = np.array(a3, dtype = np.float64)
    a4 = np.array(a4, dtype = np.float64)
    ### return data
    return time, a1, a2, a3, a4

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
    inputfile, dt, limit, ofname = read_settings()
    ### loop over the entire trajectory and read hy
    print 'reading y-position of the last bead'
    time, amp1, amp2, amp3, amp4 = read_amplitudes(inputfile)

    
    ### compute the crosscorrelation functions
    print 'computing the autocorrlation function'
    t, a1, sa1 = cross_correlate(time, amp1, amp1, dt, limit)
    t, a2, sa2 = cross_correlate(time, amp2, amp2, dt, limit)
    t, a3, sa3 = cross_correlate(time, amp3, amp3, dt, limit)
    t, a4, sa4 = cross_correlate(time, amp4, amp4, dt, limit)

    t, c12, sc12 = cross_correlate(time, amp1, amp2, dt, limit)
    t, c13, sc13 = cross_correlate(time, amp1, amp3, dt, limit)
    t, c14, sc14 = cross_correlate(time, amp1, amp4, dt, limit)
    t, c23, sc23 = cross_correlate(time, amp2, amp3, dt, limit)
    t, c24, sc24 = cross_correlate(time, amp2, amp4, dt, limit)
    t, c34, sc34 = cross_correlate(time, amp3, amp4, dt, limit)

  
    ### generate output folder
    os.system('mkdir ' + str(ofname))
       
    #### write tables
    # autocorrelation functions
    ofile = open(ofname + '/autocorrelations_with_errorbars.data', 'w')
    ofile.write('Autocorrelation functions of the amplitudes\n\n')
    ofile.write('lag time\ta1\tstd\ta2\tstd\ta3\tstd\ta4\tstd\n')
    nsteps = len(t)
    for i in range(nsteps):
        ofile.write(str(t[i]) + '\t' + str(a1[i]) + '\t' + str(sa1[i]) + '\t' + str(a2[i]) + '\t' + str(sa2[i]) + '\t' + str(a3[i]) + '\t' + str(sa3[i]) + '\t' + str(a4[i]) + '\t' + str(sa4[i]) + '\n')
    ofile.close()

    # crosscorrelation functions
    ofile = open(ofname + '/crosscorrelations_with_errorbars.data', 'w')
    ofile.write('Crosscorrelation functions of the amplitudes\n\n')
    ofile.write('lag time\tc12\tstd\tc13\tstd\tc14\tstd\tc23\tstd\tc24\tstd\tc34\t\std\n')
    nsteps = len(t)
    for i in range(nsteps):
        ofile.write(str(t[i]) + '\t' + str(c12[i]) + '\t' + str(sc12[i]) + '\t' + str(c13[i]) + '\t' + str(sc13[i]) + '\t' + str(c14[i]) + '\t' + str(sc14[i]) + '\t' + str(c23[i]) + '\t' + str(sc23[i]) + '\t' + str(c24[i]) + '\t' + str(sc24[i]) + '\t' + str(c34[i]) + '\t' + str(sc34[i]) + '\n')
    ofile.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
