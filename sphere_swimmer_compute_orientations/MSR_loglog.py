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
    inputfile = sys.argv[1]
    nmsr = int(sys.argv[2])
    limit = float(sys.argv[3])
except:
    print 'Usage: ' + sys.argv[0] + '       inupt file          #(MSR points)        length MSR / length simulation'
    exit()

##################################################################

def read_angle(inputfile):
    """ read in the angles from the input file"""
    ifile = open('orientations.data')
    t = []
    phi = []
    for i in range(3):
        ifile.readline()
    for line in ifile:
        line = line.split()
        t.append(int(line[0]))
        phi.append(float(line[-1]))
    ifile.close()
    t = np.array(t, dtype = np.int32)
    phi = np.array(phi, dtype = np.float64)
    return t, phi

##################################################################

def compute_msr(time, phi, nmsr, limit):
    """ compute and average the MSR"""
    ### compute required values
    nsteps = len(phi)
    limit = int(limit*nsteps)
    ### allocate output array to store the results
    t_msr = np.zeros((nmsr), dtype = np.int32)
    msr = np.zeros((nmsr), dtype = np.float64)
    counter = np.zeros((nmsr), dtype = np.int32)
    logval = np.zeros((nmsr), dtype = np.int32)
    ### submit MSR computation using the performance tools
    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.compute_MSR_loglog(t_msr, time, phi, msr, counter, logval, nsteps, nmsr, limit)
    ### return time and msr
    return t_msr, msr

##################################################################

def main():
    """ main function, called when the script is started"""
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    print 'reading end-to-end vector'
    time, phi = read_angle(inputfile)
    ### compute local angle phi based on overall rotation
  
    ### compute msr
    print 'computing mean square rotation'
    t_msr, msr = compute_msr(time, phi, nmsr, limit)
    
    #### write tables
    # correlation functions
    ofile = open('msr.data', 'w')
    ofile.write('Mean Square Rotation\n\n')
    ofile.write('lag time\tMSR\n')
    for i in range(len(t_msr)):
        ofile.write(str(t_msr[i]) + '\t' + str(msr[i]) + '\n')
    ofile.close()
  
    ## generate a plot of the data
    fig = plt.figure()
    plt.loglog(t_msr, msr)
    plt.savefig('msr.png')
    plt.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
