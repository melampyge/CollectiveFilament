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
    limit = float(sys.argv[2])
except:
    print 'Usage: ' + sys.argv[0] + '       inupt file               length correlation / length simulation'
    exit()

##################################################################

def read_angle(inputfile, k):
    """ read in the angles from the input file"""
    ifile = open('orientations.data')
    t = []
    phi = []
    for i in range(3):
        ifile.readline()
    for line in ifile:
        line = line.split()
        t.append(int(line[0]))
        phi.append(float(line[k]))
    ifile.close()
    t = np.array(t, dtype = np.int32)
    phi = np.array(phi, dtype = np.float64)
    return t, phi

##################################################################

def compute_correlation(time, tx, ty, limit):
    """ compute and average the MSR"""
    ### compute required values
    nsteps = len(tx)
    limit = int(limit*nsteps)
    ### allocate output array to store the results
    t_ee = np.zeros((limit), dtype = np.int32)
    c_ee = np.zeros((limit), dtype = np.float64)
    linval = np.zeros((limit), dtype = np.int32)
    std = np.zeros((limit), dtype = np.float64)
    blockdata = np.zeros((nsteps), dtype = np.float64)
    ### size of the array for the blocking method
    nblock = int(math.log(nsteps, 2)) + 2
    block_std = np.zeros((nblock), dtype = np.float64)
    block_uncert = np.zeros((nblock), dtype = np.float64)   
    ### submit MSR computation using the performance tools
    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.compute_correlation_with_errorbars(t_ee, time, tx, ty, c_ee, std, blockdata, block_std, block_uncert, linval, nsteps, limit)
    ### return time and msr
    return t_ee, c_ee, std

##################################################################

def main():
    """ main function, called when the script is started"""
    for k in [3]:
        ### loop over the entire trajectory and compute the center of mass, correct for pbc
        print 'reading end-to-end vector'
        time, phi = read_angle(inputfile, k)
        ### transform phi to tx and ty
        print 'transforming phi to tx and ty'
        n = len(phi)
        tx = np.zeros((n), dtype = np.float64)
        ty = np.zeros((n), dtype = np.float64)
        for i in range(n):
            tx[i] = np.cos(phi[i])
            ty[i] = np.sin(phi[i])
  
        ### compute autocorrelation
        print 'computing autocorrelation function'
        t_ee, c_ee, std = compute_correlation(time, tx, ty, limit)
    
        #### write tables
        # correlation functions
        ofile = open('e2e_autocorrelation' + str(k) + '_std.data', 'w')
        ofile.write('Autocorrelation of the orientation vector\n\n')
        ofile.write('lag time\tc_ee\tstd\n')
        for i in range(len(t_ee)):
            ofile.write(str(t_ee[i]) + '\t' + str(c_ee[i]) + '\t' + str(std[i]) + '\n')
        ofile.close()
  
        ## generate a plot of the data
        fig = plt.figure()
        plt.errorbar(t_ee, c_ee, 3*std)
        plt.savefig('c_ee' + str(k) + '.png')
        plt.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
