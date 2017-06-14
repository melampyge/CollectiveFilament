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
    nmsd = int(sys.argv[2])
    limit = float(sys.argv[3])
except:
    print 'Usage: ' + sys.argv[0] + '       input file          #(MSD points)         length MSD / length simulation'
    exit()


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
        comx.append(float(line[-2]))
        comy.append(float(line[-1]))
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    comx = np.array(comx, dtype = np.float64)
    comy = np.array(comy, dtype = np.float64)
    ### return data
    return time, comx, comy

##################################################################

def compute_msd(time, comx, comy, nmsd, limit):
    """ compute and average the MSR"""
    ### compute required values
    nsteps = len(comx)
    limit = int(limit*nsteps)
    ### allocate output array to store the results
    t_msd = np.zeros((nmsd), dtype = np.int32)
    msd = np.zeros((nmsd), dtype = np.float64)
    counter = np.zeros((nmsd), dtype = np.int32)
    logval = np.zeros((nmsd), dtype = np.int32)
    ### loop over all molecules, submit MSD computation, read and average
    ###   the results
    performance_tools = performance_toolsWrapper.Performance_tools()
    # submit MSR computation
    performance_tools.compute_MSD(t_msd, time, comx, comy, msd, counter, logval, nsteps, nmsd, limit)
    ### return time and msr
    return t_msd, msd

##################################################################

def main():
    """ main function, called when the script is started"""
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    print 'reading end-to-end vector'
    time, comx, comy = get_com(inputfile)
    
    ### compute msd
    print 'computing mean square displacement'
    t_msd, msd = compute_msd(time, comx, comy, nmsd, limit)
    
    #### write tables
    # correlation functions
    ofile = open('msd.data', 'w')
    ofile.write('Mean Square Displacement\n\n')
    ofile.write('lag time\tMSD\n')
    for i in range(len(t_msd)):
        ofile.write(str(t_msd[i]) + '\t' + str(msd[i]) + '\n')
    ofile.close()

    ### generate a figure
    plt.loglog(t_msd, msd)
    plt.savefig('msd.png')
    plt.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
