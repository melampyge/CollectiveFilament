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
    # compute that many points for the MSR
    line = ifile.readline()
    line = line.split()
    nmsr = int(line[-1])
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
    return inputfile, nmsr, limit, ofname


##################################################################

def neigh_min(dx,lx):
    ### compute minimum distance
    dx1 = dx + lx
    dx2 = dx - lx
    if dx**2 < dx1**2 and dx**2 < dx2**2:
        return dx
    if dx1**2 < dx2**2:
        return dx1
    return dx2

##################################################################

def get_ee(dumpfile):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = open(dumpfile, 'r')
    ### generate lists for ex and ey
    ex = []
    ey = []
    # time
    time = []
    # first run variable
    counter = 0
    ### loop over all timesteps
    while True:
        # check whether end of file has been reached
        line = ifile.readline()
        if line == '':
            break
        # timestep
        line = ifile.readline()
        line = line.split()
        time.append(int(line[0]))
        # natoms
        ifile.readline()
        line = ifile.readline()
        line = line.split()
        natoms = int(line[0])
        if counter == 0:
            x = np.zeros((natoms))
            y = np.zeros((natoms))
        # box dimensions
        ifile.readline()
        line = ifile.readline()
        line = line.split()
        xlo = float(line[0])
        xhi = float(line[1])
        lx = xhi - xlo
        line = ifile.readline()
        line = line.split()
        ylo = float(line[0])
        yhi = float(line[1])
        ly = yhi - ylo
        ifile.readline()
        # skip body headline
        ifile.readline()
        # read in positions
        for i in range(natoms):
            line = ifile.readline()
            line = line.split()
            aID = int(line[0])-1
            xs = (float(line[1]) + float(line[4]))*lx
            ys = (float(line[2]) + float(line[5]))*ly
            x[aID] = xs
            y[aID] = ys
        # first rod atom
        x1 = x[0]
        y1 = y[0]
        # second atom
        x2 = x[1]
        y2 = y[1]
        # fill arrays
        dx = x1 - x2
        dx = neigh_min(dx,lx)
        dy = y1 - y2
        dy = neigh_min(dy,ly)
        ex.append(dx)
        ey.append(dy)

        counter = counter + 1
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    ex = np.array(ex, dtype = np.float64)
    ey = np.array(ey, dtype = np.float64)
    ### return data
    return time, ex, ey

##################################################################

def compute_angle(ex,ey):
    """ compute phi"""
    ### allocate phi array
    nsteps = len(ex)
    phi = np.zeros((nsteps), dtype = np.float64)
    dp = np.zeros((nsteps), dtype = np.float64)
    ### compute angle and correct angle on the fly
    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.compute_angle(ex,ey,phi,dp,nsteps)
    return phi

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
    ### read parameters from input file
    print 'reading settings'
    inputfile, nmsr, limit, ofname = read_settings()
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    print 'reading end-to-end vector'
    time, ex, ey = get_ee(inputfile)
    ### compute local angle phi based on overall rotation
    print 'computing orientation'
    phi = compute_angle(ex,ey)
    
    ### compute msr
    print 'computing mean square rotation'
    t_msr, msr = compute_msr(time, phi, nmsr, limit)
    
    ### generate output folder
    os.system('mkdir ' + str(ofname))
       
    #### write tables
    # correlation functions
    ofile = open(ofname + '/msr.data', 'w')
    ofile.write('Mean Square Rotation\n\n')
    ofile.write('lag time\tMSR\n')
    for i in range(len(t_msr)):
        ofile.write(str(t_msr[i]) + '\t' + str(msr[i]) + '\n')
    ofile.close()
    # phi
    ofile = open(ofname + '/phi.data', 'w')
    ofile.write('Current orientation of the rod\n\n')
    ofile.write('timestep\tphi\n')
    for i in range(len(time)):
        ofile.write(str(time[i]) + '\t' + str(phi[i]) + '\n')
    ofile.close()

    ## generate a plot of the data
    fig = plt.figure()
    plt.loglog(t_msr, msr)
    plt.savefig(ofname + '/msr.png')
    plt.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
