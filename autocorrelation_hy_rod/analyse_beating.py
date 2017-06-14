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

def correct_pbc(z,lz):
    """ correct periodic boundary problem"""
    nz = len(z)
    for i in range(nz-1):
        # put z1 to initial box
        z1 = z[-1-i]
        z1 /= lz
        z1 = z1 - math.floor(z1)
        z1 *= lz
        # put z2 to initial box
        z2 = z[-2-i]
        z2 /= lz
        z2 = z2 - math.floor(z2)
        z2 *= lz
        # compute minimum distance
        dz = neigh_min(z1-z2,lz)
        z[-2-i] = z[-1-i] - dz
    return

##################################################################

def read_hy(inputfile):
    """ read the y-position of the last bead"""
    ### open files for reading
    ifile = open(inputfile, 'r')
    ### generate lists for hy
    hy = []
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

        # correct pbcs, start in reverse order
        correct_pbc(x,lx)
        correct_pbc(y,ly)
        # compute signed hy
        x0 = x[-1]
        x1 = x[0]
        x2 = x[1]
        y0 = y[-1]
        y1 = y[0]
        y2 = y[1]
        signed_dist = ((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)
        signed_dist = signed_dist / math.sqrt((y2-y1)**2 + (x2-x1)**2)
        hy.append(signed_dist)
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
    counter = np.zeros((nac), dtype = np.int32)
    linval = np.zeros((nac), dtype = np.int32)
    ### loop over all molecules, submit MSD computation, read and average
    ###   the results
    performance_tools = performance_toolsWrapper.Performance_tools()
    # submit MSR computation
    performance_tools.compute_autocorrelation(t_ac, time, x, ac, counter, linval, nsteps, nac, limit, dt)
    ### return time and msr
    return t_ac, ac

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
    t_hyac, hyac = compute_autocorrelation(time, hy, dt, limit)
    
    ### generate output folder
    os.system('mkdir ' + str(ofname))
       
    #### write tables
    # correlation functions
    ofile = open(ofname + '/correlation_functions.data', 'w')
    ofile.write('Correlation functions\n\n')
    ofile.write('lag time\tHYAC\n')
    for i in range(len(t_hyac)):
        ofile.write(str(t_hyac[i]) + '\t' + str(hyac[i]) + '\n')
    ofile.close()
    # hy
    ofile = open(ofname + '/hy.data', 'w')
    ofile.write('Current "height" of end position\n\n')
    ofile.write('timestep\thy\n')
    for i in range(len(time)):
        ofile.write(str(time[i]) + '\t' + str(hy[i]) + '\n')
    ofile.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
