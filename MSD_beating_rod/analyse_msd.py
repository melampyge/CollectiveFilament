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

def get_com(inputfile):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = open(inputfile, 'r')
    ### generate lists for ex and ey
    comx = []
    comy = []
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
        # center of mass
        comx.append(np.average(x))
        comy.append(np.average(y))
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    comx = np.array(comx, dtype = np.float64)
    comy = np.array(comy, dtype = np.float64)
    ### return data
    return time, comx, comy

##################################################################

def compute_msd(time, comx, comy, dt, limit):
    """ compute and average the MSR"""
    ### compute required values
    nsteps = len(comx)
    limit = int(limit*nsteps)
    nmsd = int(limit/dt)
    ### allocate output array to store the results
    t_msd = np.zeros((nmsd), dtype = np.int32)
    msd = np.zeros((nmsd), dtype = np.float64)
    counter = np.zeros((nmsd), dtype = np.int32)
    linval = np.zeros((nmsd), dtype = np.int32)
    ### loop over all molecules, submit MSD computation, read and average
    ###   the results
    performance_tools = performance_toolsWrapper.Performance_tools()
    # submit MSR computation
    performance_tools.compute_MSD(t_msd, time, comx, comy, msd, counter, linval, nsteps, nmsd, limit, dt)
    ### return time and msr
    return t_msd, msd

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
    t_msd, msd = compute_msd(time, comx, comy, dt, limit)
    
    ### generate output folder
    os.system('mkdir ' + str(ofname))
       
    #### write tables
    # correlation functions
    ofile = open(ofname + '/correlation_functions.data', 'w')
    ofile.write('Correlation functions\n\n')
    ofile.write('lag time\tMSD\n')
    for i in range(len(t_msd)):
        ofile.write(str(t_msd[i]) + '\t' + str(msd[i]) + '\n')
    ofile.close()
    # com
    ofile = open(ofname + '/com.data', 'w')
    ofile.write('Current center of mass of the rod\n\n')
    ofile.write('timestep\tcomx\tcomy\n')
    for i in range(len(time)):
        ofile.write(str(time[i]) + '\t' + str(comx[i]) + '\t' + str(comy[i]) + '\n')
    ofile.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
