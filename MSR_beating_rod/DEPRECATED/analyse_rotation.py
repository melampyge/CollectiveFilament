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
    # filetype (dump or xyz)
    line = ifile.readline()
    line = line.split()
    filetype = line[-1]
    if filetype not in ['dump', 'xyz']:
        print 'filetype must be dump or xyz'
        exit()
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
    return inputfile, filetype, dt, limit, ofname

##################################################################

def get_ee_xyz(xyzfile):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = open(xyzfile, 'r')
    ### generate lists for ex and ey
    ex = []
    ey = []
    # time
    time = []
    ### loop over all timesteps
    while True:
        # number of atoms
        line = ifile.readline()
        if line == '':
            break
        line = line.split()
        natoms = int(line[0])
        # timestep
        line = ifile.readline()
        line = line.split()
        time.append(int(line[-1]))
        # first rod atom
        line = ifile.readline()
        line = line.split()
        x1 = float(line[1])
        y1 = float(line[2])
        # second atom
        line = ifile.readline()
        line = line.split()
        x2 = float(line[1])
        y2 = float(line[2])
        # fill arrays
        ex.append(x1-x2)
        ey.append(y1-y2)
        # all other atoms
        for i in range(natoms-2):
            ifile.readline()
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    ex = np.array(ex, dtype = np.float64)
    ey = np.array(ey, dtype = np.float64)
    ### return data
    return time, ex, ey

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

def get_ee_dump(dumpfile):
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
        # correct pbcs, start in reverse order
        correct_pbc(x,lx)
        correct_pbc(y,ly)
        # first rod atom
        x1 = x[0]
        y1 = y[0]
        # second atom
        x2 = x[1]
        y2 = y[1]
        # fill arrays
        ex.append(x1-x2)
        ey.append(y1-y2)
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    ex = np.array(ex, dtype = np.float64)
    ey = np.array(ey, dtype = np.float64)
    ### return data
    return time, ex, ey

##################################################################

def get_ee(inputfile, filetype):
    """ get ee using either the method for xyz or dump files"""
    if filetype == 'dump':
        time, ex, ey = get_ee_dump(inputfile)
    if filetype == 'xyz':
        time, ex, ey = get_ee_xyz(inputfile)
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

def compute_msr(time, phi, dt, limit):
    """ compute and average the MSR"""
    ### compute required values
    nsteps = len(phi)
    limit = int(limit*nsteps)
    nmsr = int(limit/dt)
    ### allocate output array to store the results
    t_msr = np.zeros((nmsr), dtype = np.int32)
    msr = np.zeros((nmsr), dtype = np.float64)
    counter = np.zeros((nmsr), dtype = np.int32)
    linval = np.zeros((nmsr), dtype = np.int32)
    ### loop over all molecules, submit MSD computation, read and average
    ###   the results
    performance_tools = performance_toolsWrapper.Performance_tools()
    # submit MSR computation
    performance_tools.compute_MSR(t_msr, time, phi, msr, counter, linval, nsteps, nmsr, limit, dt)
    ### return time and msr
    return t_msr, msr

##################################################################

def compute_msr_with_derivatives(time, phi, dphi, d2phi, dt, limit):
    """ compute and average the MSR"""
    ### compute required values
    nsteps = len(phi)
    limit = int(limit*nsteps)
    nmsr = int(limit/dt)
    ### allocate output array to store the results
    t_msr = np.zeros((nmsr), dtype = np.int32)
    msr = np.zeros((nmsr), dtype = np.float64)
    dmsr = np.zeros((nmsr), dtype = np.float64)
    d2msr = np.zeros((nmsr), dtype = np.float64)
    d2msr_t1 = np.zeros((nmsr), dtype = np.float64)
    d2msr_t2 = np.zeros((nmsr), dtype = np.float64)
    counter = np.zeros((nmsr), dtype = np.int32)
    linval = np.zeros((nmsr), dtype = np.int32)
    ### loop over all molecules, submit MSD computation, read and average
    ###   the results
    performance_tools = performance_toolsWrapper.Performance_tools()
    # submit MSR computation
    performance_tools.compute_MSR_with_derivatives(t_msr, time, phi, dphi, d2phi, msr, dmsr, d2msr, d2msr_t1, d2msr_t2, counter, linval, nsteps, nmsr, limit, dt)
    ### return time and msr
    return t_msr, msr, dmsr, d2msr, d2msr_t1, d2msr_t2

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
    inputfile, filetype, dt, limit, ofname = read_settings()
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    print 'reading end-to-end vector'
    time, ex, ey = get_ee(inputfile, filetype)
    ### compute local angle phi based on overall rotation
    print 'computing orientation'
    phi = compute_angle(ex,ey)
    # compute derivative of phi
    dphi, d2phi = compute_derivatives(phi)
    
    ### compute msr
    print 'computing mean square rotation'
    t_msr, msr = compute_msr(time, phi, dt, limit)
    ### compute the mean square rotational velocity???? (actually nonsense, compute it for testing purposes)
    #print 'computing measn square rotational velocity'
    #t_msrv, msrv = compute_msr(time, dphi, dt, limit)
    ### computing the mean square rotational accelleration??? (actually nonsense, compute it for testing purposes)
    #print 'computing mean square rotational acceleration'
    #t_msra, msra = compute_msr(time, d2phi, dt, limit)
    ### compute the autocorrelation function of the rotational velocity
    print 'computing the rotational velocity autocorrelation function'
    t_vac, vac = compute_autocorrelation(time, dphi, dt, limit)
    ### compute the autocorrelation function of the rotational accelleration
    #print 'computing the rotational acceleration autocorrelation function'
    #t_aac, aac = compute_autocorrelation(time, d2phi, dt, limit)
    
    ### generate output folder
    os.system('mkdir ' + str(ofname))
       
    #### write tables
    # correlation functions
    ofile = open(ofname + '/correlation_functions.data', 'w')
    ofile.write('Correlation functions\n\n')
    ofile.write('lag time\tMSR\tVAC\n')
    for i in range(len(t_msr)):
        ofile.write(str(t_msr[i]) + '\t' + str(msr[i]) + '\t' + str(vac[i]) + '\n')
    ofile.close()
    # phi
    ofile = open(ofname + '/phi.data', 'w')
    ofile.write('Current orientation of the rod\n\n')
    ofile.write('timestep\tphi\n')
    for i in range(len(time)):
        ofile.write(str(time[i]) + '\t' + str(phi[i]) + '\n')
    ofile.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
