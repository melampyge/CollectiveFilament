#!/usr/local/bin/python2.7

# output: tx and ty of the head
#         tangent autocorrelation function
#         position of the center bead of the sphere

import numpy as np
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
    # compute tangent autocorrelation for every this many snapshots
    line = ifile.readline()
    line = line.split()
    dt = int(line[-1])
    # timestep up to which to compute the autocorrelation
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

def get_tangent_center(dumpfile):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = open(dumpfile, 'r')
    ### generate lists for ex and ey
    tx = []
    ty = []
    xc = []
    yc = []
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
        dx = x1 - x2
        dy = y1 - y2
        d = np.sqrt(dx**2 + dy**2)
        tx.append(dx/d)
        ty.append(dy/d)
        xc.append(x1)
        yc.append(y1)
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    tx = np.array(tx, dtype = np.float64)
    ty = np.array(ty, dtype = np.float64)
    xc = np.array(xc, dtype = np.float64)
    yc = np.array(yc, dtype = np.float64)
    ### return data
    return time, tx, ty, xc, yc



##################################################################

def compute_ttac(time, tx, ty, dt, limit):
    """ compute the autocorrelation function of x"""
    ### compute required values
    nsteps = len(tx)
    limit = int(limit*nsteps)
    nac = int(limit/dt)
    ### allocate output array to store the results
    t_ac = np.zeros((nac), dtype = np.int32)
    ac = np.zeros((nac), dtype = np.float64)
    linval = np.zeros((nac), dtype = np.int32)
    std = np.zeros((nac), dtype = np.float64)
    blockdata = np.zeros((nsteps), dtype = np.float64)
    ### size of the array for the blocking method
    nblock = 100;
    block_std = np.zeros((nblock), dtype = np.float64)
    block_uncert = np.zeros((nblock), dtype = np.float64)
    ### use cpp code to compute the correlation function
    performance_tools = performance_toolsWrapper.Performance_tools()
    # submit MSR computation
    performance_tools.compute_ttac_with_errorbars(t_ac, time, tx, ty, ac, std, blockdata, block_std, block_uncert, linval, nsteps, nac, limit, dt)

    ### return time and msr
    return t_ac, ac, std


##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    print 'reading settings'
    inputfile, dt, limit, ofname = read_settings()
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    print 'reading relevant information'
    time, tx, ty, xc, yc = get_tangent_center(inputfile)
    ### compute local angle phi based on overall rotation
    print 'compute tangent-tangent autocorrelation function'
    tc, ttac, ttstd = compute_ttac(time, tx, ty, dt, limit)

    
    ### generate output folder
    os.system('mkdir ' + str(ofname))
       
    #### write tables
    # tangent vectors
    ofile = open(ofname + '/tangent_vectors.data', 'w')
    ofile.write('Tangent vector of the head\n\n')
    ofile.write('t\ttx\tty\n')
    for i in range(len(time)):
        ofile.write(str(time[i]) + '\t' + str(tx[i]) + '\t' + str(ty[i]) + '\n')
    ofile.close()
    # center bead
    ofile = open(ofname + '/center_position.data', 'w')
    ofile.write('Position of the center bead (bead that connects the head and the tail)\n\n')
    ofile.write('t\tx\ty\n')
    for i in range(len(time)):
        ofile.write(str(time[i]) + '\t' + str(xc[i]) + '\t' + str(yc[i]) + '\n')
    ofile.close()
    # auto-correlation functions
    ofile = open(ofname + '/tangent_correlation.data', 'w')
    ofile.write('Tangent-tangent time autocorrelation\n\n')
    ofile.write('lag time\ttac\tstd\n')
    for i in range(len(tc)):
        ofile.write(str(tc[i]) + '\t' + str(ttac[i]) + '\t' + str(ttstd[i]) + '\n')
    ofile.close()

        
##################################################################

if __name__ == '__main__':
    main()
    
