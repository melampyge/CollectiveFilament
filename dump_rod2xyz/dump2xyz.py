#!/usr/bin/python

import numpy as np
import sys
import math
import os
import matplotlib.pyplot as plt

try:
    infilename = sys.argv[1]
    outfilename = sys.argv[2]
except:
    print 'Usage: ' + sys.argv[0] + '       input file       output file'
    exit()

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

def copy_snapshots(infilename, outfilename):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = open(infilename, 'r')
    ofile = open(outfilename, 'w')
    ### generate lists for ex and ey
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
        tstep = int(line[0])
        # natoms
        ifile.readline()
        line = ifile.readline()
        line = line.split()
        natoms = int(line[0])
        if counter == 0:
            x = np.zeros((natoms), dtype = np.float64)
            y = np.zeros((natoms), dtype = np.float64)
        # box dimensions
        ifile.readline()
        line = ifile.readline()
        line = line.split()
        xlo = np.float(line[0])
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
            x[aID] = float(line[2])*lx
            y[aID] = float(line[3])*ly
#            x[aID] = (float(line[1]) + float(line[4]))*lx
#            y[aID] = (float(line[2]) + float(line[5]))*ly
        # correct pbcs, start in reverse order
        correct_pbc(x,lx)
        correct_pbc(y,ly)
        #### write data to xyz file
        ofile.write(str(natoms) + '\n')
        ofile.write('Atoms. Timestep: ' + str(tstep) + '\n')
        for i in range(natoms):
            ofile.write('1\t' + str(x[i]) + '\t' + str(y[i]) + '\t0\n')
        counter = counter + 1
    # close the input file
    ifile.close()
    ofile.close()
   
    return

##################################################################

def main():
    """ main function, called when the script is started"""
    copy_snapshots(infilename, outfilename)
 
    
##################################################################

if __name__ == '__main__':
    main()
    
