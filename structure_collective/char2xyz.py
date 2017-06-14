
import numpy as np
import sys
import math
import codecs
import read_char
import os
import h5py

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
    
    # char-file name
    line = ifile.readline()
    line = line.split()
    charfile = line[-1]
    
    # header-file name
    line = ifile.readline()
    line = line.split()
    headerfile = line[-1]
    
    # length of the filaments
    line = ifile.readline()
    line = line.split()
    nfil = int(line[-1])
    
    # box bound
    line = ifile.readline()
    line = line.split()
    box = float(line[-1])
    
    # close file
    ifile.close()
    
    return charfile, headerfile, nfil, box

##################################################################

def neigh_min(dx,lx):
    """ compute minimum distance"""
    
    dx1 = dx + lx
    dx2 = dx - lx
    if dx**2 < dx1**2 and dx**2 < dx2**2:
        return dx
    if dx1**2 < dx2**2:
        return dx1
    return dx2

##################################################################

def write_xyz(x, y, nsteps, natoms, filename):
    """ write the data into xyz file"""
    
    print "WRITING TO XYZ"
    
    ### open files for writing
    
    f = open(filename, 'w')
    
    ### write to file
    
    for step in range(nsteps):
        
        f.write(str(natoms) + '\n\n')
        
        for j in range(natoms):
            
            f.write('S\t' + str(x[step][j]) + '\t' + str(y[step][j]) + '\t 0\n')
    
    f.close()
    
    return 
    
##################################################################

def write_hdf5(time, x, y, box, nsteps, natoms, nmol, filename):
    """ write the data into hdf5 file"""

    print "WRITING TO HDF5"
    
    ### open files for writing
    
    f = h5py.File(filename, 'w')
    
    ### write to file
    
    # position info
    pos = f.create_group('pos')
    pos.create_dataset('x', (nsteps,natoms), data=x, dtype='float64', compression='gzip')
    pos.create_dataset('y', (nsteps,natoms), data=y, dtype='float64', compression='gzip')
    
    # simulation info
    info = f.create_group('info')
    boxinfo = info.create_group('box')
    boxinfo.create_dataset('lx', data=box, dtype='float')
    boxinfo.create_dataset('ly', data=box, dtype='float')
    info.create_dataset('nsteps', data=nsteps)
    info.create_dataset('natoms', data=natoms)
    info.create_dataset('nmol', data=nmol)
    
    # time info
    f.create_dataset('time', data=time, dtype='float64', compression='gzip')
        
    f.close()
    
    return    
    
##################################################################

def read_input(charfile, headerfile, nfil):
    """ load the position data at once"""
    
    ### open files for reading
    
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    
    ### get general information from header file
    
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
        
    ### allocate arrays
    
    time = np.zeros((nsteps))
    x = np.zeros((nsteps,natoms))
    y = np.zeros((nsteps,natoms))
            
    ### read all the coordinates into a coordinate array
    
    for step in range(nsteps):
        
        ### print progress to screen
        
        print 'READING - Current Step / All Steps', step, '/', nsteps
        
        ### read in coordinates
        
        xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
        time[step] = tstep
        x[step] = xs*lx
        y[step] = ys*ly
         
    ### close the input files

    hfile.close()
    ifile.close()   
    
    return time, x, y, nsteps, natoms, nmol

##################################################################

def main():
    
    ### read parameters from input file
    
    charfile, headerfile, nfil, box = read_settings()
    
    ### read the input
    
    time, x, y, nsteps, natoms, nmol = read_input(charfile, headerfile, nfil)
    
    ### write to xyz
    
    xyzfile = 'out1.xyz'
    write_xyz(x, y, nsteps, natoms, xyzfile)
    
    ### write to hdf5
    
    hdf5file = 'out1.hdf5'
    write_hdf5(time, x, y, box, nsteps, natoms, nmol, hdf5file)
    
    return
    
    
##################################################################

if __name__ == '__main__':
    main()
    
##################################################################