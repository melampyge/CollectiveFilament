#!/usr/local/bin/python2.7

import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
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
    # sigma for gauss filter
    line = ifile.readline()
    line = line.split()
    sigma = int(line[-1])
    # output folder
    line = ifile.readline()
    line = line.split()
    ofname = line[-1]
    # close file
    ifile.close()
    # return input values
    return inputfile, sigma, ofname

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

def first_pass(inputfile):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = open(inputfile, 'r')
    # first run variable
    counter = 0
    ### loop over all timesteps
    while True:
        ### read a snapshot
        # check whether end of file has been reached
        line = ifile.readline()
        if line == '':
            break
        # timestep
        line = ifile.readline()
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
            #xs = (float(line[1]) + float(line[4]))*lx
            #ys = (float(line[2]) + float(line[5]))*ly
            xs = float(line[1])*lx
            ys = float(line[2])*ly
            x[aID] = xs
            y[aID] = ys
        # correct periodic boundaries
        correct_pbc(x, lx)
        correct_pbc(y, ly)

        ### compute the curvatures
        if counter == 0:
            curv = np.zeros((natoms - 2), dtype = np.float64)
        performance_tools = performance_toolsWrapper.Performance_tools()
        performance_tools.compute_curvature_single(curv, x, y, natoms)

        ### add the data to the correlation matrix
        if counter == 0:
            corr = np.zeros((natoms-2, natoms -2), dtype = np.float64)
        corr += np.outer(curv, curv)
        
        counter = counter + 1
    # close the input file
    ifile.close()
    # normalize the correlation matrix
    corr /= counter - 1
    ### return correlation matrix
    return corr

##################################################################

def second_pass(inputfile, eigvec, ofname):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = open(inputfile, 'r')
    # first run variable
    counter = 0
    ### loop over all timesteps
    while True:
        ### read a snapshot
        # check whether end of file has been reached
        line = ifile.readline()
        if line == '':
            break
        # timestep
        line = ifile.readline()
        line = line.split()
        t = int(line[0])
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
            #xs = (float(line[1]) + float(line[4]))*lx
            #ys = (float(line[2]) + float(line[5]))*ly
            xs = float(line[1])*lx
            ys = float(line[2])*ly
            x[aID] = xs
            y[aID] = ys
        # correct periodic boundaries
        correct_pbc(x, lx)
        correct_pbc(y, ly)

        ### compute the curvatures
        if counter == 0:
            curv = np.zeros((natoms - 2), dtype = np.float64)
        performance_tools = performance_toolsWrapper.Performance_tools()
        performance_tools.compute_curvature_single(curv, x, y, natoms)

        ### compute the amplitudes 
        amp = eigvec.T.dot(curv.T)

        ### store the amplitudes, times, and maximum amplitudes
        if counter == 0:
            of1 = open(ofname + '/amplitudes.data', 'w')
            of1.write('Amplitude of the first 10 modes\n\n')
            of1.write('timestep')
            for i in range(10):
                of1.write('\tamp_' + str(i))
            of1.write('\n')
            of2 = open(ofname + '/max_amplitudes.data', 'w')
            of2.write('Maximum amplitude\n\n')
            of2.write('timestep\tamp_id\tamp_value\n')
        # amplitudes
        of1.write(str(t))
        for i in range(10):
            of1.write('\t' + str(amp[i]))
        of1.write('\n')
        # maximum amplitudes
        idmax = np.argmax(amp**2)
        ampmax = amp[idmax]
        of2.write(str(t) + '\t' + str(idmax) + '\t' + str(ampmax) + '\n')
        counter = counter + 1
    # close the input file
    ifile.close()
    of1.close()
    of2.close()
    ### return correlation matrix
    return

##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    print '  reading settings'
    inputfile, sigma, ofname = read_settings()

    ### generate output folder
    os.system('mkdir ' + str(ofname))
    
    ### read in the coordinates for the entire trajectory, discard center of mass movement
    print '  first pass: reading coordinates, computing curvatures, computing covariance matrix'
    covar = first_pass(inputfile)

    ### compute the principal modes and eigenvalues
    print '  eigenvalue decomposition'
    eigval, eigvec = np.linalg.eigh(covar)
    
	# sort the eingen values and modes by ascending order of eigen value
    idxx = np.argsort(eigval)[::-1]
    eigval = eigval[idxx]
    eigvec = eigvec[:,idxx]

    ### compute the amplitudes in a second pass
    print '  second pass: computing the amplitudes, store on the fly: amplitudes 1 -- 10, id of the maximum amplitude'
    second_pass(inputfile, eigvec, ofname)

   
    #### write tables
    # eigenvalues, eigenvector
    ofile = open(ofname + '/eigen.data', 'w')
    ofile.write('eigenvalues and eigenvectors (sorted by column) of the correlation matrix\n\n')
    neig = len(eigval)
    ofile.write('n = ' + str(neig) + '\n')
    # eigenvalues
    ofile.write('eigval')
    for i in range(neig):
        ofile.write('\t' + str(eigval[i]))
    ofile.write('\n')
    # eigenvectors
    for i in range(neig):
        ofile.write('c' + str(i))
        for j in range(neig):
            ofile.write('\t' + str(eigvec[j,i]))
        ofile.write('\n')
    ofile.close()

    # correlation matrix
    ofile = open(ofname + '/correlation_matrix.data', 'w')
    ofile.write('Non-normalized correlation matrix\n\n')
    ofile.write('n = ' + str(neig) + '\n')
    for i in range(neig):
        for j in range(neig):
            ofile.write(str(covar[i,j]) + '\t')
        ofile.write('\n')
    ofile.close()
            
##################################################################

if __name__ == '__main__':
    main()
    
