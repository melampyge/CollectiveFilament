#!/usr/local/bin/python2.7

import numpy as np
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
import sys
import math
import os
import pma
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
    # compute crosscorrelation for every this many snapshots
    line = ifile.readline()
    line = line.split()
    dt = int(line[-1])
    # timestep up to which to compute the cross correlation
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

def get_coords(inputfile):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = open(inputfile, 'r')
    ### generate lists for coordinates x and y
    xall = []
    yall = []
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
        # add coordinates to array
        xall.append(np.copy(x))
        yall.append(np.copy(y))
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    xall = np.array(xall, dtype = np.float64)
    yall = np.array(yall, dtype = np.float64)
    ### return data
    return time, xall, yall

##################################################################

def compute_curvature(xall,yall):
    """ compute curvature"""
    ### allocate curvature array
    nsteps = len(xall)
    natoms = len(xall[0])
    curv = np.zeros((nsteps, natoms-2), dtype = np.float64)
    ### compute the local curvature
    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.compute_curvature(curv, xall, yall, nsteps, natoms)

    return curv

##################################################################

def cross_correlate(time, a1, a2, dt, limit):
    """ compute the cross correlation between signals a1 and a2"""
    ### compute required values
    nsteps = len(time)
    limit = int(limit*nsteps)
    ncc = int(limit/dt)
    ### allocate output array to store the results
    t_cc = np.zeros((ncc), dtype = np.int32)
    cc1 = np.zeros((ncc), dtype = np.float64)
    cc2 = np.zeros((ncc), dtype = np.float64)
    counter = np.zeros((ncc), dtype = np.int32)
    linval = np.zeros((ncc), dtype = np.int32)
    ### compute the cross correlation function using a c++ code
    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.compute_crosscorrelation(t_cc, time, a1, a2, cc1, cc2, counter, linval, nsteps, ncc, limit, dt)
    ### return t_cc and cc
    return t_cc, cc1, cc2


##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    print 'reading settings'
    inputfile, dt, limit, ofname = read_settings()
    ### read in the coordinates for the entire trajectory, discard center of mass movement
    print 'reading in coordinates'
    time, xall, yall = get_coords(inputfile)
    ### compute the curvature
    print 'computing curvature'
    curv = compute_curvature(xall, yall)

    ### perform principal component analysis
    print 'performing pca'
    eval, evec, correlation_matrix = pma.pma(curv)

    ### project curvatures onto 2 principal coordinates
    print 'computing amplitudes'
    evsel = evec[:,:2]
    amplitudes = pma.get_XY_fast(curv, evsel)
    amp_pca1 = amplitudes[0]
    amp_pca2 = amplitudes[1]
    amp_pca1 = np.array(amp_pca1, dtype = np.float64)
    amp_pca2 = np.array(amp_pca2, dtype = np.float64)
    
    ### compute the autocorrelation and cross-correlation of the 2 dominant amplitudes
    print 'computing cross-correlations'
    t_a1ac, a1ac, a1ac2 = cross_correlate(time, amp_pca1, amp_pca1, dt, limit)
    t_a2ac, a2ac, a2ac2 = cross_correlate(time, amp_pca2, amp_pca2, dt, limit)
    t_a1a2cc, a1a2cc, a2a1cc = cross_correlate(time, amp_pca1, amp_pca2, dt, limit)
   
    ### generate output folder
    os.system('mkdir ' + str(ofname))

    #### write tables
    # eigenvalues, eigenvector
    ofile = open(ofname + '/eigen.data', 'w')
    ofile.write('eigenvalues and eigenvectors (sorted by column) of the correlation matrix\n\n')
    neig = len(eval)
    ofile.write('n = ' + str(neig) + '\n')
    # eigenvalues
    ofile.write('eval')
    for i in range(neig):
        ofile.write('\t' + str(eval[i]))
    ofile.write('\n')
    # eigenvectors
    for i in range(neig):
        ofile.write('c' + str(i))
        for j in range(neig):
            ofile.write('\t' + str(evec[j,i]))
        ofile.write('\n')
    ofile.close()

    # correlation matrix
    ofile = open(ofname + '/correlation_matrix.data', 'w')
    ofile.write('Non-normalized correlation matrix\n\n')
    ofile.write('n = ' + str(neig) + '\n')
    for i in range(neig):
        for j in range(neig):
            ofile.write(str(correlation_matrix[i,j]) + '\t')
        ofile.write('\n')
    ofile.close()
    
    # amplitudes
    ofile = open(ofname + '/amplitudes.data', 'w')
    ofile.write('Amplitudes of the first two principal components\n\n')
    ofile.write('tstep\ta1\ta2\n')
    nsteps = len(time)
    for i in range(nsteps):
        ofile.write(str(time[i]) + '\t' + str(amp_pca1[i]) + '\t' + str(amp_pca2[i]) + '\n')
    ofile.close()

    # correlation functions
    ofile = open(ofname + '/amplitude_correlations.data', 'w')
    ofile.write('Auto- and cross-correlation functions of the amplitudes\n\n')
    ofile.write('lag time\ta1ac\ta2ac\ta1a2cc\ta2a1cc\n')
    nsteps = len(t_a1ac)
    for i in range(nsteps):
        ofile.write(str(t_a1ac[i]) + '\t' + str(a1ac[i]) + '\t' + str(a2ac[i]) + '\t' + str(a1a2cc[i]) + '\t' + str(a2a1cc[i]) + '\n')
    ofile.close()
            
##################################################################

if __name__ == '__main__':
    main()
    
