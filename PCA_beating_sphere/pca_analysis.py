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
from scipy import ndimage

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
    # beads to ignore
    line = ifile.readline()
    line = line.split()
    nhead = int(line[-1])
    # output folder
    line = ifile.readline()
    line = line.split()
    ofname = line[-1]
    # close file
    ifile.close()
    # return input values
    return inputfile, sigma, nhead, ofname

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

def get_coords(inputfile, nhead):
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
            x = np.zeros((natoms - nhead + 2), dtype = np.float64)
            y = np.zeros((natoms - nhead + 2), dtype = np.float64)
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
        # read in positions, neglected position of most of the head molecules
        for i in range(natoms):
            line = ifile.readline()
            line = line.split()
            aID = int(line[0])
            if aID > 2 and aID <= nhead:
                continue
            if aID > nhead:
                aID = aID - nhead + 2
            aID = aID - 1
            #xs = (float(line[1]) + float(line[4]))*lx
            #ys = (float(line[2]) + float(line[5]))*ly
            xs = float(line[1])*lx
            ys = float(line[2])*ly
            x[aID] = xs
            y[aID] = ys
        # correct pbc
        correct_pbc(x, lx)
        correct_pbc(y, ly)
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

def blocking_method(data):
    """ compute the std and average of some data using the blocking method"""
    n = len(data)
    blockdata = np.copy(data)
    """
    block_std = np.zeros((100), dtype = np.float64)
    block_uncert = np.zeros((100), dtype = np.float64)
    av = 0.0
    std = 0.0
    performance_tools = performance_toolsWrapper.Performance_tools()
    av, std = performance_tools.blocking_method(blockdata, block_std, block_uncert, n, av, std)
    return av, std 
    """
    block_std = []
    block_uncert = []

    av = np.average(blockdata)
    
    while n >= 2:
        # compute std for current block
        c0 = 0
        for i in range(n):
            c0 += (blockdata[i] - av)**2
        c0 /= (n-1)
        # compute the standard deviation and fill it to blocking array
        c0 = np.sqrt(c0/(n-1))
        block_std.append(c0)
        block_uncert.append(c0/np.sqrt(2*(n-1)))
        # perform block operation
        for i in range(n/2):
            blockdata[i] = 0.5*(blockdata[2*i] + blockdata[2*i+1])
        n = n/2

    # set initial value for the std
    std = 0.0
    ## find the plateau value and fill it to stdl, set to -1 if no plateau is found
    nignore = 3;
    n = len(block_std) - nignore
    smin_old = block_std[0] - 2*block_uncert[0]
    smax_old = block_std[0] + 2*block_uncert[0]
    for i in range(1,n): 
        smin = block_std[i] - 2*block_uncert[i]
        smax = block_std[i] + 2*block_uncert[i]
        if (smax_old > smin):
	        std = 0.5*(block_std[i] + block_std[i-1])
	        break
        smin_old = smin
        smax_old = smax

    return av, std

##################################################################

def curvature_averages(curv):
    """ compute the average and std of each curvature component
        using the blocking method"""
    # allocate arrays
    ncurv = len(curv[0,:])
    curv_av = np.zeros((ncurv), dtype = np.float64)
    curv_std = np.zeros((ncurv), dtype = np.float64)
    # compute average and std using the blocking method
    for i in range(ncurv):
        curv_av[i], curv_std[i] = blocking_method(curv[:,i])
        
    return curv_av, curv_std

##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    print '  reading settings'
    inputfile, sigma, nhead, ofname = read_settings()
    ### read in the coordinates for the entire trajectory, discard center of mass movement
    print '  reading in coordinates'
    time, xall, yall = get_coords(inputfile, nhead)
    ### compute the curvature
    print '  computing curvature'
    curv = compute_curvature(xall, yall)

    ### apply gauss filter to the curvature coordinates
    print '  computing the gauss filter'
    curv_filtered = ndimage.filters.gaussian_filter1d(curv, sigma, axis = 0)

    ### compute average and std of the curvature
    print '  computing the curvature averages'
    curv_av, curv_std = curvature_averages(curv_filtered) 
    
    ### perform principal component analysis
    print '  performing pca'
    eval, evec, correlation_matrix, curv_av = pma.pma_uncentered(curv)

    ### project curvatures onto 2 principal coordinates
    print '  computing amplitudes'
    evsel = evec[:,:]
    amplitudes = pma.get_XY_fast(curv_filtered, evsel)

    ### store the amplitude witht the maximum value at each timestep
    nsteps = len(time)
    amp_id = np.zeros((nsteps), dtype = np.int32)
    amp_max = np.zeros((nsteps), dtype = np.float64)
    for i in range(nsteps):
        idmax = np.argmax(amplitudes[:,i]**2)
        amp_id[i] = idmax
        amp_max[i] = amplitudes[idmax, i]

    ### compute average energy that is stored in eacha amplitude
    amp_energy = np.average(amplitudes**2, axis = 1)   
   
    ### generate output folder
    os.system('mkdir ' + str(ofname))

    #### write tables
    # average curvature
    ofile = open(ofname + '/average_shape.data', 'w')
    ofile.write('Average shape of the filament\n\n')
    ofile.write('angle\tcurvature\tstd\n')
    for i in range(len(curv_av)):
        ofile.write(str(i) + '\t' + str(curv_av[i]) + '\t' + str(curv_std[i]) + '\n')
    ofile.close()
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
    namp_store = 10
    nsteps = len(time)
    ofile = open(ofname + '/amplitudes.data', 'w')
    ofile.write('Amplitudes of all principal components\n\n')
    ofile.write('tstep')
    for i in range(namp_store):
        ofile.write('\ta' + str(i))
    ofile.write('\n')
    for i in range(nsteps):
        ofile.write(str(time[i]))
        for j in range(namp_store):
            ofile.write('\t' + str(amplitudes[j,i]))
        ofile.write('\n')
    ofile.close()

    # amplitude energy
    ofile = open(ofname + '/amplitude_energies.data', 'w')
    ofile.write('Amplitudes squared and averaged over time (= average energy stored in each mode)\n\n')
    ofile.write('mode\tamp_sq\n')
    for i in range(len(amp_energy)):
        ofile.write(str(i) + '\t' + str(amp_energy[i]) + '\n')
    ofile.close()

    # maximum ampmlitudes
    ofile = open(ofname + '/max_amplitudes.data', 'w')
    ofile.write('Amplitudes of principal components with maximum value\n\n')
    ofile.write('tstep\tamp_id\tamp_val\n')
    for i in range(nsteps):
        ofile.write(str(time[i]) + '\t' + str(amp_id[i]) + '\t' + str(amp_max[i]) + '\n')
    ofile.close()
            
##################################################################

if __name__ == '__main__':
    main()
    
