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

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    print '  reading settings'
    inputfile, sigma, ofname = read_settings()
    ### read in the coordinates for the entire trajectory, discard center of mass movement
    print '  reading in coordinates'
    time, xall, yall = get_coords(inputfile)
    ### compute the curvature
    print '  computing curvature'
    curv = compute_curvature(xall, yall)

    ### apply gauss filter to the curvature coordinates
    print '  computing the gauss filter'
    curv_filtered = ndimage.filters.gaussian_filter1d(curv, sigma, axis = 0)

    ### perform principal component analysis
    print '  performing pca'
    eval, evec, correlation_matrix = pma.pma_uncentered(curv)

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
   
    ### generate output folder
    os.system('mkdir ' + str(ofname))

    #### write tables
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
    
