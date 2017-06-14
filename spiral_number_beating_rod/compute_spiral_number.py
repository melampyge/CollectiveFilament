#!/usr/local/bin/python2.7

import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import os


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
    # output folder
    line = ifile.readline()
    line = line.split()
    ofname = line[-1]
    # close file
    ifile.close()
    # return input values
    return inputfile, ofname

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

def read_data(ifile):
    """ read in the coordinates"""
    # check whether end of file has been reached
    line = ifile.readline()
    if line == '':
        return 0, 0, 0, 0, 0, 1
    # timestep
    line = ifile.readline()
    line = line.split()
    t = int(line[0])
    # number of atoms
    ifile.readline()
    line = ifile.readline()
    line = line.split()
    natoms = int(line[0])
    # box size
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
    # atom coordinates
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    ifile.readline()
    for i in range(natoms):
        line = ifile.readline()
        line = line.split()
        aID = int(line[0]) - 1
        xs = float(line[1])
        ys = float(line[2])
        ix = int(line[4])
        iy = int(line[5])
        x[aID] = lx*(xs + ix)
        y[aID] = ly*(ys + iy)
    return t, x, y, lx, ly, 0

##################################################################

def compute_spiral_number(x,y,lx,ly):
    """ compute the absolute of the spiral number of a filament"""
    nbeads = len(x)
    # compute all bond orientations
    phi = np.zeros((nbeads - 1))
    for i in range(1,nbeads):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        dphi = math.atan2(dy,dx)
        phi[i-1] = dphi
    # corrrect fo 2*pi periodicity
    phi2 = np.copy(phi)
    nbonds = len(phi)
    for i in range(1,nbonds):
        dphi = phi[i] - phi[i-1]
        if dphi < -np.pi:
            dphi += 2*np.pi
        elif dphi > np.pi:
            dphi -= 2*np.pi
        phi2[i] = phi2[i-1] + dphi
    # compute the spiral number
    s = (phi2[-1] - phi2[0])/2/np.pi
    return s

##################################################################

def analyse_data(inputfile):
    """ read in the data and compute the spiral number"""
    # open input file for reading
    ifile = open(inputfile)
    # initialize lists
    time = []
    spiral = []
    # loop over the dumpfile
    count = 0
    while True:
        count = count + 1
        print count
        ti,x,y,lx,ly,eof = read_data(ifile)
        if eof == 1:
            break
        si = compute_spiral_number(x[-100:],y[-100:],lx,ly)
        # add values to output array
        time.append(ti)
        spiral.append(si)
    # close input file
    ifile.close()
    # transform lists to arrays and return them
    time = np.array(time)
    spiral = np.array(spiral)
    return time, spiral


##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    print '  reading settings'
    inputfile, ofname = read_settings()
    ### read in the coordinates for the entire trajectory, discard center of mass movement
    print '  computing the spiral number'
    time, spiral = analyse_data(inputfile)

 
    ### generate output folder
    os.system('mkdir ' + str(ofname))

    #### write tables
    # eigenvalues, eigenvector
    ofile = open(ofname + '/spiral.data', 'w')
    ofile.write('spiral number of the tail of the filament measured over time\n\n')
    for i in range(len(time)):
        ofile.write(str(time[i]) + '\t' + str(spiral[i]) + '\n')
    ofile.close()


            
##################################################################

if __name__ == '__main__':
    main()
    
