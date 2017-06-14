#!/usr/local/bin/python2.7

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import codecs
import read_char
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
    # char-file name
    line = ifile.readline()
    line = line.split()
    charfile = line[-1]
    # header-file name
    line = ifile.readline()
    line = line.split()
    headerfile = line[-1]
    # output folder
    line = ifile.readline()
    line = line.split()
    ofname = line[-1]
    # length of the filaments
    line = ifile.readline()
    line = line.split()
    nfil = int(line[-1])
    # number of snapshots to skip
    line = ifile.readline()
    line = line.split()
    nskip = int(line[-1])
    # close file
    ifile.close()
    # return input values
    return charfile, headerfile, ofname,  nfil, nskip

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

def compute_ee_single(x,y,lx,ly):
    """ compute the absolute of the spiral number of a filament,
        correct pbc on the fly"""
    nbeads = len(x)
    # correct pbcs
    for i in range(1,nbeads):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        x[i] = x[i-1] + neigh_min(dx,lx)
        y[i] = y[i-1] + neigh_min(dy,ly)
    # compute ex and ey
    ex = x[-1] - x[0]
    ey = y[-1] - y[0]
    return ex, ey

##################################################################

def get_ee(charfile, headerfile, nfil, nskip):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
    ### skip initial snapshots
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    ### allocate arrays, including output: spiral histogram + spiral counter
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    # center of mass coordinates
    ex = np.zeros((nmol, nsteps))
    ey = np.zeros((nmol, nsteps))
    # time
    time = np.zeros((nsteps))
    ### loop over all timesteps
    for step in range(nsteps):
        ### print progress to screen
        print 'Current Step / All Steps', step, '/', nsteps
        ### read in coordinates
        xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
        x = xs*lx
        y = ys*ly
        ### store time
        time[step] = tstep
        ### extract differences in endpoint coordinates for each filament, correct pbc on the fly
        for i in range(nmol):
            ex[i,step], ey[i,step] = compute_ee_single(x[i*nfil:(i+1)*nfil-1], y[i*nfil:(i+1)*nfil-1], lx, ly)
    ### close the input files
    hfile.close()
    ifile.close()
    return time, ex, ey

##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip = read_settings()
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    time, ex, ey = get_ee(charfile, headerfile, nfil, nskip)
    ### compute average squared end-to-end distance
    rsq = ex**2+ ey**2
    rav = np.sqrt(np.average(rsq))
    rstd = np.sqrt(np.std(rsq)/np.sqrt(len(ex)*len(ex[0])))
    ### write to file
    os.system('mkdir ' + ofname)
    ofile = open(ofname + '/e2e.data', 'w')
    ofile.write('Average RMS End-to-end distance, averaged over total time and all molecules\n\n')
    ofile.write('r_av = ' + str(rav) + '\n')
    ofile.write('r_std = ' + str(rstd) + '\n')

    
##################################################################

if __name__ == '__main__':
    main()
    
