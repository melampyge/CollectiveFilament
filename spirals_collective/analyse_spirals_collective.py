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
    return charfile, headerfile, ofname, nfil, nskip

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

def compute_spiral_number(x,y,lx,ly):
    """ compute the absolute of the spiral number of a filament,
        correct pbc on the fly"""
    nbeads = len(x)
    # correct pbcs
    for i in range(1,nbeads):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        x[i] = x[i-1] + neigh_min(dx,lx)
        y[i] = y[i-1] + neigh_min(dy,ly)
    # compute all bond orientations
    phi = np.zeros((nbeads - 1))
    for i in range(1,nbeads):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        dphi = math.atan2(dy,dx)
        phi[i-1] = dphi
    # correct for 2*pi periodicity
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

def analyse_spirals(charfile, headerfile, nfil, nskip):
    """ search for spirals"""
    ### open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
    ### skip the initial snapshots
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    ### allocate arrays, including output: spiral histogram + spiral counter
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    spiral_number = np.zeros((nmol))
    # time
    time = np.zeros((nsteps))
    # histogram values
    smax = 12.0
    smin = 0.0
    ds = 0.05    
#    smax = 5.0
#    smin = 0.0
#    ds = 0.05
    nbins = (smax - smin)/ds
    spiral_histogram = np.zeros((nbins))
    edges = np.linspace(smin, smax, nbins+1, endpoint=True) 
    # number of filaments with more than a certain spiral number
    nspiral_05 = np.zeros((nsteps))
    nspiral_10 = np.zeros((nsteps))
    nspiral_15 = np.zeros((nsteps))
    nspiral_20 = np.zeros((nsteps))
    # track the spiral number of a tagged filament in time
    stag = np.zeros((nsteps))
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
        ### compute spiral number for each filament, correct pbc on the fly
        for i in range(nmol):
            spiral_number[i] = compute_spiral_number(x[i*nfil:(i+1)*nfil-1], y[i*nfil:(i+1)*nfil-1], lx, ly)
        ### add computed spiral numbers to histogram an spiral counter
        for i in range(nmol):
            s = math.fabs(spiral_number[i])
            if s >= 0.5:
                nspiral_05[step] += 1
            if s >= 1.0:
                nspiral_10[step] += 1
            if s >= 1.5:
                nspiral_15[step] += 1
            if s >= 2.0:
                nspiral_20[step] += 1
            seg = s/(smax-smin)*nbins
            if seg >= nbins:
                seg = nbins - 1
            spiral_histogram[seg] += 1
        stag[step] = spiral_number[10]
    ### close the input files
    hfile.close()
    ifile.close()
    ### normalize spiral counts with the number of filaments
    nspiral_05 /= nmol
    nspiral_10 /= nmol
    nspiral_15 /= nmol
    nspiral_20 /= nmol
    return nspiral_05, nspiral_10, nspiral_15, nspiral_20, edges, spiral_histogram, stag, time

##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip = read_settings()
    ### loop over the entire trajectory and count spirals
    nspiral_05, nspiral_10, nspiral_15, nspiral_20, edges, spiral_histogram, stag, time = analyse_spirals(charfile, headerfile, nfil, nskip)
    ### write results to files
    # create output folder
    os.system('mkdir ' + ofname)
    # spiral evolution
    ofile = open(ofname + '/spiral_evolution.data', 'w')
    ofile.write('Fraction of the filaments with a spiral number larger than a threshold\n\n')
    ofile.write('timestep\ts>=0.5\ts>=1.0\ts>=1.5\ts>=2.0\n')
    for i in range(len(time)):
        ofile.write(str(time[i]) + '\t' + str(nspiral_05[i]) + '\t' + str(nspiral_10[i]) + '\t' + str(nspiral_15[i]) + '\t' + str(nspiral_20[i]) + '\n')
    ofile.close()
    # tagged mol spiral number
    ofile = open(ofname + '/spiral_tag.data', 'w')
    for i in range(len(time)):
        ofile.write(str(time[i]) + '\t\t' + str(stag[i]) + '\n')
    ofile.close()
    # spiral histogram
    ofile = open(ofname + '/spiral_histogram.data', 'w')
    ofile.write('Histogram of the spiral numbers\n\n')
    ofile.write('lower_edge\tupper_edge\tnumber\n')
    for i in range(len(spiral_histogram)):
        ofile.write(str(edges[i]) + '\t' + str(edges[i+1]) + '\t' + str(spiral_histogram[i]) + '\n')
    ofile.close()
    
##################################################################

if __name__ == '__main__':
    main()
    
