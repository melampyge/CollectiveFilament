#!/usr/users/iff_th2/isele/Applications/Anaconda/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import sys
import math
import os
import codecs
from scipy import optimize
import read_char
import performance_toolsWrapper

try:
    ifname = sys.argv[1]

except:
    print 'Usage: ' + sys.argv[0] + '     parameter file'
    exit()
    
##################################################################

def read_settings():
    """ read in the settings from the parameter file"""
    # open file for reading
    ifile = open(ifname, 'r')
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
    # output-folder name
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
    # critical distance for mol/mol neighbor search
    line = ifile.readline()
    line = line.split()
    dcrit = float(line[-1])
    # critical overlap for mol/mol neibhor search
    line = ifile.readline()
    line = line.split()
    lcrit = float(line[-1])
    # critical angle for mol/mol neighbor search
    line = ifile.readline()
    line = line.split()
    pcrit = float(line[-1])
    pcrit *= np.pi/180.0
    ccrit = math.cos(pcrit)
    # critical width for the histogram
    line = ifile.readline()
    line = line.split()
    wcrit = float(line[-1])
    # close file
    ifile.close()
    # return input values
    return charfile, headerfile, ofname, nfil, nskip,dcrit, lcrit, ccrit, wcrit

############################################################################

def gen_mol_info(natoms, nfil):
    """ generate information about molecular ID"""
    mol = np.zeros((natoms), dtype = np.int32)
    nmol = natoms/nfil
    k = 0
    for i in range(nmol):
        for j in range(nfil):
            mol[k] = i
            k = k + 1
    return mol,nmol

############################################################################

def nearest_neighbor(x1,x2,lx):
    """ compute vector to nearest neighbor"""
    dx1 = x1 - x2
    dx2 = x1 - x2 + lx
    dx3 = x1 - x2 - lx
    if dx1**2 < dx2**2 and dx1**2 < dx3**2:
        return dx1
    if dx2**2 < dx3**2:
        return dx2
    return dx3

############################################################################

def compute_orientation(x,y,lx,ly,nfil):
    """ compute orientation of all beads from bond vectores"""
    # number of molecules
    natoms = len(x)
    nmol = natoms/nfil
    # allocate aray for results
    phi = np.zeros((natoms), dtype = np.float64)
    tx = np.zeros((natoms), dtype = np.float64)
    ty = np.zeros((natoms), dtype = np.float64)
    # loop over all polymers
    k = 0
    for i in range(nmol):
        for j in range(nfil):
            if j == 0:
                x1 = x[k]
                y1 = y[k]
                x2 = x[k+1]
                y2 = y[k+1]
            elif j == nfil-1:
                x1 = x[k-1]
                y1 = y[k-1]
                x2 = x[k]
                y2 = y[k]
            else:
                x1 = x[k-1]
                y1 = y[k-1]
                x2 = x[k+1]
                y2 = y[k+1]
            # compute nearest neighbor
            dx = nearest_neighbor(x1,x2,lx)
            dy = nearest_neighbor(y1,y2,ly)
            # compute angle using atan2
            pi = math.atan2(dy,dx)
            phi[k] = pi
            tx[k] = dx / np.sqrt(dx**2 + dy**2)
            ty[k] = dy / np.sqrt(dx**2 + dy**2)
            # increment k
            k = k + 1
    return phi, tx, ty

############################################################################

def gen_linked_list(x,y,lx,ly,dcrit):
    """ generate a linked list"""
    # determine the number of cells in each direction
    nsegx = int(lx/dcrit)
    nsegy = int(ly/dcrit)

    # allocate head and llist
    ncells = nsegx*nsegy
    natoms = len(x)
    head = np.zeros((ncells), dtype = np.int32)
    llist = np.zeros((natoms), dtype = np.int32)

    # fill list and head
    for i in range(natoms):
        segx = int(x[i]/lx*nsegx)
        segy = int(y[i]/ly*nsegy)
        cell = segx*nsegy + segy
        llist[i] = head[cell]
        head[cell] = i

    return nsegx,nsegy,head,llist

##################################################################

def find_clusters(x,y,tx,ty,mol,nmol,nfil,lx,ly,dcrit,lcrit,ccrit):
    """ search for clusters; two molecules are defined as being
        part of the same cluster if lcrit of their body length is
        within a distance of dcrit and when the difference in orientation
        is less than phi"""
    natoms = len(x)
    # allocate required arrays
    cl = np.zeros((nmol), dtype = np.int32)
    cl -= 1
    clusters = np.zeros((natoms), dtype = np.int32)
    neighs_mol = np.zeros((nmol,nmol), dtype = np.int32)
    # generate a linked list
    nsegx, nsegy, head, llist = gen_linked_list(x,y,lx,ly,dcrit)
    # find neighborhood information
    neighs_molf = neighs_mol.ravel()
    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.fill_neigh_matrix(neighs_molf,llist,head,nsegx,nsegy,x,y,tx,ty,mol,nmol,lx,ly,dcrit,ccrit)
    # recursive search for clusters in neighbor matrix
    performance_tools.cluster_search(neighs_molf,cl,lcrit,nfil,nmol)
    # fill cluster results to per atom array
    k = 0
    for i in range(nmol):
        for j in range(nfil):
            clusters[k] = cl[i]
            k = k + 1
    ### return results
    return clusters

##################################################################

def cluster_analysis(charfile, headerfile, ofname, nfil, nskip, dcrit, lcrit, ccrit, wcrit):
    """ perform cluster analysis"""
    ### open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    nsteps = nsteps - nskip
    if nsteps <= 0:
        print 'Error, nsteps < 0'
        exit()
    ### generate information about molecule IDs
    mol, nmol = gen_mol_info(natoms, nfil)
    ### allocate arrays required for processing the data


    ### allocat arrays to store the cluster statistics
    ncl = np.zeros((nsteps), dtype = np.int32)
    clsize = [] # unfortunately a list    
    t = np.zeros((nsteps), dtype = np.int32)
    clmax = np.zeros((nsteps), dtype = np.int32)

    ### skip some initial snapshots
    read_char.skip_snapshots(hfile, ifile, nskip)
    
    ### loop over all steps
    for i in range(nsteps):
        ### print stats
        print 'Progress:',i,'/',nsteps
        ### read in the data
        xs,ys,lx,ly,tstep,natoms = read_char.read_snapshot(hfile, ifile)
        x = xs*lx
        y = ys*ly
        ### compute orientation of each bead, incorporate pbc
        phi,tx,ty = compute_orientation(x,y,lx,ly,nfil)
        ### find clusters
        clusters_atm = find_clusters(x,y,tx,ty,mol,nmol,nfil,lx,ly,dcrit,lcrit,ccrit)

        ### store atom scalars
        os.system('mkdir ' + ofname + '/cluster_atm')
        ofile = open(ofname + '/cluster_atm/t_' + str(tstep) + '.data', 'w')
        ofile.write('Cluster value for each atom\n\n')
        ofile.write('aID\tclID\n')
        for i in range(natoms):
            ofile.write(str(i) + '\t' + str(clusters_atm[i]) + '\n')
        ofile.close()

    # close input files
    ifile.close()
    hfile.close()
    # return list with cluster sizes, maybe more in the future
    return 

##################################################################

def main():
    """ main function, called when script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip, dcrit, lcrit, ccrit, wcrit = read_settings()
    ### analyse the clusters of the simulation
    cluster_analysis(charfile, headerfile, ofname, nfil, nskip, dcrit, lcrit, ccrit, wcrit)
    
    return

##################################################################

if __name__ == '__main__':
    main()
