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
import equal_time_correlationWrapper

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
    # number of bins
    line = ifile.readline()
    line = line.split()
    nbin = int(line[-1])
    # maximum cutoff
    line = ifile.readline()
    line = line.split()
    rmax = float(line[-1])
    # overall number of snapshots to use in the analysis
    line = ifile.readline()
    line = line.split()
    ntotal = int(line[-1])
    # close file
    ifile.close()
    # return input values
    return charfile, headerfile, ofname, nfil, nskip, nbin, rmax, ntotal

##################################################################

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
            dx = neigh_min(x2-x1,lx)
            dy = neigh_min(y2-y1,ly)
            # compute angle using atan2
            pi = math.atan2(dy,dx)
            phi[k] = pi
            tx[k] = dx / np.sqrt(dx**2 + dy**2)
            ty[k] = dy / np.sqrt(dx**2 + dy**2)
            # increment k
            k = k + 1
    return phi, tx, ty

############################################################################

def gen_linked_list(x,y,lx,ly,dcrit,natoms):
    """ generate a linked list"""
    # determine the number of cells in each direction
    nsegx = int(lx/dcrit)
    nsegy = int(ly/dcrit)

    # allocate head and llist
    ncells = nsegx*nsegy
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

def compute_pol_correlation_single(x, y, nbin, rmax, lx, ly, natoms, mol, tx, ty):
    """ compute the radial distribution function of the center of mass"""
    ### allocate array to store the results
    pol_correlation = np.zeros((nbin), dtype = np.float64)
    counter = np.zeros((nbin), dtype = np.float64)
    ### generate a linked list of the com
    nsegx, nsegy, head, llist = gen_linked_list(x,y,lx,ly,rmax, natoms)
    ### loop over the linked list
    etcorrelation = equal_time_correlationWrapper.Equal_time_correlation()
    etcorrelation.compute(nsegx, nsegy, natoms, head, llist, mol, x, y, pol_correlation, lx, ly, rmax, nbin, tx, ty, counter)
    return pol_correlation, counter

##################################################################

def compute_polarity_correlation(charfile, headerfile, nfil, nskip, nbin, rmax, ntotal):
    """ compute the radial distribution function
        of the center of mass of all filaments"""
    ### open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    ### skip initial steps
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    ### update nsteps based on the number of overall snapshots to use in the analysis
    nomit = nsteps/ntotal -1
    if nomit <= 0:
        nomit = 0
    else:
        nsteps = ntotal
    ### generate information about molecules
    mol, nmol = gen_mol_info(natoms, nfil)
    ### allocate arrays, including output: radial distribution function
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    # time
    time = np.zeros((nsteps))
    # polarity_correlation
    pol_correlation = np.zeros((nbin, nsteps), dtype = np.float64)
    counter_correlation = np.zeros((nbin, nsteps), dtype = np.float64)
    # distance array
    r = np.linspace(0.,rmax,nbin+1)
    # averaged and averaged squared polarity
    tx_av = 0
    ty_av = 0
    tsq_av = 0
    ### loop over all timesteps
    for step in range(nsteps):
        ### print progress to screen
        print 'Current Step / All Steps', step, '/', nsteps
        ### skip nomit nsapshots
        read_char.skip_snapshots(hfile, ifile, nomit)
        ### read in coordinates
        xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
        x = xs*lx
        y = ys*ly
        ### compute the local orientation
        phi, tx, ty = compute_orientation(x,y,lx,ly,nfil)
        ### store time
        time[step] = tstep
        ### compute the radial distribution function for the current snapshot
        pol_corri, counteri = compute_pol_correlation_single(x, y, nbin, rmax, lx, ly, natoms, mol, tx, ty)
        ### copy current polarity correlation to all correlations and current counter to all counters
        for i in range(nbin):
            pol_correlation[i,step] = pol_corri[i]
            counter_correlation[i,step] = counteri[i]
        ### compute averages for the polarity and squared polarity
        tx_av += np.average(tx)
        ty_av += np.average(ty)
        tsq_av += np.average(tx**2 + ty**2)
    ### close the input files
    hfile.close()
    ifile.close()
    ### correct the velocity correlation function
    pol_corr_av = np.average(pol_correlation, axis = 1)
    counter_correlation_av = np.average(counter_correlation, axis = 1)
    tx_av /= nsteps
    ty_av /= nsteps
    tav_sq = tx_av**2 + ty_av**2
    tsq_av /= nsteps
    pol_corr_av = (pol_corr_av - counter_correlation_av*tav_sq) / (counter_correlation_av*(tsq_av - tav_sq))
    return r, pol_corr_av


##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip, nbin, rmax, ntotal = read_settings()
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    r, pol_corr = compute_polarity_correlation(charfile, headerfile, nfil, nskip, nbin, rmax, ntotal)

    ### write results to file and generage a plot
    # generate folder structure
    os.system('mkdir ' + ofname)
    # write data of the averaged polarity correlation
    ofile = open(ofname + '/pol_correlation.data', 'w')
    ofile.write('polarity correlation function\n\n')
    ofile.write('r_min\tr_max\tpol_corrlation\n')
    for i in range(nbin):
        ofile.write(str(r[i]) + '\t' + str(r[i+1]) + '\t' + str(pol_corr[i]) + '\n')
    ofile.close()
    # gen figure of the averaged polarity correlation
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(0.5*(r[:-1]+r[1:]), pol_corr)
    ax.set_xlabel(r'r [$\sigma]')
    ax.set_ylabel(r'g_p(r)')
    ax.set_title('Polarity Correlation Function')
    plt.savefig(ofname + '/pol_correlation.png')
    plt.close()
    return
    
##################################################################

if __name__ == '__main__':
    main()
    
