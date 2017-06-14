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
    # close file
    ifile.close()
    # return input values
    return charfile, headerfile, ofname, nfil, nskip, nbin, rmax

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

##################################################################

def compute_velocity(x1,y1,x3,y3,lx,ly):
    """ compute the particle velocity based on finite difference approach"""
    natoms = len(x1)
    vx = np.zeros((natoms), dtype = np.float64)
    vy = np.zeros((natoms), dtype = np.float64)
    for i in range(natoms):
        vx[i] = neigh_min(x3[i] - x1[i], lx)
        vy[i] = neigh_min(y3[i] - y1[i], ly)
    return vx, vy

##################################################################

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

def compute_vel_correlation_single(x, y, nbin, rmax, lx, ly, natoms, mol, vx, vy):
    """ compute the radial distribution function of the center of mass"""
    ### allocate array to store the results
    vel_correlation = np.zeros((nbin), dtype = np.float64)
    counter = np.zeros((nbin), dtype = np.float64)
    ### generate a linked list of the com
    nsegx, nsegy, head, llist = gen_linked_list(x,y,lx,ly,rmax, natoms)
    ### loop over the linked list
    etcorrelation = equal_time_correlationWrapper.Equal_time_correlation()
    etcorrelation.compute(nsegx, nsegy, natoms, head, llist, mol, x, y, vel_correlation, lx, ly, rmax, nbin, vx, vy, counter)
    return vel_correlation, counter

##################################################################

def compute_velocity_correlation(charfile, headerfile, nfil, nskip, nbin, rmax):
    """ compute the radial distribution function
        of the center of mass of all filaments"""
    ### open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file, correct nsteps by 2
    natoms, nsteps = read_char.read_first(hfile)
    ### skip the initial snapshots
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - 2 - nskip
    ### generate information about molecules
    mol, nmol = gen_mol_info(natoms, nfil)
    ### allocate arrays, including output: radial distribution function
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    # polarity_correlation
    vel_correlation = np.zeros((nbin, nsteps), dtype = np.float64)
    counter_correlation = np.zeros((nbin, nsteps), dtype = np.float64)
    # distance array
    r = np.linspace(0.,rmax,nbin+1)
    # averaged and averaged squared polarity
    vx_av = 0
    vy_av = 0
    vsq_av = 0
    ### read the coordinates of the first two steps
    xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
    x1 = xs*lx
    y1 = ys*ly
    xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
    x2 = xs*lx
    y2 = ys*ly
    ### loop over all timesteps
    for step in range(nsteps):
        ### print progress to screen
        print 'Current Step / All Steps', step, '/', nsteps
        ### read in coordinates
        xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
        x3 = xs*lx
        y3 = ys*ly
        ### compute the local velocity
        vx, vy = compute_velocity(x1,y1,x3,y3,lx,ly)
        ### compute the radial distribution function for the current snapshot
        vel_corri, counteri = compute_vel_correlation_single(x2, y2, nbin, rmax, lx, ly, natoms, mol, vx, vy)
        ### copy current polarity correlation to all correlations and current counter to all counters
        for i in range(nbin):
            vel_correlation[i,step] = vel_corri[i]
            counter_correlation[i,step] = counteri[i]
        ### compute averages for the polarity and squared polarity
        vx_av += np.average(vx)
        vy_av += np.average(vy)
        vsq_av += np.average(vx**2 + vy**2)
        ### shift data
        x1 = np.copy(x2)
        y1 = np.copy(y2)
        x2 = np.copy(x3)
        y2 = np.copy(y3)
    ### close the input files
    hfile.close()
    ifile.close()
    ### correct the velocity correlation function
    vel_corr_av = np.average(vel_correlation, axis = 1)
    counter_correlation_av = np.average(counter_correlation, axis = 1)
    vx_av /= nsteps
    vy_av /= nsteps
    vav_sq = vx_av**2 + vy_av**2
    vav = np.sqrt(vav_sq)
    vsq_av /= nsteps
    vel_corr_av = (vel_corr_av - counter_correlation_av*vav_sq) / (counter_correlation_av*(vsq_av - vav_sq))
    return r, vel_corr_av, vav, vsq_av


##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip, nbin, rmax = read_settings()
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    r, vel_corr, v, vsq = compute_velocity_correlation(charfile, headerfile, nfil, nskip, nbin, rmax)

    ### write results to file and generage a plot
    # generate folder structure
    os.system('mkdir ' + ofname)
    # write data of the averaged polarity correlation
    ofile = open(ofname + '/vel_correlation.data', 'w')
    ofile.write('velocity correlation function\n\n')
    ofile.write('v_av = ' + str(v) + '\n')
    ofile.write('v_sq_av = ' + str(vsq) + '\n')
    ofile.write('r_min\tr_max\tvel_corrlation\n')
    for i in range(nbin):
        ofile.write(str(r[i]) + '\t' + str(r[i+1]) + '\t' + str(vel_corr[i]) + '\n')
    ofile.close()
    # gen figure of the averaged polarity correlation
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(0.5*(r[:-1]+r[1:]), vel_corr)
    ax.set_xlabel(r'r [$\sigma]')
    ax.set_ylabel(r'g_v(r)')
    ax.set_title('Velocity Correlation Function')
    plt.savefig(ofname + '/vel_correlation.png')
    plt.close()
    return
    
##################################################################

if __name__ == '__main__':
    main()
    
