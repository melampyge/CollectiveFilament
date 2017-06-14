#!/usr/local/bin/python2.7

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys, os
import math
import read_char
import codecs

try:
    infile = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '     parameter file'
    exit()
    
#################################################################

def read_settings():
    """ read in the settings from the parameter file"""
    # open file for reading
    ifile = open(infile, 'r')
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
    # number of histogram bins
    line = ifile.readline()
    line = line.split()
    nbins = int(line[-1])
    # minimum volume
    line = ifile.readline()
    line = line.split()
    vmin = float(line[-1])
    # maximum volume
    line = ifile.readline()
    line = line.split()
    vmax = float(line[-1])
    # close file
    ifile.close()
    # return input values
    return charfile, headerfile, ofname, nfil, nskip, nbins, vmin, vmax

##################################################################

def gen_mol_info(natoms, nfil):
    """ generate information about molecular ID"""
    mol = np.zeros((natoms), dtype = int)
    nmol = natoms/nfil
    k = 0
    for i in range(nmol):
        for j in range(nfil):
            mol[k] = i
            k = k + 1
    return mol,nmol

##################################################################

def voronoi_tessellation(x,y,mol,lx,ly,natoms,nmol):
    """ compute the molecular volumes using Voronoi tessellation"""
    ### generate input script for voro++
    ofile = open('voro.data', 'w')
    for i in range(natoms):
        ofile.write(str(i) + ' ' + str(x[i]) + ' ' + str(y[i]) + ' 0.5\n')
    ofile.close()
    ### perform Voronoi tessellation using voro++
    os.system('/usr/users/iff_th2/isele/Applications/voro++-0.4.6/src/voro++ -p 0.0 ' + str(lx) + ' 0.0 ' + str(ly) + ' 0.0 1.0 voro.data')
    ### read in the results
    vol = np.zeros((nmol))
    ifile = open('voro.data.vol')
    for i in range(natoms):
        line = ifile.readline()
        line = line.split()
        idx = int(line[0])
        v = float(line[4])
        vol[mol[idx]] += v
    ifile.close()
    ### remove voro++ files
    os.system('rm voro.data voro.data.vol')
    return vol

##################################################################

def analyse_data(charfile, headerfile, nfil, nskip, nbins, vmin, vmax):
    """ analyse the data using the Voronoi approach"""
    ### open input files, get first line from header
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile)
    natoms, nsteps = read_char.read_first(hfile)
    ### skip initial snapshots
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    ### allocate arrays to store the results
    hist_vol = np.zeros((nsteps, nbins))
    ### compute information about molecules
    mol, nmol = gen_mol_info(natoms, nfil)
    
    ### loop over all snapshots
    for i in range(nsteps):
        # print stats
        print 'Progress:',i,'/',nsteps
        # read in the data
        x,y,lx,ly,tstep,natoms = read_char.read_snapshot(hfile, ifile)
        # scale x and y coordinates
        x *= lx
        y *= ly
        # compute molecular volume using Voronoi tessellation
        vol = voronoi_tessellation(x,y,mol,lx,ly,natoms,nmol)
        # bin molecular volumes to histogram
        hist_vol[i], edges = np.histogram(np.log(vol), bins = nbins, range = (vmin,vmax))
    ### close input files
    ifile.close()
    hfile.close()
    ### compute vol axis, averages and std
    vol = 0.5*(edges[:-1] + edges[1:])
    hist_vol_av = np.average(hist_vol, axis = 0)
    hist_vol_std = np.std(hist_vol, axis = 0) / np.sqrt(nsteps)


    return edges, vol, hist_vol_av, hist_vol_std


##################################################################

def main():
    """ main function, called when script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip, nbins, vmin, vmax = read_settings()
    ### start Voronoi analysis
    edges, vol, hist_vol_av, hist_vol_std = analyse_data(charfile, headerfile, nfil, nskip, nbins, vmin, vmax)

    ### write results to files and create a plot
    # create output folder
    os.system('mkdir ' + ofname)
    ### generate a plot of the histogram
    fig = plt.figure()
    ax = plt.subplot(111)
    # add data
    ax.errorbar(vol, hist_vol_av, hist_vol_std)
    # modify the axes
    ax.set_xlabel(r'$ln(a/\sigma^2)$')
    ax.set_ylabel(r'$p(a)$')
    # save, store, and close the plot
    plt.savefig(ofname + '/voronoi_area.png')
    plt.close()

    ### write data to file
    ofile = open(ofname + '/voronoi_area.data', 'w')
    ofile.write('Average Voronoi area per molecule\n\n')
    ofile.write('left binsize\tright binsize\tp(area)\tstd\n')
    for i in range(nbins):
        ofile.write(str(edges[i]) + '\t' + str(edges[i+1]) + '\t' + str(hist_vol_av[i]) + '\t' + str(hist_vol_std[i]) + '\n')
    ofile.close()

    return

##################################################################

if __name__ == '__main__':
    main()
