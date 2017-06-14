#!/usr/local/bin/python2.7

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import os
import codecs
import read_char

try:
    infilename = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '      paraneter file'
    exit()
    

##################################################################
#    global variables

# number of fractions into which to decompose the box
frac = [100., 70., 50., 30., 20., 10., 7., 5., 3., 2.]
frac = np.array(frac)
nfrac = len(frac)

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

def compute_number_fluctuations(x,y,n):
    """ compute number fluctuations"""
    ### compute histogram of the number of bins in each cell
    hist, xedges, yedges = np.histogram2d(x,y,bins=n, range = [[0,1],[0,1]])
    ### compute the standard deviation
    std = np.std(hist)
    dn = std
    return dn

##################################################################

def analyse_data(charfile, headerfile, nfil, nskip):
    """ analyse the simulation data"""
    ### open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
    ### skip the first snapshots
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    ### generate output array
    dn = np.zeros((nsteps,nfrac))
    ### loop over the file
    for step in range(nsteps):
        ### print progress to screen
        print 'Current Step / All Steps', step, '/', nsteps
        ### read in coordinates
        xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
        # compute the number fluctuations for all bin widths
        #   and add results to output array
        for j in range(nfrac):
            dnj = compute_number_fluctuations(xs,ys,frac[j])
            dn[step,j] = dnj
    ### close input files
    ifile.close()
    hfile.close()
    ### compute averages and standard deviations
    dn_av = np.average(dn, axis = 0)
    dn_std = np.average(dn, axis = 0)/np.sqrt(nsteps)
    ### compute an array with the average number of beads in each bin
    nbeads = float(natoms)/frac**2

    return nbeads, dn_av, dn_std


##################################################################

def main():
    """ main function, called when script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip = read_settings()
    ### loop over trajectory and analyse the data
    nbeads, dn_av, dn_std = analyse_data(charfile, headerfile, nfil, nskip)

    ### write results to files
    os.system('mkdir -p ' + ofname)
#    ### generate a log-log plot of the data
#    fig = plt.figure()
#    ax = plt.subplot(111)
#    # add measured values
#    ax.errorbar(nbeads, dn_av, yerr = dn_std, ls = '', marker = 'o', label = 'measured')
#    # generate lines with slope 1 and 2
#    k05 = dn_av[0]/np.sqrt(nbeads[0])
#    k1 = dn_av[0] / nbeads[0]
#    ax.plot(nbeads, k05*np.sqrt(nbeads), color = '0.3', label = 'slope = 0.5')
#    ax.plot(nbeads, k1*nbeads, color = '0.6', label = 'slope = 1')
#    # add legend
#    ax.legend(loc = 'lower right')
#    # modify the axes
#    ax.set_xscale('log')
#    ax.set_yscale('log')
#    ax.set_xlabel(r'$n$')
#    ax.set_ylabel(r'$\Delta n$')
#    # save and show the results
#    plt.savefig(ofname + '/number_fluc.png')
#    plt.close()
    ### write data to a table
    ofile = open(ofname + '/number_fluc.data', 'w')
    ofile.write('Measured Number Fluctuations\n\n')
    ofile.write('n\tdelta_n\tstd\n')
    for i in range(nfrac):
        ofile.write(str(nbeads[i]) + '\t' + str(dn_av[i]) + '\t' + str(dn_std[i]) + '\n')
    ofile.close()

    
    return

##################################################################

if __name__ == '__main__':
    main()
