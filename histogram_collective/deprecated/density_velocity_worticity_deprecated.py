#!/usr/local/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import math

try:
    fname = sys.argv[1]
    nsteps = int(sys.argv[2])
    nbins = int(sys.argv[3])
    lbox = float(sys.argv[4])
    nhist = int(sys.argv[5])
except:
    print 'Usage: ' + sys.argv[0] + '      infilename       nsteps        nbins       lbox      nhist'
    exit()
    

##################################################################

def read_coordinates(ifile):
    """ read in the coordinates"""
    # read the number of atoms from the header, skip rest of header
    for i in range(3):
        ifile.readline()
    line = ifile.readline()
    line = line.split()
    natoms = int(line[0])
    # read in the box dimension
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
    # skip rest of header
    ifile.readline()
    ifile.readline()

    # allocate arrays
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    mol = np.zeros((natoms), dtype = int)

    # read in the coordinates
    for i in range(natoms):
        line = ifile.readline()
        line = line.split()
        aID = int(line[0]) - 1
        xs = float(line[2])
        ys = float(line[3])
        xs = xs - math.floor(xs)
        ys = ys - math.floor(ys)
        x[aID] = xs*lx
        y[aID] = ys*ly

    # define the molecule
    nmol = natoms/npol
    k = 0
    for i in range(nmol):
        for j in range(npol):
            mol[k] = i
            k = k + 1
    return x,y,mol,lx,ly

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

def compute_orientation(x,y,lx,ly):
    """ compute orientation of all beads from bond vectores"""
    # number of molecules
    natoms = len(x)
    nmol = natoms/npol
    # allocate aray for results
    phi = np.zeros((natoms))
    # loop over all polymers
    k = 0
    for i in range(nmol):
        for j in range(npol):
            if j == 0:
                x1 = x[k]
                y1 = y[k]
                x2 = x[k+1]
                y2 = y[k+1]
            elif j == npol-1:
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
            # increment k
            k = k + 1
    return phi

##################################################################

def compute_histogram(x,y, rhomin, rhomax):
    """ compute the density in each segment"""
    # create a list of bins
    bins = np.zeros((nbins,nbins))
    binwidth = lbox/nbins
    # loop over all particles, increment bin counter
    natoms = len(x)
    for i in range(natoms):
        xi = x[i]
        yi = y[i]
        j = math.floor(xi/binwidth)
        k = math.floor(yi/binwidth)
        bins[j,k] += 1
    # create a colormap of the bins
    #plt.imshow(bins.transpose(), interpolation = 'None')
    #plt.plot(x/binwidth-0.5,y/binwidth-0.5,ls = '', marker = 'o', markersize = 1)
    #plt.axis('equal')
    #plt.colorbar()
    #plt.show()
    #plt.close()
    # counter to density
    bins /= rhomax
    # flatten density
    bins = bins.flatten()
    # generate histogram
    hist, edges = np.histogram(bins, bins = nhist, range = (0,1))
    #plt.hist(bins, bins = nhist, range = (0,1))
    #plt.show()
    #plt.close()
    return hist

##################################################################

def compute_density_histogram():
    """ compute a histogram of the density"""
    # allocate an array to store the data
    binwidth = lbox/nbins
    rhomin = 0.0
    rhomax = binwidth**2
    hist = np.zeros((nsteps,nhist))
    # open the file for reading
    ifile = open(fname, 'r')
    # loop over all snapshots 
    for i in range(nsteps):
        print i
        # read in the coordinates
        x,y = read_coordinates(ifile)
        # compute the histogram for a given snapshot
        h = compute_histogram(x,y,rhomin,rhomax)
        # fill values to target array
        for j in range(nhist):
            hist[i,j] = h[j]
    # close file
    ifile.close()
    # return histogram
    return hist

##################################################################

def block_av(data):
    """ compute block averages"""
    n1 = len(data)
    n2 = len(data[0])
    nblock = 100
    n3 = n1/nblock
    data_av = np.zeros((n3,n2))
    # compute averages
    m = 0
    for i in range(n3):
        for k in range(nblock):
            for j in range(n2):
                data_av[i,j] += data[m,j]
            m = m + 1
    data_av /= nblock
    return data_av

##################################################################

def analyse_data():
    """ analyse the simulation data"""
    # open input file
    ifile = open(fname)
    # loop over all steps
    for i in range(nsteps):
        # read in the coordinates
        x,y,mol,lx,ly = read_coordinates(ifile)
        # compute orientation of all beads
        phi = compute_orientation(x,y,lx,ly)
        # generate values for histogram edges
        edges_print, edges_measure = gen_edges(x,lx,ly)
        # generate a histogram of the density
        hist_rho = hist_density(x,y,edges_measure)
        # generate velocity field
        hist_vel = hist_velocity(x,y,edges_measure)

##################################################################

def gen_plots(rho_hist):
    """ genereate plots of the density"""
    for i in range(len(rho_hist)):
        plt.plot(rho_hist[i])
        plt.show()
        plt.close()
    return

##################################################################

def main():
    """ main function, called when script is started"""
    # analyse the data
    results = analyse_data()
    return

##################################################################

if __name__ == '__main__':
    main()
