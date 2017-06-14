#!/usr/local/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import math

try:
    fname = sys.argv[1]
    nsteps = int(sys.argv[2])
    dcrit = float(sys.argv[3])
    lcrit = float(sys.argv[4])
    pcrit = float(sys.argv[5])
    npol = int(sys.argv[6])

except:
    print 'Usage: ' + sys.argv[0] + '      infilename       nsteps       dcrit       lcrit      circular fraction       polymer length'
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

############################################################################

def gen_linked_list(x,y,lx,ly,rc):
    """ generate a linked list"""
    # determine the number of cells in each direction
    nsegx = int(lx/rc)
    nsegy = int(ly/rc)

    # allocate head and llist
    ncells = nsegx*nsegy
    natoms = len(x)
    head = np.zeros((ncells), dtype = int)
    llist = np.zeros((natoms), dtype = int)

    # fill list and head
    for i in range(natoms):
        segx = int(x[i]/lx*nsegx)
        segy = int(y[i]/ly*nsegy)
        cell = segx*nsegy + segy
        llist[i] = head[cell]
        head[cell] = i
    
    return nsegx,nsegy,head,llist

##################################################################

def fill_neigh_matrix((neighs,llist,head,nsegx,nsegy,x,y,phi,mol,lx,ly,dcrit)):
    """ count the number of neighboring beads of all molecules"""
    for i in range(nsegx):
        for j in range(nsegy):
           # store header of current cell
            sv1 = head[i*nsegy + j]
            # loop over neighboring cells
            for a in range(3):
                i2 = (i-1+a)%nsegx
                for b in range(3):
                    j2 = (j-1+b)%nsegy
                    # store header of neighbor cell
                    sv2 = head[i2*nsegy + j2]
                    
                    # restore head values at for each new cell
                    val1 = sv1
                    val2 = sv2
                    while val1 != 0:
                        x1 = x[val1]/lx
                        y1 = y[val1]/ly
                        p1 = phi[val1]
                        mol1 = mol[val1]
                        while val2 != 0:
                            if val1 != val2:
                                x2 = x[val2]/lx
                                y2 = y[val2]/ly
                                p2 = phi[val2]
                                mol2 = mol[val2]
                                if mol1!= mol2:
                                    dx = x2-x1
                                    dx = dx - math.floor(dx + 0.5)
                                    dx = dx*lx
                                    dy = y2-y1
                                    dy = dy - math.floor(dy + 0.5)
                                    dy = dy*ly
                                    dphi = p1 - p2
                                    if dx**2 + dy**2 < dcrit**2 and min([dphi**2, (dphi+2*np.pi)**2, (dphi-2*np.pi)**2])  < (pcrit*2*np.pi)**2:
                                        neighs[mol1,mol2] += 1
                                        neighs[mol2,mol1] += 1
                            val2 = llist[val2]
                        val1 = llist[val1]
                        val2 = sv2
    return

##################################################################

def recursion(neighs,cl,i,nmol,ncrit):
    """ recursiely find clusters"""
    for j in range(nmol):
        if neighs[i,j] > ncrit:
                if cl[j] == 0:
                    cl[j] = cl[i]
                    recursion(neighs,cl,j,nmol,ncrit)

##################################################################

def cluster_search(neighs,cl):
    """ recursively search for clusters"""
    # define critical overlap
    ncrit = lcrit*npol
    # number of molecules
    nmol = len(neighs)
    # loop over all colums of the matrix:
    for i in range(nmol):
        # assign a cluster id to the current molecule
        if cl[i] != 0:
            continue
        cl[i] = max(cl) + 1
        recursion(neighs,cl,i,nmol,ncrit)

    return
                    


##################################################################

def find_clusters(x,y,phi,mol,lx,ly):
    """ search for clusters; two molecules are defined as being
        part of the same cluster if lcrit of their body length is
        within a distance of dcrit and when the difference in orientation
        is less than phi"""
    # allocate clusters array
    natoms = len(x)
    nmol = natoms/npol
    cl = np.zeros((nmol), dtype = int)
    clusters = np.zeros((natoms), dtype = int)
    # generate a linked list
    print 'Generating linked list'
    nsegx, nsegy, head, llist = gen_linked_list(x,y,lx,ly,dcrit)
    # allocate array to store neighbors
    print ' filling neighbor matrix'
    nmol = natoms/npol
    neighs = np.zeros((nmol,nmol), dtype = int)
    fill_neigh_matrix((neighs,llist,head,nsegx,nsegy,x,y,phi,mol,lx,ly,dcrit))
    # recursive search for clusters in neighbor matrix
    print ' recursive neighbor search'
    cluster_search(neighs,cl)
    # fill cluster results to per atom array
    k = 0
    for i in range(nmol):
        for j in range(npol):
            clusters[k] = cl[i]
            k = k + 1
    return clusters

##################################################################

def gen_cluster_plot(x,y,clusters):
    """ generate a plot of the atoms color coded with the cluster
        to which they belong"""
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.scatter(x,y,s=5,c=clusters, linewidths = 0)
    ax.axis('equal')
    plt.show()
    plt.close()
    return

##################################################################

def cluster_analysis():
    """ perform cluster analysis"""
    # open file for reading
    ifile = open(fname)
    # loop over all steps
    for i in range(nsteps):
        # read in the coordinates and mol
        x,y,mol,lx,ly = read_coordinates(ifile)
        # compute orientation of each bead
        phi = compute_orientation(x,y,lx,ly)
        # find clusters
        clusters = find_clusters(x,y,phi,mol,lx,ly)
        # generate a plot of the clusters
        gen_cluster_plot(x,y,clusters)
    ifile.close()
    # still need to decide what exactly will be returned
    return 0

##################################################################

def main():
    """ main function, called when script is started"""
    # start cluster analysis
    clusters = cluster_analysis()
    return

##################################################################

if __name__ == '__main__':
    main()
