#!/usr/local/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import os

try:
    fname = sys.argv[1]
    npol = int(sys.argv[2])

except:
    print 'Usage: ' + sys.argv[0] + '      infilename       polymer length'
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

def bond_overlap(x,y,lx,ly,rc,nsegx,nsegy,head,llist):
    """ loop over the linked list and search for bond overlap"""
    counter = 0
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
                        x1 = x[val1]
                        y1 = y[val1]
                        x2 = x[val1-1]
                        y2 = y[val1-1]
                        dx = nearest_neighbor(x1,x2,lx)
                        x2 = x1 - dx
                        dy = nearest_neighbor(y1,y2,ly)
                        y2 = y1 - dy
                        while val2 != 0:
                            if math.fabs(val1 -val2) > 1.5:
                                x3 = x[val2]
                                y3 = y[val2]
                                x4 = x[val2-1]
                                y4 = y[val2-1]
                                dx = nearest_neighbor(x3,x4,lx)
                                x4 = x3 - dx
                                dy = nearest_neighbor(y3,y4,ly)
                                y4 = y3 - dy
                                ### check for overlap
                                # check whether one of the atoms is a startomg atom
                                checkflag = 1
                                if val1 % npol == 0 or val2 % npol == 0:
                                    checkflag = 0
                                # loop over all bonds and check whether they

                                #  cross, if yes, print out a warning message
                                #  and increase the counter
                                if checkflag == 1:
                                    # compute intersection
                                    d = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
                                    xi = ((x3-x4)*(x1*y2-y1*x2)-(x1-x2)*(x3*y4-y3*x4))/d
                                    yi = ((y3-y4)*(x1*y2-y1*x2)-(y1-y2)*(x3*y4-y3*x4))/d
                                    if xi > min(x1,x2) and xi < max(x1,x2) and xi > min(x3,x4) and xi < max(x3,x4):
                                        counter = counter + 1
                                        print 'Found intersection for atoms ', val1, val2,' close to x1, y1 =', x1,y1
 

                            val2 = llist[val2]
                        val1 = llist[val1]
                        val2 = sv2
    print 'Found ' + str(counter/2) + ' overlaping bonds.'
    return

##################################################################

def search_bond_overlap():
    """ search for bond crossing"""
    # read in the data
    ifile = open(fname)
    # read coordinates
    x,y,nmol,lx,ly = read_coordinates(ifile)
    # generate a linked list
    rc = 3.0
    nsegx, nsegy, head, llist = gen_linked_list(x,y,lx,ly,rc)
    # find overlap
    bond_overlap(x,y,lx,ly,rc,nsegx,nsegy,head,llist)

##################################################################

def main():
    """ main function, called when script is started"""
    # analyse the clusters of the simulation
    search_bond_overlap()
    return

##################################################################

if __name__ == '__main__':
    main()
