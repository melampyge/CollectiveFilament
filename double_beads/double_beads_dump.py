#!/usr/local/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import math


try:
    ifname = sys.argv[1]
    ofname = sys.argv[2]
    nfil = int(sys.argv[3])
    lscale = float(sys.argv[4])
    
except:
    print 'Usage: ' + sys.argv[0] + '      input dump file       output data file       nfil    lscale'
    exit()

##################################################################

def read_dump(ifname):
    """ read information from dump file"""
    ifile = open(ifname)
    ### read header
    # timestep
    ifile.readline()
    line = ifile.readline()
    line = line.split()
    tstep = int(line[0])
    # natoms
    ifile.readline()
    line = ifile.readline()
    line = line.split()
    natoms = int(line[0])
    # box dimensions
    ifile.readline()
    line = ifile.readline()
    line = line.split()
    xlo = float(line[0])
    xhi = float(line[1])
    line = ifile.readline()
    line = line.split()
    ylo = float(line[0])
    yhi = float(line[1])
    line = ifile.readline()
    line = line.split()
    zlo = float(line[0])
    zhi = float(line[1])
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo
    # last header line
    ifile.readline()
    ### allocate memory
    xs = np.zeros((natoms))
    ys = np.zeros((natoms))
    ### read the body
    for i in range(natoms):
        line = ifile.readline()
        line = line.split()
        aID = int(line[0]) - 1
        xi = float(line[2])
        yi = float(line[3])
        xi = xi - math.floor(xi)
        yi = yi - math.floor(yi)
        xs[aID] = xi
        ys[aID] = yi
    ### close file and return data
    ifile.close()
    return tstep, natoms, lx, ly, xs, ys

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

def correct_pbcs(xs, ys, nfil, natoms):
    """ correct periodic boundary conditions"""
    ### allocate distance arrays
    dx = np.zeros((nfil - 1))
    dy = np.zeros((nfil - 1))
    ### loop over all molecules
    nmol = natoms / nfil
    for i in range (nmol):
        ### compute the minimum distances in x and y direction
        for j in range(nfil-1):
            d = xs[i*nfil + j+1] - xs[i*nfil + j]
            dx[j] = neigh_min(d,1.0)
            d = ys[i*nfil + j+1] - ys[i*nfil + j]
            dy[j] = neigh_min(d,1.0)
        for j in range(nfil-1):
            xs[i*nfil + j+1] = xs[i*nfil + j] + dx[j]
            ys[i*nfil + j+1] = ys[i*nfil + j] + dy[j]
    return

##################################################################

def add_beads(xs,ys,nfil,natoms):
    """ add new beads in between"""
    ### determine array sizes
    nmol = natoms/nfil
    nfild = 2*nfil - 1
    natomsd = nmol*nfild
    ### allocate output arrays
    xsd = np.zeros((natomsd))
    ysd = np.zeros((natomsd))
    ### loop over all molecules
    k = 0
    for i in range(nmol):
        # add first bead
        xsd[k] = xs[i*nfil]
        ysd[k] = ys[i*nfil]
        k = k + 1
        # add remaining beads
        for j in range(nfil - 1):
            # new intermediate bead
            xsd[k] = 0.5*(xs[i*nfil + j] + xs[i*nfil + j + 1])
            ysd[k] = 0.5*(ys[i*nfil + j] + ys[i*nfil + j + 1])
            k = k + 1
            # copy existing bead
            xsd[k] = xs[i*nfil + j + 1]
            ysd[k] = ys[i*nfil + j + 1]
            k = k + 1
    return xsd, ysd, natomsd, nfild

##################################################################

def add_per_atom_info(natomsd, nfild):
    """ add per atom properties"""
    tpe = np.ones((natomsd), dtype = int)
    molid = np.ones((natomsd), dtype = int)
    for i in range(natomsd):
        molid[i] = i / nfild + 1
    return tpe, molid

##################################################################

def add_bond_information(natomsd,nfild):
    """ add bond information"""
    ### allocate data
    nmol = natomsd/nfild
    nbonds = nmol*(nfild-1)
    bid = np.zeros((nbonds), dtype = int)
    btpe = np.zeros((nbonds), dtype = int)
    b1 = np.zeros((nbonds), dtype = int)
    b2 = np.zeros((nbonds), dtype = int)
    ### add new bonds
    k = 0
    for i in range(nmol):
        for j in range(nfild-1):
            bid[k] = k + 1
            btpe[k] = 1
            b1[k] = i*nfild + j + 1
            b2[k] = i*nfild + j + 2
            k = k + 1
    return bid, btpe, b1, b2

##################################################################

def add_angle_information(natomsd, nfild):
    """ add bond information"""
    ### allocate data
    nmol = natomsd/nfild
    nangles = nmol*(nfild-2)
    aid = np.zeros((nangles), dtype = int)
    atpe = np.zeros((nangles), dtype = int)
    a1 = np.zeros((nangles), dtype = int)
    a2 = np.zeros((nangles), dtype = int)
    a3 = np.zeros((nangles), dtype = int)
    ### add new bonds
    k = 0
    for i in range(nmol):
        for j in range(nfild-2):
            aid[k] = k + 1
            atpe[k] = 1
            a1[k] = i*nfild + j + 1
            a2[k] = i*nfild + j + 2
            a3[k] = i*nfild + j + 3
            k = k + 1
    return aid, atpe, a1, a2, a3

##################################################################

def write_xyz(xs,ys):
    """ write xyz"""
    natoms = len(xs)
    ofile = open('test.xyz', 'w')
    ofile.write(str(natoms) + '\n')
    ofile.write('This is a comment line\n')
    for i in range(natoms):
        ofile.write('C ' + str(xs[i]) + ' ' + str(ys[i]) + ' 0.0\n')
    ofile.close()

##################################################################

def write_data(natoms, lx,ly, tpe, molid, x, y, bid, btpe, b1, b2, aid, atpe, a1, a2, a3,ofname):
    """ write LAMMPS type data file"""
    ofile = open(ofname, 'w')
    ### write down header information
    ofile.write('LAMMPS data file filaments in 2D\n\n')
    ofile.write(str(natoms) + ' atoms\n')
    ofile.write('1 atom types\n')
    ofile.write(str(max(bid)) + ' bonds\n')
    ofile.write('1 bond types\n')
    ofile.write(str(max(aid)) + ' angles\n')
    ofile.write('1 angle types\n\n')
    ofile.write('0.0 ' + str(lx) + ' xlo xhi\n')
    ofile.write('0.0 ' + str(ly) + ' ylo yhi\n')
    ofile.write('-2.5 2.5 zlo zhi\n\n')
    ofile.write('Masses\n\n')
    ofile.write('1 1\n\n')
    ### Atoms section
    ofile.write('Atoms\n\n')
    for i in range(natoms):
        ofile.write(str(i+1) + ' ' + str(molid[i]) + ' ' + str(tpe[i]) + ' ' + str(x[i]) + ' ' + str(y[i]) + ' 0.0\n')
    ofile.write('\n')
    ### Bonds section
    ofile.write('Bonds\n\n')
    for i in range(len(bid)):
        ofile.write(str(bid[i]) + ' ' + str(btpe[i]) + ' ' + str(b1[i]) + ' ' + str(b2[i]) + '\n')
    ofile.write('\n')
    ### Angles section
    ofile.write('Angles\n\n')
    for i in range(len(aid)):
        ofile.write(str(aid[i]) + ' ' + str(atpe[i]) + ' ' + str(a1[i]) + ' ' + str(a2[i]) + ' ' + str(a3[i]) + '\n')
    ofile.write('\n')
    ofile.close()
    return
        
##################################################################

def main():
    """ main function, called when script is started"""
    ### read dump
    timestep, natoms, lx, ly, xs, ys = read_dump(ifname)
    ### correct pbcs for all molecules
    correct_pbcs(xs, ys, nfil, natoms)
    ### add new beads
    xsd, ysd, natomsd, nfild = add_beads(xs,ys,nfil,natoms)
    ### add new per/atom information
    tpe, molid = add_per_atom_info(natomsd, nfild)
    ### add new bond information
    bid, btpe, b1, b2 = add_bond_information(natomsd, nfild)
    ### add new angle information
    aid, atpe, a1, a2, a3 = add_angle_information(natomsd, nfild)
    ### scale xs and ys
    lx *= lscale
    ly *= lscale
    xsd *= lx
    ysd *= ly
    ### write results to xyz file
    write_xyz(xsd,ysd)
    ### write results to lammps data file
    write_data(natomsd, lx,ly, tpe, molid, xsd, ysd, bid, btpe, b1, b2, aid, atpe, a1, a2, a3,ofname)
    return

##################################################################

if __name__ == '__main__':
    main()
