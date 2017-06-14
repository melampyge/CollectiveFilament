#!/usr/local/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import os

try:
    npol = int(sys.argv[1])
except:
    print 'Usage: ' + sys.argv[0] + '        polymer length'
    exit()
 
##################################################################

def gen_coords():
    """ generate some example clusters"""
    x = []
    y = []
    ### two parallel rods, L = npol
    x0 = 0.
    y0 = 0.
    for i in range(npol):
        x.append(x0 + i)
        y.append(y0 + 0)
    for i in range(npol):
        x.append(x0 + i)
        y.append(y0 + 1)
    ### two parrallel rods, L = npol/2
    x0 = 0
    y0 = 10
    for i in range(npol):
        x.append(x0 + 0.5*i)
        y.append(y0 + 0)
    for i in range(npol):
        x.append(x0 + 0.5*i)
        y.append(y0 + 1)
    ### two by two parrallel rods, Li = npol/2
    x0 = 0
    y0 = 20
    for i in range(npol*2):
        x.append(x0 + 0.5*i)
        y.append(y0 + 0)
    for i in range(npol*2):
        x.append(x0 - 0.25*npol + 0.5*i)
        y.append(y0 +1)


    ### two small circles, 2 rods each circle
    # define inner ring
    ui = 2*npol
    ri = 0.5*ui/np.pi
    ti = np.linspace(0,2*np.pi,num=2*npol,endpoint = False)
    cosi = np.cos(ti)*ri
    sini = np.sin(ti)*ri
    # define outer ring
    ro = ri + 1
    to = np.linspace(0,2*np.pi,num=2*npol,endpoint = False) + np.pi/2
    coso = np.cos(to)*ro
    sino = np.sin(to)*ro
    x0 = ro
    y0 = y0 + ro + 10
    # add values
    for i in range(2*npol):
        x.append(x0 + cosi[i])
        y.append(y0 + sini[i])
    for i in range(2*npol):
        x.append(x0 + coso[i])
        y.append(y0 + sino[i])
    # updated y0
    y0 += ro
    ### two small circles, 4 rods each circle
     # define inner ring
    ui = 2*npol
    ri = 0.5*ui/np.pi
    ti = np.linspace(0,2*np.pi,num=4*npol,endpoint = False)
    cosi = np.cos(ti)*ri
    sini = np.sin(ti)*ri
    # define outer ring
    ro = ri + 1
    to = np.linspace(0,2*np.pi,num=4*npol,endpoint = False) + np.pi/2
    coso = np.cos(to)*ro
    sino = np.sin(to)*ro
    # add values
    x0 = ro
    y0 = y0 + ro + 10
    for i in range(4*npol):
        x.append(x0 + cosi[i])
        y.append(y0 + sini[i])
    for i in range(4*npol):
        x.append(x0 + coso[i])
        y.append(y0 + sino[i])
    # updated y0
    y0 += ro

    ### two larger circles, 2 rods each circle
    # define inner ring
    ui = 4*npol
    ri = 0.5*ui/np.pi
    ti1 = np.linspace(0,0.5*np.pi,num=npol,endpoint = False)
    cosi1 = np.cos(ti1)*ri
    sini1 = np.sin(ti1)*ri
    ti2 = np.linspace(0,0.5*np.pi,num=npol,endpoint = False) + np.pi
    cosi2 = np.cos(ti2)*ri
    sini2 = np.sin(ti2)*ri 
    # define outer ring  
    ro = ri + 1
    to1 = np.linspace(0,0.5*np.pi,num=npol,endpoint = False) + 0.5*np.pi
    coso1 = np.cos(to1)*ro
    sino1 = np.sin(to1)*ro
    to2 = np.linspace(0,0.5*np.pi,num=npol,endpoint = False) + 1.5*np.pi
    coso2 = np.cos(to2)*ro
    sino2 = np.sin(to2)*ro
    # add values
    x0 = ro
    y0 = y0 + ro + 10
    for i in range(npol):
        x.append(x0 + cosi1[i])
        y.append(y0 + sini1[i])
    for i in range(npol):
        x.append(x0 + cosi2[i])
        y.append(y0 + sini2[i])
    for i in range(npol):
        x.append(x0 + coso1[i])
        y.append(y0 + sino1[i])
    for i in range(npol):
        x.append(x0 + coso2[i])
        y.append(y0 + sino2[i])

    ### transform lists to arrays and return
    x = np.array(x)
    y = np.array(y)
    return x,y


##################################################################

def write_dump(x,y):
    """ write LAMMPS dumpfile-like output"""
    # define header values
    natoms = len(x)
    xlo = min(x) - 10
    xhi = max(x) + 10
    lx = xhi - xlo
    ylo = min(y) - 10
    yhi = max(y) + 10
    ly = yhi - ylo
    # open file for writing
    ofile = open('test.dump', 'w')
    # write header
    ofile.write('ITEM: TIMESTEP\n')
    ofile.write('0\n')
    ofile.write('ITEM: NUMBER OF ATOMS\n')
    ofile.write(str(natoms) + '\n')
    ofile.write('ITEM: BOX BOUNDS pp pp pp\n')
    ofile.write(str(xlo) + ' ' + str(xhi) + '\n')
    ofile.write(str(ylo) + ' ' + str(yhi) + '\n')
    ofile.write('-5 5\n')
    ofile.write('ITEM: ATOMS id type xs ys zs\n')
    # write body
    for i in range(natoms):
        ofile.write(str(i+1) + ' 1 ' + str((x[i]-xlo)/lx) + ' ' + str((y[i]-ylo)/ly) + ' 0.5\n')
    # close file
    ofile.close()
    return

##################################################################

def main():
    """ main function, called when script is started"""
    # generate an example to test the cluster analysis
    x,y = gen_coords()
    # write LAMMPS dump file
    write_dump(x,y)
    return

##################################################################

if __name__ == '__main__':
    main()
