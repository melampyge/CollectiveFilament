#!/usr/local/bin/python2.7

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import os

try:
    inputfile = sys.argv[1]
    ofname = sys.argv[2]
    dhead = int(sys.argv[3])
except:
    print 'Usage: ' + sys.argv[0] + '       dumpfile       output folder     d'
    exit()

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

def correct_pbc(z,lz):
    """ correct periodic boundary problem"""
    nz = len(z)
    for i in range(nz-1):
        # put z1 to initial box
        z1 = z[-1-i]
        z1 /= lz
        z1 = z1 - math.floor(z1)
        z1 *= lz
        # put z2 to initial box
        z2 = z[-2-i]
        z2 /= lz
        z2 = z2 - math.floor(z2)
        z2 *= lz
        # compute minimum distance
        dz = neigh_min(z1-z2,lz)
        z[-2-i] = z[-1-i] - dz
    return

##################################################################

def get_positions(dumpfile, dhead):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = open(dumpfile, 'r')
    ### generate lists for ex and ey
    xt = []
    yt = []
    xc = []
    yc = []
    xe = []
    ye = []
    xcom = []
    ycom = []
    # time
    time = []
    # first run variable
    counter = 0
    ### loop over all timesteps
    while True:
        # check whether end of file has been reached
        line = ifile.readline()
        if line == '':
            break
        # timestep
        line = ifile.readline()
        line = line.split()
        time.append(int(line[0]))
        # natoms
        ifile.readline()
        line = ifile.readline()
        line = line.split()
        natoms = int(line[0])
        if counter == 0:
            x = np.zeros((natoms))
            y = np.zeros((natoms))
        # box dimensions
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
        ifile.readline()
        # skip body headline
        ifile.readline()
        # read in positions
        for i in range(natoms):
            line = ifile.readline()
            line = line.split()
            aID = int(line[0])-1
            xs = (float(line[1]) + float(line[4]))*lx
            ys = (float(line[2]) + float(line[5]))*ly
            x[aID] = xs
            y[aID] = ys
        # correct pbcs, start in reverse order
        correct_pbc(x,lx)
        correct_pbc(y,ly)
        # first sphere atom
        xt.append(x[dhead-1])
        yt.append(y[dhead-1])
        # last sphere atom
        xc.append(x[0])
        yc.append(y[0])
        # last tail atom
        xe.append(x[-1])
        ye.append(y[-1])
        # center of mass
        xcom.append(np.average(x))
        ycom.append(np.average(y))
    
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    xt = np.array(xt, dtype = np.float64)
    yt = np.array(yt, dtype = np.float64)
    xc = np.array(xc, dtype = np.float64)
    yc = np.array(yc, dtype = np.float64)
    xe = np.array(xe, dtype = np.float64)
    ye = np.array(ye, dtype = np.float64)
    xcom = np.array(xcom, dtype = np.float64)
    ycom = np.array(ycom, dtype = np.float64)
    
    ### return data
    return time, xt, yt, xc, yc, xe, ye, xcom, ycom

##################################################################

def main():
    """ main function, called when the script is started"""
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    print 'reading end-to-end data'
    time, xt, yt, xc, yc, xe, ye, xcom, ycom = get_positions(inputfile, dhead)
   
    ### generate output folder
    os.system('mkdir ' + str(ofname))
       
    #### write tables
    # orientations
    ofile = open(ofname + '/positions.data', 'w')
    ofile.write('Current position of characteristic beads (t = tip, c = connection, e = end, com = center of mass)\n\n')
    ofile.write('timestep\txt\tyt\txc\tyc\txe\tye\txcom\tycom\n')
    for i in range(len(time)):
        ofile.write(str(time[i]) + '\t')
        ofile.write(str(xt[i]) + '\t' + str(yt[i]) + '\t')
        ofile.write(str(xc[i]) + '\t' + str(yc[i]) + '\t')
        ofile.write(str(xe[i]) + '\t' + str(ye[i]) + '\t')
        ofile.write(str(xcom[i]) + '\t' + str(ycom[i]) + '\n')
        
    ofile.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
