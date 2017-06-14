#!/usr/local/bin/python2.7

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import os
import performance_toolsWrapper

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

def get_ee(dumpfile, dhead):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = open(dumpfile, 'r')
    ### generate lists for ex and ey
    exrr = []
    eyrr = []
    extt = []
    eytt = []
    exrt = []
    eyrt = []
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
        x1 = x[dhead-1]
        y1 = y[dhead-1]
        # last sphere atom
        x2 = x[0]
        y2 = y[0]
        # last tail atom
        x3 = x[-1]
        y3 = y[-1]
        # fill arrays
        exrr.append(x1-x2)
        eyrr.append(y1-y2)
        extt.append(x2-x3)
        eytt.append(y2-y3)
        exrt.append(x1-x3)
        eyrt.append(y1-y3)
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    exrr = np.array(exrr, dtype = np.float64)
    eyrr = np.array(eyrr, dtype = np.float64)
    extt = np.array(extt, dtype = np.float64)
    eytt = np.array(eytt, dtype = np.float64)
    exrt = np.array(exrt, dtype = np.float64)
    eyrt = np.array(eyrt, dtype = np.float64)
    
    ### return data
    return time, exrr, eyrr, extt, eytt, exrt, eyrt

##################################################################

def compute_angle(ex,ey):
    """ compute phi"""
    ### allocate phi array
    nsteps = len(ex)
    phi = np.zeros((nsteps), dtype = np.float64)
    dp = np.zeros((nsteps), dtype = np.float64)
    ### compute angle and correct angle on the fly
    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.compute_angle(ex,ey,phi,dp,nsteps)
    return phi

##################################################################

def main():
    """ main function, called when the script is started"""
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    print 'reading end-to-end data'
    time, exrr, eyrr, extt, eytt, exrt, eyrt = get_ee(inputfile, dhead)
    ### compute local angle phi based on overall rotation
    print 'computing orientation'
    phirr = compute_angle(exrr, eyrr)
    phitt = compute_angle(extt, eytt)
    phirt = compute_angle(exrt, eyrt)
   
    ### generate output folder
    os.system('mkdir ' + str(ofname))
       
    #### write tables
    # orientations
    ofile = open(ofname + '/orientations.data', 'w')
    ofile.write('Current orientation of the rod, tail, and entire swimmer\n\n')
    ofile.write('timestep\tphirr\tphitt\rphirt\n')
    for i in range(len(time)):
        ofile.write(str(time[i]) + '\t' + str(phirr[i]) + '\t' + str(phitt[i]) + '\t' + str(phirt[i]) + '\n')
    ofile.close()
        
##################################################################

if __name__ == '__main__':
    main()
    
