#!/usr/local/bin/python2.7

import numpy as np
import sys
import math
import os
try:
    inputfile = sys.argv[1]
    idstore = int(sys.argv[2])
except:
    print 'Usage: ' + sys.argv[0] + '       dumpfile        idstore'
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

def get_com(inputfile, idstore):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = open(inputfile, 'r')
    ### generate lists for coordinates x and y
    comx = []
    comy = []
    # time
    time = []
    # first run variable
    counter = 0
    ### loop over all timesteps
    while True:
        print counter
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
            x = np.zeros((natoms), dtype = np.float64)
            y = np.zeros((natoms), dtype = np.float64)
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
            #xs = float(line[1])*lx
            #ys = float(line[2])*ly
            x[aID] = xs
            y[aID] = ys
        # correct periodic boundaries
        correct_pbc(x, lx)
        correct_pbc(y, ly)
        # add coordinates to array
        #comx.append(np.average(x))
        #comy.append(np.average(y))
        comx.append(x[idstore])
        comy.append(y[idstore])
        counter = counter + 1
    # close the input file
    ifile.close()
    # transform lists to numpy arrays
    time = np.array(time, dtype = np.int32)
    comx = np.array(comx, dtype = np.float64)
    comy = np.array(comy, dtype = np.float64)
    ### return data
    return time, comx, comy

##################################################################

def compute_curvature(xall,yall):
    """ compute curvature"""
    ### allocate curvature array
    nsteps = len(xall)
    natoms = len(xall[0])
    curv = np.zeros((nsteps, natoms-2), dtype = np.float64)
    ### compute the local curvature
    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.compute_curvature(curv, xall, yall, nsteps, natoms)

    return curv


##################################################################

def main():
    """ main function, called when the script is started"""
    ### read in the coordinates for the entire trajectory, discard center of mass movement
    print '  reading in coordinates'
    t, x, y = get_com(inputfile, idstore)

    ### write information to file
    os.system('mkdir com_head')
    ofile = open('com_head/com_head.data', 'w')
    ofile.write('coordinates of the com\n\n')
    ofile.write('t\tx\ty\n')
    for i in range(len(t)):
        ofile.write(str(t[i]) + '\t' + str(x[i]) + '\t' + str(y[i]) + '\n')
    ofile.close()
            
##################################################################

if __name__ == '__main__':
    main()
    
