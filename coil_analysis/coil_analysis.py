#!/usr/local/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys, math

try:
    fname = sys.argv[1]
    nsteps = int(sys.argv[2])
except:
    print 'Usage: ' + sys.argv[0] + '      infilename         #(steps)'
    exit()
    
##################################################################

def read_snapshot(ifile):
    """ read in the coordinates from the xyz file"""
    # get the number of atoms
    line = ifile.readline()
    line = line.split()
    natoms = int(line[0])
    # allocate array to store data
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    # read in the timestep
    line = ifile.readline()
    line = line.split()
    t = int(line[-1])
    # read in the coordinates
    for i in range(natoms):
        line = ifile.readline()
        line = line.split()
        x[i] = float(line[-3])
        y[i] = float(line[-2])
    # return the data array
    return x,y,t

##################################################################

def get_distance(x,y):
    """ compute the distance between the first atom and all
        other atoms"""
    dist = (x-x[0])**2 + (y-y[0])**2
    return dist

##################################################################

def contact_atoms(dist):
    """ find contact atoms"""
    # search for the first minimum
    c1 = -1
    c2 = -1
    # search critical distance for evaluation
    dcrit = 3.0**2
    for i in range(1,len(dist)-1):
        # check for minimum
        if dist[i-1] > dist[i] and dist[i+1] > dist[i] and dist[i] < dcrit:
            if dist[i-1] < dist[i+1]:
                c1 = i-1
                c2 = i
            else:
                c1 = i
                c2 = i+1
            break
    return c1, c2

##################################################################

def write_results(t,c1):
    """ write results to file"""
    ofile = open('contact_atom.data', 'w')
    ofile.write('# atom against which the first atom pushes\n\n')
    ofile.write('t\taID\n')
    for i in range(len(t)):
        ofile.write(str(t[i]) + '\t' + str(c1[i]) + '\n')
    ofile.close
    return

##################################################################

def run_analysis():
    """ detect whether coils are existent and how compact they are"""
    # open the file for reading
    ifile = open(fname)
    c1 = []
    c2 = []
    t = []
    # loop over all snapshots
    for i in range(nsteps):
        # read data of a single snapshot
        x,y,ti = read_snapshot(ifile)
        # check the distance of the first atom to all other atoms
        dist = get_distance(x,y)
        # search contact atoms
        c1i,c2i = contact_atoms(dist)
        # append values to results arrays
        t.append(ti)
        c1.append(c1i)
        c2.append(c2i)
    # close the input file
    ifile.close()
    # write results for c1 to file
    write_results(t,c1)
    return

##################################################################

def main():
    """ main function, called when script is started"""
    # run the analysis
    run_analysis()
    return

##################################################################

if __name__ == '__main__':
    main()
