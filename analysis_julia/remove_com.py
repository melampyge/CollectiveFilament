#!/usr/users/iff_th2/isele/Applications/Anaconda/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10

@author: isele-holder

Remove the com from an xyz file
"""

import sys, os
import numpy as np

try:
    xyzfile = sys.argv[1]                    # file with coordinates
    outfile = sys.argv[2]                     # folder to store the results
    
except:
    print 'Usage: ' + sys.argv[0] + '      xyz-file     output file'
    exit()

###############################################################################

def remove_com(xyzfile, outfile):
    """ remove the com of an xyz file"""
    ifile = open(xyzfile, 'r')
    ofile = open(outfile, 'w')
    while True:
        ### read snapshot
        line = ifile.readline()
        if line == '':
            break
        line = line.split()
        natoms = int(line[0])
        line = ifile.readline()
        x = np.zeros((natoms))
        y = np.zeros((natoms))
        for i in range(natoms):
            line = ifile.readline()
            line = line.split()
            x[i] = float(line[1])
            y[i] = float(line[2])
        x -= np.average(x)
        y -= np.average(y)
        ### writre snapshot
        ofile.write(str(natoms) + '\n')
        ofile.write('This is a comment\n')
        for i in range(natoms):
            ofile.write('1 ' + str(x[i]) + ' ' + str(y[i]) + ' 0\n')
    ifile.close()
    ofile.close()
    return

###############################################################################
 
def main():
    """ analyse principal components"""
    remove_com(xyzfile, outfile)
    return

###############################################################################

if __name__ == '__main__':
    main()
