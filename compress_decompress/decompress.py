#!/usr/local/bin/python

import numpy as np
import sys, math, os
import codecs

try:
    fname = sys.argv[1]
    ofname = sys.argv[2]
    hname = sys.argv[3]
except:
    print 'Usage: ' + sys.argv[0] + '      infilename         outfilename       header filename'
    exit()
    
##################################################################

def decompress_data():
    """ decompress unicode representation to dump file"""
    # read in the header
    hfile = open(hname)
    line = hfile.readline()
    line = line.split()
    nsnap = int(line[3])
    natoms = int(line[4])
    # start to uncompress the data
    ifile = codecs.open(fname, 'r', 'UTF-8')
    ofile = open(ofname, 'w')
    for i in range(nsnap):
	# print stats to screen
	print 'current snapshot / all snapshots:', i, nsnap
        # write down the header
        for j in range(9):
            line = hfile.readline()
            ofile.write(line)
        # read in the coords
        for j in range(natoms):
            b1 = ifile.read(1)
            b2 = ifile.read(1)
            b3 = ifile.read(1)
            b4 = ifile.read(1)

            x1 = ord(b1)
            x2 = ord(b2)
            y1 = ord(b3)
            y2 = ord(b4)

            x = float(x1*256 + x2)/256**2
            y = float(y1*256 + y2)/256**2
            
            ofile.write(str(j+1) + '\t1\t' + str(x) + '\t' + str(y) + '\t0.5\n')

    ifile.close()
    ofile.close()
    hfile.close()
    return

##################################################################

def main():
    """ main function, called when script is started"""
    # decompress the data
    decompress_data()
    return

##################################################################

if __name__ == '__main__':
    main()
