#!/usr/local/bin/python

import numpy as np
import sys, math, os
import codecs
import read_char

try:
    charfile = sys.argv[1]
    dumpfilename = sys.argv[2]
    headerfile = sys.argv[3]
    snap = int(sys.argv[4])
except:
    print 'Usage: ' + sys.argv[0] + '      charfilename         outfilename       header filename         frame to decompress (-1 == last snapshot)'
    exit()

##################################################################

def main(snap):
    """ main function, called when script is started"""
    ### open file for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    ### adjust snap
    if snap >= nsteps or snap == -1:
        snap = nsteps - 1
    ### skip the not required files
    read_char.skip_snapshots(hfile, ifile, snap)
    ### read the desired snapshot
    xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
    ### write a dump file
    ofile = open(dumpfilename, 'w')
    # write down the header
    ofile.write('ITEM: TIMESTEP\n')
    ofile.write(str(tstep) + '\n')
    ofile.write('ITEM: NUMBER OF ATOMS\n')
    ofile.write(str(natoms) + '\n')
    ofile.write('ITEM: BOX BOUNDS pp pp pp\n')
    ofile.write('0 ' + str(lx) + '\n')
    ofile.write('0 ' + str(ly) + '\n')
    ofile.write('-5 5\n')
    ofile.write('ITEM: ATOMS id type xs ys zs\n')
    # write down the body
    for i in range(natoms):
        ofile.write(str(i+1) + ' ') # ID
        ofile.write('1 ') # Atomtype
        ofile.write(str(xs[i]) + ' ') # xs
        ofile.write(str(ys[i]) + ' ') # ys
        ofile.write('0.5\n') # zs
    # close all files
    ofile.close()
    ifile.close()
    hfile.close()
    return

##################################################################

if __name__ == '__main__':
    main(snap)
