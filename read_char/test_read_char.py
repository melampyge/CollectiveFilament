#!/usr/users/iff_th2/isele/Applications/Anaconda/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import math
import codecs
import read_char

try:
    fname = sys.argv[1]
    hname = sys.argv[2]
except:
    print 'Usage: ' + sys.argv[0] + '      char file        header file'
    exit()
    

##################################################################

def write(x,y,f):
    ofile = open(f, 'w')
    for i in range(len(x)):
        ofile.write(str(x[i]) + '\t' + str(y[i]) + '\n')
    ofile.close()

##################################################################

def main():
    """ test the functions in read_char"""
    ### define nskip
    nskip = 100

    ### open input files for reading
    ifile = codecs.open(fname, 'r', 'UTF-8')
    hfile = open(hname, 'r')
    natoms, nsteps = read_char.read_first(hfile)
    ### read in nskip snapshots
    for i in range(nskip):
        xi,yi,lx,ly,tstep,natoms = read_char.read_snapshot(hfile, ifile)
    ### read in snapshot nskip+1 and print results to file
    xi,yi,lx,ly,tstep,natoms = read_char.read_snapshot(hfile, ifile)
    write(xi,yi,'no_skip.data')
    plt.plot(xi,yi,ls = '', marker = 'o', markersize = 1)
    plt.show()
    plt.close()
    ### close files
    ifile.close()
    hfile.close()

    ### open input files for reading
    ifile = codecs.open(fname, 'r', 'UTF-8')
    hfile = open(hname, 'r')
    natoms, nsteps = read_char.read_first(hfile)
    ### skip nskip snapshots
    read_char.skip_snapshots(hfile,ifile,nskip)
    ### read in snapshot nskip+1 and print results to file
    xi,yi,lx,ly,tstep,natoms = read_char.read_snapshot(hfile, ifile)
    write(xi,yi,'skip.data')
    ### close the files
    ifile.close()
    hfile.close()
    return

##################################################################

if __name__ == '__main__':
    main()
