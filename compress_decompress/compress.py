#!/usr/local/bin/python

import numpy as np
import sys, math, os
import codecs

try:
    fname = sys.argv[1]
    ofname = sys.argv[2]
    hname = sys.argv[3]
    nsnap = int(sys.argv[4])
except:
    try:
        fname = sys.argv[1]
        ofname = sys.argv[2]
        hname = sys.argv[3]
        nsnap = 0
    except:
        print 'Usage: ' + sys.argv[0] + '      infilename         outfilename     header filename      nsnap (optional)'
        exit()


##################################################################

def compress_data(nsnap):
    """ compress the data into a two byte file"""
    # determine the number of atoms
    ifile = open(fname)
    ifile.readline()
    ifile.readline()
    ifile.readline()
    line = ifile.readline()
    line = line.split()
    natoms = int(line[0])
    ifile.close()
    # determine the number of snapshots
    if nsnap == 0:
        os.system('wc -l ' + fname + ' > tmp.txt')
        ifile = open('tmp.txt')
        line = ifile.readline()
        line = line.split()
        nlines = int(line[0])
        nsnap = nlines/(natoms+9)
        ifile.close()
        os.system('rm tmp.txt')
    # allocate an array to store the x and y coords and modified coords
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    xi = np.zeros((natoms), dtype = int)
    yi = np.zeros((natoms), dtype = int)
    x1 = np.zeros((natoms), dtype = int)
    x2 = np.zeros((natoms), dtype = int)
    y1 = np.zeros((natoms), dtype = int)
    y2 = np.zeros((natoms), dtype = int)
    # open files, start to compress the data
    ifile = open(fname)
    hfile = open(hname, 'w')
    hfile.write('nsnapshots, natoms = ' + str(nsnap) + ' ' + str(natoms) + '\n')
    ofile = codecs.open(ofname, 'w', 'UTF-8')
    for i in range(nsnap):
	print 'current snapshot / all snapshots:',i, nsnap
        # copy the header
        for j in range(9):
            line = ifile.readline()
            hfile.write(line)
        # read in the coords
        for j in range(natoms):
            line = ifile.readline()
            line = line.split()
            aID = int(line[0]) - 1
            xs = float(line[2])
            ys = float(line[3])
            xs = xs - math.floor(xs)
            ys = ys - math.floor(ys)
            x[aID] = xs
            y[aID] = ys
        # transfrom the coords to integers and bytes to store
        for j in range(natoms):
            xi[j] = x[j]*256**2
            yi[j] = y[j]*256**2
            x1[j] = xi[j]/256
            y1[j] = yi[j]/256
            x2[j] = xi[j]%256
            y2[j] = yi[j]%256
        # write down the bytes as utf8 code to a file
        for j in range(natoms):
            u1 = unichr(x1[j])
            u2 = unichr(x2[j])
            u3 = unichr(y1[j])
            u4 = unichr(y2[j])
            ofile.write(u1)
            ofile.write(u2)
            ofile.write(u3)
            ofile.write(u4)
    # close files and return
    ifile.close()
    ofile.close()
    hfile.close()
    return

##################################################################

def main():
    """ main function, called when script is started"""
    # compress the data
    compress_data(nsnap)
    return

##################################################################

if __name__ == '__main__':
    main()
