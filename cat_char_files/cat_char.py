#!/usr/users/iff_th2/isele/Applications/Anaconda/bin/python2.7

import numpy as np
import codecs
import os

################################################################
# collection of functions to read in the data of a char file
################################################################

def read_first(hfile):
    line = hfile.readline()
    line = line.split()
    nsteps = int(line[3])
    natoms = int(line[4])
    return natoms, nsteps

################################################################
  
def skip_snapshots(hfile, ifile, nskip):
    """ skip some snapshots"""
    if nskip < 1:
        return
    # get number of atoms from the first header
    for i in range(3):
        hfile.readline()
    line = hfile.readline()
    line = line.split()
    natoms = int(line[0])
    for i in range(5):
        hfile.readline()
    # skip the remaining header lines
    for i in range(9*(nskip-1)):
        hfile.readline()
    # skip the body
    for i in range(nskip):
        ifile.reader.read(size = natoms*4, chars = natoms*4)
    return

################################################################

def main():
    """ main function, called when the script is started"""
    # generate a list with filenames
    os.system('ls out*.char > tmp.char')
    ifile = open('tmp.char')
    charfiles = []
    for line in ifile:
        line = line.split()
        charfiles.append(line[0])
    # check the number of filenames, leave if nfiles = 1
    nfiles = len(charfiles)
    if nfiles == 1:
        print 'NO CONCENTATION NECESSARY'
        exit()
    # compute the number of total steps
    print 'Determining the number of total steps'
    natoms_all = 0
    nsteps_all = 0
    for i in range(nfiles):
        hfile = open('out' + str(i+1) + '.header')
        natoms, nsteps = read_first(hfile)
        hfile.close()
        natoms_all = natoms
        nsteps_all = nsteps_all + nsteps
    nsteps_all = nsteps_all - nfiles + 1
    # removing previous files if there are any
    os.system('rm out.header.all out.char.all')
    # cat header files to a master file
    print 'Combining header files'
    ofile = open('out.header.all', 'w')
    ofile.write('nsnapshots, natoms = ' + str(nsteps_all) + ' ' + str(natoms_all) + '\n')
    ofile.close()
    for i in range(nfiles):
        if i == 0:
            os.system('tail -n +2 out1.header >> out.header.all')
        else:
            os.system('tail -n +11 out' + str(i+1) + '.header >> out.header.all')
    # cat char files to a master file
    print 'Combining char files'
    ofile = codecs.open('out.char.all', 'w', 'UTF-8')
    for i in range(nfiles):
        hfile = open('out' + str(i+1) + '.header')
        cfile = codecs.open('out' + str(i+1) + '.char', 'r', 'UTF-8')
        natoms, nsteps = read_first(hfile)
        ### skip first snapshot if i > 0
        if i > 0:
            skip_snapshots(hfile, cfile, 1)
            nsteps = nsteps - 1
        ### copy the rest of the snapshot
        for j in range(nsteps):
            print 'i/nfiles; j/nsteps', i,nfiles,j,nsteps
            b = cfile.reader.read(natoms*4,natoms*4)
            ofile.write(b)
        hfile.close()
        cfile.close()
    ofile.close()
    return

################################################################

if __name__ == '__main__':
    main()
