#!/usr/local/bin/python2.7
### read in xyz atom coordinates
### compute radius of gyration

import numpy as np
import sys, os, math
import matplotlib.pyplot as plt
from scipy import signal

try:
    fname = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '       filename'
    exit()



############################################################################

def read_data(ifile):
    """ read in coordinates of a single snapshot"""
    line = ifile.readline()
    line = line.split()
    try:
        natoms = int(line[0])
    except:
        return 0,0,0,1
    ifile.readline()
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    z = np.zeros((natoms))
    for i in range(natoms):
        line = ifile.readline()
        line = line.split()
        x[i] = float(line[1])
        y[i] = float(line[2])
        z[i] = float(line[3])
    return x,y,z,0

############################################################################

def compute_rgyr(x,y,z):
    """ compute the radius of gyration"""
    natoms = len(x)
    xav = np.average(x)
    yav = np.average(y)
    zav = np.average(z)
    rgyr = np.sqrt(np.average((x-xav)**2 + (y-yav)**2 + (z-zav)**2))
    return rgyr

############################################################################

def main():
    """ main function, called when the scrit is started"""
    # allocate rgyr list
    rgyr = []
    # open the xyz file
    ifile = open(fname, 'r')
    # compute rgyr till no further frames in file
    while 1:
        # try to read in the new snapshot
        x,y,z,endflag = read_data(ifile)
        if endflag:
            break
        # compute rgyr
        rgyr.append(compute_rgyr(x,y,z))
    # transform to numpy array
    rgyr = np.array(rgyr)
    # compute average
    rav = np.average(rgyr)
    # compute histogram
    rhist = np.histogram(rgyr, bins = 500, range = (3,30), density = True)
    # compute the spectral density using the Welch and Periodogram method
    freq_w, spec_w = signal.welch(rgyr)
    freq_q, spec_q = signal.periodogram(rgyr)
    # output rgyr
    ofile = open('rgyr.data', 'w')
    ofile.write('# Radius of gyration\n#\n')
    ofile.write('# Average: ' + str(rav) + '\n#\n')
    ofile.write('# frame\trgyr\n')
    for i in range(len(rgyr)):
        ofile.write(str(i) + '\t' + str(rgyr[i]) + '\n')
    ofile.close()
    # output histogram
    ofile = open('rgyr_hist.data', 'w')
    ofile.write('# histogram of the radius of gyration\n\n')
    ofile.write('# rgyr\tp\n')
    for i in range(len(rhist[0])):
        ofile.write(str(rhist[1][i]*0.5 + rhist[1][i+1]*0.5) + '\t' + str(rhist[0][i]) + '\n')
    ofile.close()
    # output power spectras
    ofile = open('rgyr_spectra_welch.data', 'w')
    ofile.write('# rgyr power spectra from the Welche method\n#\n')
    ofile.write('# omega\tP\n')
    for i in range(len(freq_w)):
        ofile.write(str(freq_w[i]) + '\t' + str(spec_w[i]) + '\n')
    ofile.close()
    ofile = open('rgyr_spectra_periodogram.data', 'w')
    ofile.write('# rgyr power spectra from the periodogram method\n#\n')
    ofile.write('# omega\tP\n')
    for i in range(len(freq_q)):
        ofile.write(str(freq_q[i]) + '\t' + str(spec_q[i]) + '\n')
    ofile.close()
    
    return

############################################################################

if __name__ == '__main__':
    main()
