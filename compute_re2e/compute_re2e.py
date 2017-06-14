#!/usr/local/bin/python2.7
### read in xyz atom coordinates
### compute end-to-end distance

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

def compute_re2e(x,y,z):
    """ compute the radius of gyration"""
    x1 = x[0]
    y1 = y[0]
    z1 = z[0]
    x2 = x[-1]
    y2 = y[-1]
    z2 = z[-1]
    re2e = np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
    return re2e

############################################################################

def main():
    """ main function, called when the scrit is started"""
    # allocate re2e list
    re2e = []
    # open the xyz file
    ifile = open(fname, 'r')
    # compute re2e till no further frames in file
    while 1:
        # try to read in the new snapshot
        x,y,z,endflag = read_data(ifile)
        if endflag:
            break
        # compute re2e
        re2e.append(compute_re2e(x,y,z))
    # transform to numpy array
    re2e = np.array(re2e)
    # compute average
    rav = np.average(re2e)
    # compute histogram
    rhist = np.histogram(re2e, bins = 500, range = (0,50), density = True)
    # compute the spectral density using the Welch and Periodogram method
    freq_w, spec_w = signal.welch(re2e)
    freq_q, spec_q = signal.periodogram(re2e)
    # output re2e
    ofile = open('re2e.data', 'w')
    ofile.write('# End-to-end distance\n#\n')
    ofile.write('# Average: ' + str(rav) + '\n#\n')
    ofile.write('# frame\tre2e\n')
    for i in range(len(re2e)):
        ofile.write(str(i) + '\t' + str(re2e[i]) + '\n')
    ofile.close()
    # output histogram
    ofile = open('re2e_hist.data', 'w')
    ofile.write('# histogram of the end-to-end distance\n\n')
    ofile.write('# re2e\tp\n')
    for i in range(len(rhist[0])):
        ofile.write(str(rhist[1][i]*0.5 + rhist[1][i+1]*0.5) + '\t' + str(rhist[0][i]) + '\n')
    ofile.close()
    # output power spectras
    ofile = open('re2e_spectra_welch.data', 'w')
    ofile.write('# re2e power spectra from the Welche method\n#\n')
    ofile.write('# omega\tP\n')
    for i in range(len(freq_w)):
        ofile.write(str(freq_w[i]) + '\t' + str(spec_w[i]) + '\n')
    ofile.close()
    ofile = open('re2e_spectra_periodogram.data', 'w')
    ofile.write('# re2e power spectra from the periodogram method\n#\n')
    ofile.write('# omega\tP\n')
    for i in range(len(freq_q)):
        ofile.write(str(freq_q[i]) + '\t' + str(spec_q[i]) + '\n')
    ofile.close()
    
    return

############################################################################

if __name__ == '__main__':
    main()
