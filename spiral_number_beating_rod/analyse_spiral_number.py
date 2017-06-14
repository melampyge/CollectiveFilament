#!/usr/local/bin/python2.7

import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import os


try:
    infilename = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '       spiral file'
    exit()

##################################################################

def read_data(infilename):
    """ read in the spiral number"""
    ifile = open(infilename)
    ifile.readline()
    ifile.readline()
    t = []
    s = []
    for line in ifile:
        line = line.split()
        t.append(int(line[0]))
        s.append(float(line[1]))
    ifile.close()
    t = np.array(t)
    s = np.array(s)
    return t, s
    
##################################################################

def main():
    """ main function, called when the script is started"""
    ### read in the spiral number
    t, spiral = read_data(infilename)

    fig = plt.figure()
    ax = plt.subplot(111)
    ax.hist(spiral, bins = 51)
    plt.show()
    plt.close()

 


            
##################################################################

if __name__ == '__main__':
    main()
    
