#!/usr/local/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys
import math

##################################################################
#    global variables

# number of fractions into which to decompose the box
frac = [100., 70., 50., 30., 20., 10., 7., 5., 3., 2.]
frac = np.array(frac)
nfrac = len(frac)

class NumberFluctuation:


    ##############################################################

    def __init__(self, nsteps, natoms,line):
        """ initialize: allocate arrays to store results"""
        ### array to store number fluctuations
        self.dn = np.zeros((nsteps, nfrac))
        ### compute ``x-axis'' values
        self.nav = float(natoms)/frac**2
        return

    ##############################################################

    def compute(self, step,xs,ys,plot='False'):
        """ compute number fluctuations"""
        for j in range(nfrac):
            ### compute histogram of the number of bins in each cell
            hist, xedges, yedges = np.histogram2d(xs,ys,bins=frac[j], range = [[0,1],[0,1]])
            ### compute the standard deviation
            std = np.std(hist)
            self.dn[step,j] = std
        if plot == 'True':
            plt.loglog(self.nav, self.dn[step],marker = 'o')
            plt.show()
            plt.close()
        return
