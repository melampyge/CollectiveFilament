#!/usr/local/bin/python2.7

import numpy as np
import matplotlib.pyplot as plt
import sys, os
import math


class VoronoiDensity:

    #############################################################

    def gen_mol_info(self,natoms,npol):
        """ generate information on molecules"""
        mol = np.zeros((natoms), dtype = int)
        nmol = natoms/npol
        k = 0
        for i in range(nmol):
            for j in range(npol):
                mol[k] = i
                k = k + 1
        return mol,nmol


    #############################################################

    def __init__(self, nsteps, natoms, npol, line):
        """ initialize, allocate required arrays and precompute
            misc required information"""
        ### set the options
        line = line.split()
        self.nbins = int(line[1])
        self.vmin = float(line[2])
        self.vmax = float(line[3])
        ### initialize some information on molecules
        self.mol, self.nmol = self.gen_mol_info(natoms,npol)
        ### allocate results array
        self.hist_vol = np.zeros((nsteps,self.nbins))
        return

    #############################################################

    def compute(self,step,x,y,lx,ly,natoms,plot = 'False'):
        """ compute Voronoi tessellation"""
        ### generate input script for voro++
        ofile = open('voro.data', 'w')
        for i in range(natoms):
            ofile.write(str(i) + ' ' + str(x[i]) + ' ' + str(y[i]) + ' 0.5\n')
        ofile.close()
        ### perform Voronoi tessellation using voro++
        os.system('/usr/users/iff_th2/isele/Applications/voro++-0.4.6/src/voro++ -p 0.0 ' + str(lx) + ' 0.0 ' + str(ly) + ' 0.0 1.0 voro.data')
        ### read in the results
        vol = np.zeros((self.nmol))
        ifile = open('voro.data.vol')
        for i in range(natoms):
            line = ifile.readline()
            line = line.split()
            idx = int(line[0])
            v = float(line[4])
            vol[self.mol[idx]] += v
        ifile.close()
        ### remove voro++ files
        os.system('rm voro.data voro.data.vol')
        ### store results in histogram
        self.hist_vol[step], self.edges = np.histogram(np.log(vol), bins = self.nbins, range = (self.vmin, self.vmax))
        if plot == 'True':
            plt.plot(self.edges[1:],self.hist_vol[step])
            plt.show()
            plt.close()
        return

