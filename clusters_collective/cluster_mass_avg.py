#!/usr/bin/python


# Load needed libraries and necessary files
import argparse
import numpy as np
#import matplotlib as mpl
#mpl.use('Agg')
#import matplotlib.pyplot as plt
import sys
import os

#
# Function definitions
#
 
# Load initial simulation data      
def loadSimData(datafile):
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
    dtSamp, T, box_area, nt, body_length, Pe, persistence, flexure 

    datafile = open(datafile,"r")
    for line in datafile:
        A = line.split()
        if A[0] == "dt":                    # Time interval between MD steps
            dt = float(A[-1])
        elif A[0] == "ti":                  # Beginning time for data acquisition
            ti = float(A[-1])
        elif A[0] == "Lx":                  # Box size in x
            Lx = float(A[-1])            
        elif A[0] == "Ly":                  # Box size in y
            Ly = float(A[-1])
        elif A[0] == "totalStep":           # Total MD steps
            totalStep = float(A[-1])
        elif A[0] == "nsamp":               # Data sampling frequency
            nsamp = float(A[-1])
        elif A[0] == "nfil":                # Number of particles per polymer
            N = float(A[-1])
        elif A[0] == "L":                   # Number of particles
            L = float(A[-1])
        elif A[0] == "B":                   # Bond length between particles of a body
            B = float(A[-1])
        elif A[0] == "kT":                  # Boltzmann constant*Temperature
            kT = float(A[-1])
        elif A[0] == "Fmc":                 # Self propulsion force constant
            Fmc = float(A[-1])     
        elif A[0] == "Kbend":               # Bending constant
            Kbend = float(A[-1])
    
    Lx /= B
    Ly /= B
    M = L/N
    dtSamp = dt*nsamp
    T = totalStep - ti
    nt = T/nsamp
    box_area = Lx*Ly
    body_length = B*N
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    flexure = Pe/persistence
        

#
# Class definitions
#

# Particle data
class Particles:
    
    def __init__(self, path):
        file = np.transpose(np.loadtxt(path, dtype=float))
        self.xi = file[0]/B                 # Image particle positions in x
        self.yi = file[1]/B                 # Image particle positions in y 
        self.phi = file[2]                  # Bead orientation 
        self.cidx = file[3]                 # Cluster index
        
# Cluster sizes
class ClusterSize:
    
    def __init__(self, path):
        file = np.loadtxt(path, dtype=float)
        self.cs = file                     # Cluster sizes
 

# Argument parsing (command line options)
parser = argparse.ArgumentParser()
parser.add_argument("initfile", help="File containing initial simulation data")
args = parser.parse_args()

# Load saved preliminary simulation data into relevant variables
loadSimData(args.initfile+"/init_info.txt")
datafile = args.initfile+"/CLUSTER"

cs_cnt = np.zeros((M+1,1))
cnt = 0

# Time frame loop
for frame in np.arange(int(ti),int(totalStep),int(nsamp)):
    
    
    print frame , ' of ', totalStep
    cnt += 1.
    
    #
    # Load and set the data
    #  
    
    # Load particle data  -- 1xL --
    #path = datafile + '/beads_' + str(frame) + '.txt'
    #p = Particles(path)                    

    
    # Load cluster sizes
    path = datafile + '/cluster_sizes_' + str(frame) + '.txt'
    if os.path.exists(path):
        clsizes = ClusterSize(path)
        for siz in clsizes.cs:
            cs_cnt[siz] += 1
    
savefile = datafile + "/avg_cl_size_cnt.txt"
sfile = open(savefile, 'w')
for siz in np.arange(1,M+1):
    cs_cnt[siz] /= cnt
    sfile.write(str(siz) + "\t" + str(cs_cnt[siz]) + "\n")
    
sfile.close()