#!/usr/bin/python


# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
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
    box_area = Lx*Ly
    body_length = B*N
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    flexure = Pe/persistence
    T = totalStep - ti
    nt = T/nsamp
                    

# Argument parsing (command line options)
parser = argparse.ArgumentParser()
parser.add_argument("initfile", help="File containing initial simulation data")
parser.add_argument("savefile", help="File in which figures should be saved")
args = parser.parse_args()

# Load saved preliminary simulation data into relevant variables
loadSimData(args.initfile+"/init_info.txt")
clusterfile = args.initfile+"/CLUSTER"
histofile = args.initfile+"/HISTOGRAMS/tables"

# Plot properties
downlim = -5
uplim = max(Lx,Ly)+5
ax_len = 0.38                         # Length of one subplot square box
ax_b = 0.1                            # Beginning/offset of the subplot in the box
ax_sep = 0.15                          # Separation length between two subplots
total_subplots_in_y = 2               # Total number of subplots
tick_interval = int(uplim/5)
downlim_zoom = 800
uplim_zoom = 1600
tick_interval_zoom = int((uplim_zoom - downlim_zoom)/5)

# Read the preliminary data of each file   
rho_path = histofile + '/rho.data'
rho_file = open(rho_path, 'r')

wor_path = histofile + '/worticity.data'
wor_file = open(wor_path, 'r')
for i in np.arange(7):
    wor_file.readline()

# Read total number of steps
line = rho_file.readline()   
line = line.split()
nsteps = int(line[-1])

# Read number of bins in x and y directions
line = rho_file.readline()   
line = line.split()
nx = int(line[-1])

line = rho_file.readline()   
line = line.split()
ny = int(line[-1])

xedges = np.zeros((nx, ny))
yedges = np.zeros((nx, ny))

# Read the first frame
rho_file.readline()

# Read the bin edges
line = rho_file.readline()   
line = line.split()
for i in range(nx):
    xedges[i] = float(line[i])

rho_file.readline()

line = rho_file.readline()   
line = line.split()
for i in range(ny):
    yedges[i] = float(line[i])


xlin = np.linspace(0., Lx, nx)
ylin = np.linspace(0., Ly, ny)
xgrid, ygrid = np.meshgrid(xlin, ylin)


# Time averaged histograms
rho = []
#wor = []

# Time frame loop
print "nsteps = ", nsteps
for frame in np.arange(nsteps):
    
    
    print frame 
    
    # Read the lines corresponding to the current frame 
#    rho = np.zeros((nx, ny))
#    wor = np.zeros((nx, ny))
    
    line = rho_file.readline()
    line = line.split()
    timestep = float(line[-1])

#    wor_file.readline()
    
    
    for j in np.arange(nx):
        
        rho_line = rho_file.readline()
        rho_line = rho_line.split()
        
#        wor_line = wor_file.readline()
#        wor_line = wor_line.split()
        
        for k in np.arange(ny):
            rho.append(float(rho_line[k]))
#            rho[j, k] = float(rho_line[k])
#            wor[j, k] = float(wor_line[k])
            

savepath = args.savefile + '/tables/rho_hist.data'
savefile = open(savepath, 'w')
for i in rho:
    savefile.write(str(i) + '\n') 
        

