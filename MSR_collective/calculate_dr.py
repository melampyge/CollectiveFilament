##############################################################################

## 
## Calculate rotational diffusion coefficient from MSR of the end-to-end vector
##

##############################################################################

## load needed libraries and necessary files

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import h5py
import os
import pandas as pd
from scipy.optimize import curve_fit

##############################################################################

def loadSimData(datafile):
    """ load initial simulation data"""
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
        dtSamp, T, box_area, nt, body_length, Pe, persistence, flexure, tau_D, tau_A, \
            nbeads, nmol, nfil

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
    body_length = B*(N-1)
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    flexure = Pe/persistence
    T = totalStep - ti
    nt = T/nsamp
    nsamp = int(nsamp)
    tau_D = body_length**2*(N+1)/4/kT
    if Fmc != 0.:
        tau_A = (N+1)/Fmc
    else:
        tau_A = 0.0000001
    nmol = int(M)
    nfil = int(N)
    nbeads = int(L)
    
    print '\n\n*** SIMULATION PARAMETERS FOR ', datafile
    print 'dt = ', dt
    print 'L = ', body_length
    print 'Pe = ', Pe
    print 'xi_p/L = ', persistence
    print 'T = ', T
    print 'nfil = ', nfil
    print 'nmol = ', nmol
    print 'nbeads = ', nbeads
    print 'lx = ly = ', Lx
    print "Diffusive time scale is ", tau_D
    print "Advective time scale is ", tau_A
    print '***\n\n'
    
    return

##############################################################################

def theoretical_msr(t,a):
    """ Theoretical MSR"""
    
    return 2*a*t
    
##############################################################################    

def read_data(datafile):
    """ Read the data with black magic"""
    
    data = pd.read_csv(datafile, sep='\t', skiprows=1, header=0) 
    x = data['Timestep']
    y = data['MSR']

    return x, y
  
##############################################################################
    
def format_data(datafile, density, kappa, fp):
    """ Format the data"""

    ## Continue indexing the data 

    folders = []
    for d in density:
        for k in kappa:
            for f in fp:
                folders.append( datafile + 'density_' + str(d) + \
                    '/kappa_' + str(k) + '/fp_' + str(f) + '/' )
                    
    return folders    
       
##############################################################################
   
def main(): 
    
    ## Index the data
    
    density = [0.2]
    kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 400.0]
    fp = [0.0, 0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
    ## Format the data
    
    database = "/local/duman/SIMULATIONS/many_polymers_5/"
    folders = format_data(database, density, kappa, fp)    
       
    ## Save the data
       
    ofname = 'dr.data'
       
    for f in folders:
        
        ## Load initial simulation data
        
        initfile = f + "init_info.txt"
        loadSimData(initfile)
        
        ## Load analysis data
        
        datafile = f + 'MSR/msr.data'
        
        if os.path.exists(datafile):
            
            ## Use black magic to read the text file
            
            time, msr = read_data(datafile)
            
            ## Perform fitting and save the fit
            
            time *= dt
            #x = np.log(time)
            #y = np.log(msr)
            
            popt, pcov = curve_fit(theoretical_msr, time, msr)
            dr = popt[0]
            print 'Rotational diffusion for file ', f, dr
            
            savefilename = f + ofname
            savefile = open(savefilename, 'w')
            savefile.write(str(dr) + '\n')
            savefile.close()
   
        else:
            print "Analysis file could not be loaded for : ", f

##############################################################################            

if __name__ == '__main__':
    main()

##############################################################################
