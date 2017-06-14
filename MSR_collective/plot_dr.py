##############################################################################

## 
## Plot rotational diffusion coeff. of end-to-end vector for short filaments from MSR
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
        dtSamp, T, box_area, nt, body_length, pe, xil, flexure, tau_D, tau_A, \
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
    pe = Fmc*body_length**2/kT
    xil = Kbend/(kT*body_length)
    flexure = pe/xil
    T = totalStep - ti
    nt = T/nsamp
    nsamp = int(nsamp)
    tau_D = body_length**2*(N+1)/4/kT
    if Fmc != 0:
        tau_A = (N+1)/Fmc
    else:
        tau_A = 0.000001
    nmol = int(M)
    nfil = int(N)
    nbeads = int(L)
    
    
    print '\n\n*** SIMULATION PARAMETERS FOR ', datafile
    print 'dt = ', dt
    print 'L = ', body_length
    print 'Pe = ', pe
    print 'xi_p/L = ', xil
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

class Subplots:
    """ plot structure"""
    
    totcnt = -1             # Total number of subplots 
    
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction
        
    def addSubplot(self):
        """ add a subplot in the grid structure"""
        
        ## increase the number of subplots in the figure
        
        self.totcnt += 1
        
        ## get indices of the subplot in the figure
        
        self.nx = self.totcnt%(self.tot)
        self.ny = self.totcnt/(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])

##############################################################################

def read_data(path):
    """ read data"""
  
    ## read data
    
    data = np.loadtxt(path, dtype=np.float64)
    
    return data
        
##############################################################################

def plot_data(folderbase, savebase):
    """ plot data"""

    ## set plot properties

    ax_len = 0.9                          # Length of one subplot square box
    ax_b = 0.1                            # Beginning/offset of the subplot in the box
    ax_sep = 0.1                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots    
    fig = plt.figure()
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
    
    ## set more plot properties
    
    downlim = -5
    #uplim = max(Lx,Ly)+5
    num_ticks = 5
    
    ## data access structure
    
    datapath = folderbase + 'density_0.2/'
    savepath1 = savebase + 'dr.png'
    savepath1e = savebase + 'dr.eps'           
    textpath = 'dr.data'
    
    ## data points in the parameter space
    
    kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 400.0]
    fp = [0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
    ## load all the dr values as a fnc. of pe
    
    pe_number = []
    dr = []
        
    for k in kappa:
        
        dr_per_kappa = []
        
        for f in fp:
                        
            path = datapath + 'kappa_' + str(k) + '/fp_' + str(f) + '/'
            
            initpath = path + 'init_info.txt'
            assert os.path.exists(initpath), 'init_info.txt does not exist for ' + path
            loadSimData(initpath)
            
            ## save the pe
            
            if k == kappa[-1]:
                pe_number.append(pe)
                
            ## save the dr
                
            datafile = path + textpath
            if os.path.exists(datafile):
                assert os.path.exists(datafile), 'data does not exist for ' + datafile
                dr_value = read_data(datafile)
            else:
                dr_value = 0.
                print path
                exit()
            
            print 'for kappa = ' + str(k) + ' and fp = ' + str(f) + ' Dr = ' + str(dr_value) + '\n'
            
            dr_per_kappa.append(dr_value)
            
        dr.append(dr_per_kappa)
        
        
    ## plot data with pe in the x axis and xil in the legend

    color=iter(plt.cm.jet(np.linspace(0,1,len(kappa))))
        
    for j, k in enumerate(kappa):
           
        ## plot the data
        
        c = next(color)
        ax0.loglog(pe_number, dr[j], 'o', label='_nolegend_', markersize=7, color=c)
        #ax0.loglog(pe_number, dr[j], 'o', label=str(k), markersize=7, color=c)        
        ax0.loglog(pe_number, dr[j], '--', label='_nolegend_', linewidth=2., color=c)
        
        ## limits
        
        #ax0.set_xlim((1e-6, 1e-1))
        #ax0.set_ylim((1e-6,1e+1))
        
        ## ticks
        
        #ax0.xaxis.set_ticks([1e-6, 1e-4, 1e-2])
        #ax0.yaxis.set_ticks([1e-6, 1e-4, 1e-2, 1e+0])
        ax0.tick_params(axis='both', which='major', labelsize=40)

    ## labels 
    
    ax0.set_xlabel("$Pe$", fontsize=40)    
    ax0.set_ylabel("$D_{r}\\tau_{D}$", fontsize=40)   
    
    ## title
    
    #ax0.set_title("Spiral number of a tagged filament as a fnc. of time", fontsize=20)
        
    ## legend

#    ax0.legend(bbox_to_anchor=(0.02,0.,0.45,1.), loc=2, borderaxespad=0., \
#        prop={'size': 20}, mode="expand", frameon=False)

    ## save figure

    plt.savefig(savepath1, dpi=200, bbox_inches='tight', pad_inches=0.08)
    plt.savefig(savepath1e, dpi=200, bbox_inches='tight', pad_inches=0.08)    
    plt.clf()
 
    return

##############################################################################  

def main():

    ## do argument parsing
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-l", "--length", help="Filament length, like many_polymers_5/ or long_filaments/")
    args = parser.parse_args()
    
    ## data access structure
    
    folderbase = "/local/duman/SIMULATIONS/many_polymers_5/"     # the data is in iff416 (assumption!)
    
    ## set save folder paths 

    savebase = "/usr/users/iff_th2/duman/RolfData/OTHER_PLOTS/"
    os.system("mkdir -p " + savebase)
    
    ## plot certain time frames to generate a movie later on
               
    plot_data(folderbase, savebase)    
    
    return
    
##############################################################################

if __name__ == "__main__":
    main()  
    
##############################################################################   
    
    