##############################################################################

## 
## Plot MSR of end-to-end vector for short filaments
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
    tau_A = (N+1)/Fmc
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
    
    data = pd.read_csv(path, sep='\t', skiprows=1, header=0)
    x = data['Timestep']
    y = data['MSR'] 
    
    return x, y
        
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
    savepath1 = savebase + 'msr.png'
    savepath1e = savebase + 'msr.eps'           
    textpath = '/MSR/msr.data'

    ## GAS PHASE
    # K = 2.5, xi/L = 0.1
    # fp = 0.08, Pe = 50 
    
    data_gas_path = datapath + 'kappa_2.5/fp_0.024'
    fullpath = data_gas_path + textpath
    time, theta_g = read_data(fullpath)

    init_gas_path = data_gas_path + '/init_info.txt'
    loadSimData(init_gas_path)
    
    time *= dt
    popt, pcov = curve_fit(theoretical_msr, time, theta_g)
    fit_g = theoretical_msr(time, popt[0])
    print 'Rotational diffusion for the gas phase ', popt[0]

    ## CLUSTER PHASE
    # K = 62.5, xi/L = 2.5
    # fp = 0.08, Pe = 50 
    
    data_cluster_path = datapath + 'kappa_62.5/fp_0.024'
    fullpath = data_cluster_path + textpath
    time, theta_c= read_data(fullpath)

    time *= dt
    popt, pcov = curve_fit(theoretical_msr, time, theta_c)
    fit_c = theoretical_msr(time, popt[0])
    print 'Rotational diffusion for the cluster phase ', popt[0]
    
    init_cluster_path = data_cluster_path + '/init_info.txt'
    loadSimData(init_cluster_path)

    ## GIANT CLUSTERS PHASE
    # K = 400.0, xi/L = 16.0
    # fp = 0.08, Pe = 50 
    
    data_giant_path = datapath + 'kappa_400.0/fp_0.024'
    fullpath = data_giant_path + textpath
    time, theta_gc = read_data(fullpath)

    time *= dt
    popt, pcov = curve_fit(theoretical_msr, time, theta_gc)
    fit_gc = theoretical_msr(time, popt[0])
    print 'Rotational diffusion for the giant cluster phase ', popt[0]

    init_giant_path = data_giant_path + '/init_info.txt'
    loadSimData(init_giant_path)  

    ## GAS HIGH PE
    # K = 2.5, xi/L = 0.1
    # fp = 7.0, Pe = 4300 
    
    data_gas_path = datapath + 'kappa_2.5/fp_7.0'
    fullpath = data_gas_path + textpath
    time, theta_g_high_pe = read_data(fullpath)

    init_gas_path = data_gas_path + '/init_info.txt'
    loadSimData(init_gas_path)
    
    time *= dt
    popt, pcov = curve_fit(theoretical_msr, time, theta_g_high_pe)
    fit_g_high_pe = theoretical_msr(time, popt[0])
    print 'Rotational diffusion for the gas phase ', popt[0]

    ## CLUSTER HIGH PE
    # K = 62.5, xi/L = 2.5
    # fp = 7.0, Pe = 4300 
    
    data_cluster_path = datapath + 'kappa_62.5/fp_7.0'
    fullpath = data_cluster_path + textpath
    time, theta_c_high_pe = read_data(fullpath)

    time *= dt
    popt, pcov = curve_fit(theoretical_msr, time, theta_c_high_pe)
    fit_c_high_pe = theoretical_msr(time, popt[0])
    print 'Rotational diffusion for the cluster phase ', popt[0]
    
    init_cluster_path = data_cluster_path + '/init_info.txt'
    loadSimData(init_cluster_path)

    ## GIANT CLUSTERS HIGH PE
    # K = 400.0, xi/L = 16.0
    # fp = 7.0, Pe = 4300 
    
    data_giant_path = datapath + 'kappa_400.0/fp_7.0'
    fullpath = data_giant_path + textpath
    time, theta_gc_high_pe = read_data(fullpath)

    time *= dt
    popt, pcov = curve_fit(theoretical_msr, time, theta_gc_high_pe)
    fit_gc_high_pe = theoretical_msr(time, popt[0])
    print 'Rotational diffusion for the giant cluster phase ', popt[0]

    init_giant_path = data_giant_path + '/init_info.txt'
    loadSimData(init_giant_path)      
    time = time/tau_D     
    
    ## plot the data
   
    ax0.loglog(time, theta_g, 'o', label="_nolegend_", color='blue', linewidth=2.0)
    ax0.loglog(time, theta_c, 'o', label="_nolegend_", color='green', linewidth=2.0)
    ax0.loglog(time, theta_gc, 'o', label="_nolegend_", color='red', linewidth=2.0)
    ax0.loglog(time, theta_g_high_pe, 'o', label="_nolegend_", color='cyan', linewidth=2.0)
    ax0.loglog(time, theta_c_high_pe, 'o', label="_nolegend_", color='yellow', linewidth=2.0)
    ax0.loglog(time, theta_gc_high_pe, 'o', label="_nolegend_", color='magenta', linewidth=2.0)
    
    ax0.loglog(time, fit_g, '--', label="_nolegend_", color='blue', linewidth=2.0)
    ax0.loglog(time, fit_c, '--', label="_nolegend_", color='green', linewidth=2.0)
    ax0.loglog(time, fit_gc, '--', label="_nolegend_", color='red', linewidth=2.0)   
    ax0.loglog(time, fit_g_high_pe, '--', label="_nolegend_", color='cyan', linewidth=2.0)
    ax0.loglog(time, fit_c_high_pe, '--', label="_nolegend_", color='yellow', linewidth=2.0)
    ax0.loglog(time, fit_gc_high_pe, '--', label="_nolegend_", color='magenta', linewidth=2.0) 
    
    ax0.loglog(time, time, '--', label='_nolegend_', linewidth=1.0, color='grey')
    #ax0.loglog(time*100, (time*100)**2, '--', label='_nolegend_', linewidth=1.0, color='grey')
    
    ## labels 
    
    ax0.set_xlabel("$t/\\tau_{D}$", fontsize=40)    
    ax0.set_ylabel("$\\Delta \\theta^{2}$", fontsize=40)   
    
    ## title
    
    #ax0.set_title("Spiral number of a tagged filament as a fnc. of time", fontsize=20)
    
    ## limits
    
    #ax0.set_xlim((1e-6, 1e-1))
    #ax0.set_ylim((1e-6,1e+1))
    
    ## ticks
    
    #ax0.xaxis.set_ticks([1e-6, 1e-4, 1e-2])
    #ax0.yaxis.set_ticks([1e-6, 1e-4, 1e-2, 1e+0])
    ax0.tick_params(axis='both', which='major', labelsize=40)

    ## legend

    ax0.legend(bbox_to_anchor=(0.02,0.,0.45,1.), loc=2, borderaxespad=0., \
        prop={'size': 20}, mode="expand", frameon=False)

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
    
    folderbase = "/local/duman/SIMULATIONS/" + args.length     # the data is in iff416 (assumption!)
    
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
    
    