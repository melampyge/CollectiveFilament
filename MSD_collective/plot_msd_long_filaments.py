##############################################################################

## 
## Plot MSD of long filaments
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

##############################################################################

def loadSimData(datafile):
    """ load initial simulation data"""
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
    dtSamp, T, box_area, nt, body_length, Pe, persistence, flexure, tau_D, tau_A 

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
    
    dt = 0.0001
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
    print "Diffusive time scale is ", tau_D
    print "Advective time scale is ", tau_A
    
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
    
    #data = pd.read_csv(path, sep='\t', skiprows=1, header=0)
    #x = data['Timestep']
    #y = data['MSD'] 
    
    data = np.transpose(np.loadtxt(path, dtype=np.float64))
    x = data[0]
    y = data[1]
    
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
    savepath1 = savebase + 'msd_long_filaments.png'
    savepath1e = savebase + 'msd_long_filaments.eps'    
    savepath2 = savebase + 'msd_exponents_long_filaments.png'    
    savepath2e = savebase + 'msd_exponents_long_filaments.eps'        
    textpath = '/MSD/msd_alternative.data'

    ## GAS PHASE
    
    data_gas_path = datapath + 'kappa_5.0/fp_0.075'
    fullpath = data_gas_path + textpath
    time, msd_gas = read_data(fullpath)
    
    init_gas_path = data_gas_path + '/init_info.txt'
    loadSimData(init_gas_path)
    pe_gas = Pe
    xil_gas = persistence
    
    msd_gas /= body_length**2

    ## LOOP PHASE
    
    data_loop_path = datapath + 'kappa_1600.0/fp_0.075'
    fullpath = data_loop_path + textpath
    time, msd_loop = read_data(fullpath)

    init_loop_path = data_loop_path + '/init_info.txt'
    loadSimData(init_loop_path)
    pe_loop = Pe
    xil_loop = persistence
    
    msd_loop /= body_length**2
    
    ## SPIRALS PHASE
    
    data_spiral_path = datapath + 'kappa_5.0/fp_1.0'
    fullpath = data_spiral_path + textpath
    time, msd_spiral = read_data(fullpath)

    init_spiral_path = data_spiral_path + '/init_info.txt'
    loadSimData(init_spiral_path)
    pe_spiral = Pe
    xil_spiral = persistence  
    
    msd_spiral /= body_length**2
    
    time = time/tau_D
    
    logx = np.log(time[::3])
    logy = np.log(msd_gas[::3])
    
    alpha_gas = np.diff(logy)/np.diff(logx)
    
    logy = np.log(msd_loop[::3])
    alpha_loop = np.diff(logy)/np.diff(logx)
    
    logy = np.log(msd_spiral[::3])
    alpha_spiral = np.diff(logy)/np.diff(logx)    
    
    print '** GAS PHASE **\n\n'
    print alpha_gas, '\n'

    print '** LOOP PHASE **\n\n'
    print alpha_loop, '\n'

    print '** SPIRAL PHASE **\n\n'
    print alpha_spiral, '\n'    
    
    ## plot the data
    
    ax0.loglog(time, msd_gas, label="_nolegend_", linewidth=2.0)
    ax0.loglog(time, msd_loop, label="_nolegend_", linewidth=2.0)
    ax0.loglog(time, msd_spiral, label="_nolegend_", linewidth=2.0)
    ax0.loglog(time, time, '--', label='_nolegend_', linewidth=1.0, color='grey')
    ax0.loglog(time*100, (time*100)**2, '--', label='_nolegend_', linewidth=1.0, color='grey')
    
    ## labels 
    
    ax0.set_xlabel("$t/\\tau_{D}$", fontsize=30)    
    ax0.set_ylabel("$\\Delta r^{2}/L^{2}$", fontsize=30)   
    
    ## title
    
    #ax0.set_title("Spiral number of a tagged filament as a fnc. of time", fontsize=20)
    
    ## limits
    
    ax0.set_xlim((1e-6, 1e-1))
    ax0.set_ylim((1e-6,1e+1))
    
    ## ticks
    
    ax0.xaxis.set_ticks([1e-6, 1e-4, 1e-2])
    ax0.yaxis.set_ticks([1e-6, 1e-4, 1e-2, 1e+0])
    ax0.tick_params(axis='both', which='major', labelsize=20)

    ## legend

    ax0.legend(bbox_to_anchor=(0.02,0.,0.45,1.), loc=2, borderaxespad=0., \
        prop={'size': 20}, mode="expand", frameon=False)

    ## save figure

    plt.savefig(savepath1, dpi=200, bbox_inches='tight', pad_inches=0.08)
    plt.savefig(savepath1e, dpi=200, bbox_inches='tight', pad_inches=0.08)    
    plt.clf()
    
    
    ## plot the data
    
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
    
    print len(time[::3][1:-1]), len(alpha_gas[1:])
    
    ax0.plot(time[::3][1:-1], alpha_gas[1:], label="_nolegend_", linewidth=2.0)
    ax0.plot(time[::3][1:-1], alpha_loop[1:], label="_nolegend_", linewidth=2.0)
    ax0.plot(time[::3][1:-1], alpha_spiral[1:], label="_nolegend_", linewidth=2.0)
    #ax0.set_xscale('log')
    #ax0.set_yscale('log')
    
    ## labels 
    
    ax0.set_xlabel("$t/\\tau_{D}$", fontsize=60)    
    ax0.set_ylabel("$\\alpha$", fontsize=60)   
    
    ## title
    
    #ax0.set_title("Spiral number of a tagged filament as a fnc. of time", fontsize=20)
    
    ## ticks
    
    ax0.xaxis.set_ticks([0.000, 0.004,  0.008])
    ax0.yaxis.set_ticks([0.6, 1.0, 1.5, 2.0])
    ax0.tick_params(axis='both', which='major', labelsize=40)

    ## legend

    ax0.legend(bbox_to_anchor=(0.02,0.,0.45,1.), loc=2, borderaxespad=0., \
        prop={'size': 20}, mode="expand", frameon=False)    

    ## save figure

    plt.savefig(savepath2, dpi=200, bbox_inches='tight', pad_inches=0.08)
    plt.savefig(savepath2e, dpi=200, bbox_inches='tight', pad_inches=0.08)
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

    savebase = "/usr/users/iff_th2/duman/RolfData/SPIRAL_FORMATION/"
    os.system("mkdir -p " + savebase)
    
    ## plot certain time frames to generate a movie later on
               
    plot_data(folderbase, savebase)    
    
    return
    
##############################################################################

if __name__ == "__main__":
    main()  
    
##############################################################################   
    
    