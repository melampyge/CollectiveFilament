##############################################################################

## 
## Plot the normalized distribution of absolute value of spiral number
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
    
    data = pd.read_csv(path, sep='\t', skiprows=1, header=0)
    x = (data['upper_edge'] + data['lower_edge'])/2
    y = data['number'] 
    y = np.asarray(y)

    ## normalize distribution

    ysum = np.sum(y)
    for el in range(len(y)):
        y[el] /= (0.05*ysum)   
        
    ## check normalization
        
    print "Normalization of ", path, " is : ", np.sum(y[:-1]*np.diff(x))
    
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
    
    downlim = 1e-6
    uplim = 1e+4
    num_ticks = 9
    
    ## data access structure
    
    datapath = folderbase + 'density_0.2/'
    savepath1 = savebase + 'spiral_dist.eps'
    savepath2 = savebase + 'spiral_dist.png'    
    textpath = '/SPIRALS/spiral_histogram.data'
            
    ## GAS PHASE
    
    data_gas_path = datapath + 'kappa_5.0/fp_0.015'
    fullpath = data_gas_path + textpath
    s_gas, ps_gas = read_data(fullpath)
    
    init_gas_path = data_gas_path + '/init_info.txt'
    loadSimData(init_gas_path)
    pe_gas = Pe
    xil_gas = persistence

    ## LOOP PHASE
    
    data_loop_path = datapath + 'kappa_1600.0/fp_0.075'
    fullpath = data_loop_path + textpath
    s_loop, ps_loop = read_data(fullpath)

    init_loop_path = data_loop_path + '/init_info.txt'
    loadSimData(init_loop_path)
    pe_loop = Pe
    xil_loop = persistence

    ## SPIRALS PHASE
    
    data_spiral_path = datapath + 'kappa_5.0/fp_1.0'
    fullpath = data_spiral_path + textpath
    s_spiral, ps_spiral = read_data(fullpath)

    init_spiral_path = data_spiral_path + '/init_info.txt'
    loadSimData(init_spiral_path)
    pe_spiral = Pe
    xil_spiral = persistence    
        
    ## plot the data
    
    #ax0.semilogy(s_gas, ps_gas, label="Gas phase", linewidth=2.0, color='blue')
    ax0.semilogy(s_gas, ps_gas, 'o', label="_nolegend_", color='blue', markersize=7)
    
    #ax0.semilogy(s_loop, ps_loop, label="Loop phase", linewidth=2.0, color='green')
    ax0.semilogy(s_loop, ps_loop, 'o', label="_nolegend_", color='green', markersize=7)

    #ax0.semilogy(s_spiral, ps_spiral, label="Spirals phase", linewidth=2.0, color='red')    
    ax0.semilogy(s_spiral, ps_spiral, 'o', label="_nolegend_", color='red', markersize=7)
    
    ## labels 
    
    ax0.set_xlabel("$|s|$", fontsize=30)    
    ax0.set_ylabel("$P(|s|)$", fontsize=30)   
    
    ## title
    
    #ax0.set_title("Probability distribution of spiral number", fontsize=20)
    
    ## limits

    ax0.set_xlim(0,8)
    ax0.set_ylim(downlim,uplim)    
    
    ## ticks
    
    #ax0.xaxis.set_ticks([0, 1, 2, 3, 4])
    ax0.yaxis.set_ticks([1e-6, 1e-4, 1e-2, 1e+0, 1e+2, 1e+4])
    ax0.tick_params(axis='both', which='major', labelsize=20)

    ## legend

    ax0.legend(bbox_to_anchor=(0.60,0.,0.45,1.), loc=2, borderaxespad=0., \
        prop={'size': 20}, mode="expand", frameon=False)

    ## save figure

    plt.savefig(savepath1, dpi=200, bbox_inches='tight', pad_inches=0.08)
    plt.savefig(savepath2, dpi=200, bbox_inches='tight', pad_inches=0.08)
    
    print ps_gas
    print ps_loop
    print ps_spiral
    #print ps_gas-ps_loop

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
    
    