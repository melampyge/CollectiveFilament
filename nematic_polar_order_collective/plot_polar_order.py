##########################################################################

## 
## Plot the normalized distribution of cluster masses
##

##########################################################################

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

##########################################################################

def power_law(x, a, b):
    return a * x**(-b)

##########################################################################

def exp_law(x, a, b):
    return a * np.exp(-x/b)

##########################################################################

def power_law_with_exp_tail(x, a, b):
    return (x**(-a))*np.exp(-x/b)
    
##########################################################################

def power_law_with_up_exp_tail(x, a, b):
    return (x**(-a))*np.exp(x/b)    

##########################################################################

def mixed_power_law_with_exp_tail(x, a, b, c, d, f):
    return (x**(-a))*np.exp(-x/b) + c*(x**(d))*np.exp(-x/f)

##########################################################################
    
def mixed_power_law_with_exp_tail_2(x, a, b, c, d):
    return (x**(b))*np.exp(-x/c) + a*(x**(-d))  

##########################################################################

def stretch_exp_law(x, a, b, c):
    return a * np.exp(-(x/c)**b)  
    
##########################################################################

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

##########################################################################

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

##########################################################################

def read_data(path):
    """ read data"""
    
    ## read data
    
    data = pd.read_csv(path, sep='\t', skiprows=4, header=None) 
    cols = np.arange(1,len(data.columns))
    t = data[0]
    y = data[cols].mean(axis=0)
    y.index = cols-1
    ystd = y/np.sqrt(len(t))
    ystd.index = cols-1
    
    order_data = pd.read_csv(path, sep='\t', skiprows=3, nrows=1, header=None)
    x = np.asarray(order_data[cols].transpose())
    
    xdata = np.ones((16,1))
    ydata = np.ones((16,1))
    xdata2 = np.transpose(np.array(x, dtype=float))
    ydata2 = np.array(y, dtype=float)
    for iii in np.arange(16):
        xdata[iii] = xdata2[0][iii]
        ydata[iii] = ydata2[iii]
    xdata = xdata[:,0]
    ydata = ydata[:,0]  

    
    return xdata, ydata, x, y, ystd
       
##########################################################################

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
    savepath1 = savebase + 'polar_order.eps'
    savepath2 = savebase + 'polar_order.png'    
    textpath = '/ORDER_PARAMETER/polar_order_parameter.data'
            
    ## GAS PHASE
    # K = 2.5, xi/L = 0.1
    # fp = 0.08, Pe = 15 
    
    data_gas_path = datapath + 'kappa_2.5/fp_0.08'
    fullpath = data_gas_path + textpath
    xfit, yfit, w, Sg, std_g = read_data(fullpath)
    
    popt, pcov = curve_fit(stretch_exp_law, xfit, yfit)
    fit_g = stretch_exp_law(xfit, popt[0], popt[1], popt[2])
    print 'Stretched exponential fit for the gas phase with ', popt[0], popt[1], popt[2]
    
    init_gas_path = data_gas_path + '/init_info.txt'
    loadSimData(init_gas_path)

    ## CLUSTER PHASE
    # K = 62.5, xi/L = 2.5
    # fp = 0.08, Pe = 15 
    
    data_cluster_path = datapath + 'kappa_62.5/fp_0.08'
    fullpath = data_cluster_path + textpath
    xfit, yfit, w, Sc, std_c = read_data(fullpath)

    print np.shape(xfit)
    popt, pcov = curve_fit(stretch_exp_law, xfit[2:15], yfit[2:15])
    fit_c = stretch_exp_law(xfit[2:15], popt[0], popt[1], popt[2])
    print 'Stretched exponential fit for the cluster phase with ', popt[0], popt[1], popt[2]
    
    init_cluster_path = data_cluster_path + '/init_info.txt'
    loadSimData(init_cluster_path)

    ## GIANT CLUSTERS PHASE
    # K = 400.0, xi/L = 16.0
    # fp = 0.08, Pe = 15 
    
    data_giant_path = datapath + 'kappa_400.0/fp_0.08'
    fullpath = data_giant_path + textpath
    xfit, yfit, w, Sgc, std_gc = read_data(fullpath)

    popt, pcov = curve_fit(stretch_exp_law, xfit[5:10], yfit[5:10])
    fit_gc = stretch_exp_law(xfit[5:10], popt[0], popt[1], popt[2])
    print 'Stretched exponential fit for the giant cluster phase with ', popt[0], popt[1], popt[2]

    init_giant_path = data_giant_path + '/init_info.txt'
    loadSimData(init_giant_path)  
        
    ## plot the data
        
    w /= B
    
    ax0.errorbar(w, Sg, yerr=std_g, label="_nolegend_", linewidth=2.0, color='blue')
    #ax0.plot(w, fit_g, '--', label="_nolegend_", linewidth=2.0, color='blue')
    
    ax0.errorbar(w, Sc, yerr=std_c, label="_nolegend_", linewidth=2.0, color='green')
    #ax0.plot(w[2:15], fit_c, '--', label="_nolegend_", color='green', linewidth=2.0)

    ax0.errorbar(w, Sgc, yerr=std_gc, label="_nolegend_", linewidth=2.0, color='red')    
    #ax0.plot(w[5:10], fit_gc, '--', label="_nolegend_", color='red', linewidth=2.0)
    
    ax0.set_xscale('log')
    
    ## labels 
    
    ax0.set_xlabel("$w/r_{0}$", fontsize=40)    
    ax0.set_ylabel("$S_{1}$", fontsize=40)   
    
    ## title
    
    #ax0.set_title("Probability distribution of spiral number", fontsize=20)
    
    ## limits

    ax0.set_xlim(1e+0,3*1e+3)
    ax0.set_ylim(0,1.1)    
    
    ## ticks
    
    ax0.xaxis.set_ticks([1e+0, 1e+1, 1e+2, 1e+3])
    ax0.yaxis.set_ticks([0.0, 0.5, 1.0])
    ax0.tick_params(axis='both', which='major', labelsize=40)

    ## legend

    ax0.legend(bbox_to_anchor=(0.60,0.,0.45,1.), loc=2, borderaxespad=0., \
        prop={'size': 20}, mode="expand", frameon=False)

    ## save figure

    plt.savefig(savepath1, dpi=200, bbox_inches='tight', pad_inches=0.08)
    plt.savefig(savepath2, dpi=200, bbox_inches='tight', pad_inches=0.08)
    
    print w
    print Sg
    print Sc
    print Sgc

    return

##########################################################################  

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
    
##########################################################################

if __name__ == "__main__":
    main()  
##########################################################################   
    
    