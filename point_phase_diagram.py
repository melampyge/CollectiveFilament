
## Load needed libraries and necessary files
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import pandas as pd
import argparse
import pylab

################################################################################

#
# Function definitions
#
 
##################################################################
 
def loadSimData(datafile):
    """ Load initial simulation data"""
    
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
    
    L = int(L)
    N = int(N)
    Lx /= B
    Ly /= B
    M = L/N
    dtSamp = dt*nsamp
    T = totalStep - ti
    nt = T/nsamp
    box_area = Lx*Ly
    body_length = B*(N-1)
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    flexure = Pe/persistence
    
    return

##################################################################

class Subplots:
    """ Arrange subplot grid structure (square box is assumed)"""
    
    totcnt = -1             # Total number of subplots 
    
    ## constructor
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction
        
    ## add a subplot in the grid structure
    def addSubplot(self):
        
        # Increase the number of subplots in the figure
        self.totcnt += 1
        
        # Indices of the subplot in the figure
        self.nx = self.totcnt%(self.tot)
        self.ny = self.totcnt/(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])

##################################################################

class Phase:
    """ Properties of the point in phase space"""        
        
        
    def __init__(self, x, p, k, f, l):
        """ Constructor"""
        
        self.xi_L = x               # Persistence length of point
        self.pe = p                 # Peclet number of point
        self.k = k                  # Bending rigidty
        self.f = f                  # Pe number       
        self.l = l                  # Aspect ratio of the filament
        
        return

    def analyse_type(self):
        """ Determine the type of the point by analysing average cluster size"""
        
        t = " "                     # Holder string
        
        data_base = '/local/duman/SIMULATIONS/many_polymers_5/density_0.2/kappa_'
        path = data_base + str(self.k) + '/fp_' + str(self.f) + '/CLUSTER/avg_size.txt'
        if os.path.exists(path):
            data = np.loadtxt(path, dtype=float)
        else:
            data = 10.
        self.cs = data
        print path
        print data
        if data < 12.:
            t = "gas"
        elif data > 200.:
            t = "giant"
        else:
            t = "cluster"

        self.type = t               # Type of point
            
        return 
            
       
    def determine_type(self):
        """ Determine the type of the point"""
        
        t = " "                      # Holder string
        self.cs = 0.        
        ## On the lowest rigidities
        if self.xi_L < 0.3:
            if self.pe > 700.:
                t = "spiral"
            elif self.pe < 700. and self.pe > 100.:
                if self.xi_L > 0.1:
                    t = "cluster"
                else:
                    t = "gas"
            else:
                t = "gas"
        else:
            t = "loop"
            
        self.type = t               # Type of point
        
        if self.pe < 7500.:
            if self.pe < 10.:
                self.pe += 1.
            elif self.pe < 700. and self.pe > 100.:  
                self.pe += 50.
            else:
                self.pe += 200.
            
        return 
        
        
    def set_plot_props(self):
        """ Determine plot properties based on the type"""
        
        if self.type == "gas":
            self.marker = "v"
            self.color = "cyan"
            
        elif self.type == "cluster":
            self.marker = "o"
            self.color = "maroon"
            
        elif self.type == "spiral":
            self.marker = "*"
            self.color = "green"
            
        elif self.type == "loop":
            self.marker = "o"
            self.color = "maroon"
            
        elif self.type == "giant":
            self.marker = "s"
            self.color = "red"
        
        return
  

##################################################################

def plot_data(points, long_points, sfolder, xil, pe):
    """ Plot phase diagram"""
      
    ## Latexify 
    fig_width_pt = 255.22124  
    inches_per_pt = 1.0/72.27
    golden_mean = (np.sqrt(5)-1.0)/2.0
    fig_width = fig_width_pt*inches_per_pt
    fig_height = fig_width*golden_mean
    fig_size = [fig_width, fig_height]
#    params = {'backend': 'ps',
#              'axes.labelsize':10,
#              'text.fontsize': 10,
#              'legend.fontsize': 10,
#              'xtick.labelsize': 8,
#              'yick.labelsize': 8,
#              'tex.usetex': True,
#              'figure.figsize': fig_size}
#    pylab.rcParams.update(params)

    fig = plt.figure()
    ax_len = 0.9                          # Length of one subplot square box
    ax_b = 0.1                            # Beginning/offset of the subplot in the box
    ax_sep = 0.1                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots in the x direction  
    tick_num = 4                          # Number of ticks
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x)  
    downlim_x = 1
    uplim_x = 3e+4
    downlim_y = 0.02
    uplim_y = 30
    
    ax = subp.addSubplot()
    
    for point in points:
        ax.loglog(point.pe, point.xi_L, marker=point.marker, c=point.color, markersize=15)
    for point in long_points:
        ax.loglog(point.pe, point.xi_L, marker=point.marker, markerfacecolor="None", markeredgecolor=point.color, markersize=15)
    #ax.axvspan(6000,13000,facecolor='yellow',alpha=0.3)
    #ax.pcolor(pe, xil, )
    ax.set_xlabel(r'Pe',fontsize=40)        
    ax.set_ylabel(r'$\xi_{p} /L$',fontsize=40)
    ax.set_xlim((downlim_x,uplim_x))
    ax.set_ylim((downlim_y,uplim_y))
#    ax.xaxis.set_ticks( np.linspace(downlim, uplim, num=tick_num, endpoint=True) )
#    ax.yaxis.set_ticks( np.linspace(downlim, uplim, num=tick_num, endpoint=True) )
    ax.tick_params(axis='both', which='major', labelsize=40)   
    
    plt.savefig(sfolder + '/point_phase_digram.eps',dpi=200,bbox_inches='tight',pad_inches=0.08)
    plt.savefig(sfolder + '/point_phase_digram.png',dpi=200,bbox_inches='tight',pad_inches=0.08)    
    plt.clf()
    
    return
      
##################################################################

def main():
    """ main function, called when the script is started"""
      
    data_base = '/local/duman/SIMULATIONS/many_polymers_5/'
    save_base = '/usr/users/iff_th2/duman/RolfData/many_polymers_5'
    
    
    ## Index the data
#    density = [0.08, 0.2, 0.4]
#    xi_L = [0.05, 0.1, 0.2, 1.0, 2.5, 5.0, 8.0, 16.0]
#    Pe = [3.0, 7.0, 15.0, 50.0, 150.0, 500.0, 750.0, 1500.0, 4375.0, 8000.0, 10000.0]
#    kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 400.0]
#    fp = [0.0, 0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
    density = [0.2]
    kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 400.0]
    xi_L = [0.05, 0.1, 0.2, 1.0, 2.5, 5.0, 8.0, 16.0]
    fp = [0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    Pe = [3.0, 7.0, 15.0, 50.0, 150.0, 500.0, 750.0, 1500.0, 4375.0]
    
    ## Create points
    points = []
    for i, x in enumerate(xi_L):
        for j, p in enumerate(Pe):
            points.append( Phase(x, p, kappa[i], fp[j], 'short') ) 
          
    for point in points:
        point.analyse_type()
        point.set_plot_props()
     
    long_xil = [0.05, 0.2, 1.0, 2.5, 16.0]
    long_pe = [3.0, 150.0, 750.0, 8000.0, 10000.0]
    long_kappa = [5.0, 20.0, 100.0, 250.0, 1600.0]
    long_fp = [0.0003, 0.015, 0.075, 0.0, 0.0]
    long_points = []
    for i, x in enumerate(long_xil):
        for j, p in enumerate(long_pe):
            long_points.append( Phase(x, p, long_kappa[i], long_fp[j], 'long') ) 
        
    for point in long_points:
        point.determine_type()
        point.set_plot_props()
        
    plot_data(points, long_points, save_base, xi_L, Pe)
    


##################################################################

if __name__ == '__main__':
    main()
    



