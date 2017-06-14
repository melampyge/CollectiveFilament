# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os.path
import pandas as pd
from string import atof


##########################################################################

#
# Function definitions
#

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
    nbeads = int(L)
    nfil = int(N)
    nmol = int(M)
    tau_D = body_length**2*(N+1)/4/kT
    tau_D /= 10.
    #tau_A = (N+1)/Fmc
    print "Diffusive time scale is ", tau_D
    #print "Advective time scale is ", tau_A
    
    return


##########################################################################
 
def checkMinMax(mini, minprev, maxi, maxprev):
    """ Check minimum and maximum values"""
    if mini < minprev:
        minprev = mini
    if maxi > maxprev:
        maxprev = maxi
        
    return minprev, maxprev

##########################################################################
    
def convert_to_dimensionless(param, ptype):
    """ Convert to dimensionless numbers"""
    
    if ptype == 'kappa':
        return param/kT/body_length
    else:
        return param*body_length**2/kT
    

##########################################################################

def generate_path(dens, kap, fpr, fn):
    """ Generate file path"""
    
    return "/usr/users/iff_th2/duman/RolfData/GraphData/" + fn + "/density_" + \
        dens + "_kappa_" + kap + "_fp_" + fpr + ".data"
    
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
    
def linearFunc(x,a):
    """ Linear function"""
    return a * x
    
##########################################################################    
    
def quadFunc(x,a):
    """ Quadratic function"""
    return a * x**2  

##########################################################################
     
#
# Class definitions
#

##########################################################################

class Subplots:
    """ Arrange subplot grid structure (square box is assumed)"""
    
    totcnt = -1             # Total number of subplots 
    
    # Constructor
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in x direction
        
    # Add a subplot in the grid structure
    def addSubplot(self):
        
        # Increase the number of subplots in the figure
        self.totcnt += 1
        
        # Indices of the subplot in the figure
        self.nx = self.totcnt%(self.tot)
        self.ny = self.totcnt/(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])    
   
##########################################################################
   
def main(): 
    
    ## Index the singled out data
    density = [0.2]
    kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 300.0, 400.0]
    #kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 400.0]
    fp = [0.0, 0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
    ## Plot properties
    ax_len = 0.8                      # Length of one subplot square box
    ax_b = 0.1                        # Beginning/offset of the subplot in the box
    ax_sep = 0.1                      # Separation length between two subplots
    total_subplots_in_x = 3           # Total number of subplots
    fig = plt.figure()
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x)
    ax = []
    ax.append( subp.addSubplot() )
       
    ## Labels to be used in plotting 
    labels={'Pe':r"Pe = ",
       'xi_L':r"$\xi_p/L = $",
       'density':'Density'}   
       
    data_base = '/local/duman/SIMULATIONS/many_polymers_5/'
    save_base = '/usr/users/iff_th2/duman/RolfData/many_polymers_5'
    data_path = data_base + 'density_0.08/kappa_5.0/fp_0.024'
    loadSimData(data_path + '/init_info.txt')     
    
    x = []
    for f in fp:
        pe = convert_to_dimensionless(f, 'fp')
        x.append(pe)     
        
    for d in density:
        
        ## Load all the data
        ## y as a fnc. of fp
        yf = []
        for f in fp:
        
            ## y as a fnc. of k
            yk = []
            for k in kappa:
                
                ## Load initial simulation data
                data_path = data_base + 'density_' + str(d) + '/kappa_' + \
                    str(k) + '/fp_' + str(f)
                loadSimData(data_path + '/init_info.txt')     
                    
                ## Load data to be plotted
                file_path = data_path + '/CLUSTER/avg_lifetime.txt'
                if os.path.exists(file_path) == False:    
                    file_path = data_path + '/CLUSTER/avg_lifetime.txt.gz'
                if os.path.exists(file_path) == False:
                    #data += 0.2
                    data = 0.0
                else:
                    data = np.loadtxt(file_path, dtype=float)
#                if f == 0.08 and k == 400.0:
#                    data += 1.6
#                elif f == 0.24 and k == 400.0:
#                    data += 0.9
#                elif f == 0.24 and k == 25.0:
#                    data += 0.2
#                elif f == 0.0112 and k == 5.0:
#                    data -= 0.09
                print k, f, data    
                yk.append(data)
 
            ## Save as a function of propulsion force
            yf.append(yk)
    
        yf = np.transpose(yf)
        
        print np.shape(x), np.shape(yf)
        
        first_time = True 
        ## Plot all the data
        for i, f in enumerate(fp):
                
            ## Color code for the plots
            color = iter(plt.cm.rainbow(np.linspace(0,1,len(kappa))))
            
            if first_time:
                labelon = True
                first_time = False
            else:
                labelon = False

            
            for j, k in enumerate(kappa):
                
                ## Convert to dimensionless parameters
                xil = convert_to_dimensionless(k, 'kappa')
                c = next(color)
                                
                ## Plot the data
                if labelon:      
                    ax[0].loglog(x,yf[j]*50./tau_D,'o', \
                        label=labels['xi_L'] + '{0:.2f}'.format(xil),linewidth=2.,color=c,markersize=13)
                        
                    ax[0].loglog(x,yf[j]*50./tau_D, '--', \
                        label='_nolegend_',linewidth=2.,color=c)    
                        
                else:
                    ax[0].loglog(x,yf[j]*50./tau_D,'o', \
                        label='_nolegend_',linewidth=2.,color=c,markersize=15)
                    ax[0].loglog(x,yf[j]*50./tau_D, '--', \
                        label='_nolegend_',linewidth=2.,color=c)                         
                        
        ## Plot properties
        pe = convert_to_dimensionless(f, 'fp')
        #ax[0].set_title(r'Rotational force of clusters', weight='bold', size=20)
        ax[0].set_xlabel(r'$Pe$', weight='bold', size=40)
        #ax[0].set_xlim((-100,4600))
        #ax[0].set_ylim((4*10**-2, 10**0)) 
        #ax[0].xaxis.set_ticks( [0, 2000, 4000] )
        #ax[0].yaxis.set_ticks( [0.5, 2.0, 3.5])
        ax[0].tick_params(axis='both', which='major', labelsize=40)
        ax[0].set_ylabel(r'$\Delta t/\tau_{D}$', weight='bold', size=40)

        #ax[0].legend(bbox_to_anchor=(0.71,0.,0.45,1.), loc=2, borderaxespad=0., prop={'size': 14}, mode="expand", frameon=False)
        
#        plt.figtext(0.05, 1.1, \
#            'Rotational Force, ' + labels['density'] + ' = ' + "{0:.2f}".format(d), size=35, weight='bold') 
            
        ## Save the previous figure to this file
        path1 = save_base + '/plots'
        path2 = save_base + '/plots/phase/CLUSTER'
        if os.path.exists(path1) == False:  
            os.mkdir(path1)
        if os.path.exists(path2) == False: 
            os.mkdir(path2)                
        plt.savefig(path2 + '/paper_LIFETIME_set_' + \
            '.png', dpi=200, bbox_inches='tight', pad_inches=0.08) 
        plt.savefig(path2 + '/paper_LIFETIME_set_' + \
            '.eps', dpi=200, bbox_inches='tight', pad_inches=0.08)               
        plt.clf()
                
        


if __name__ == '__main__':
    main()

