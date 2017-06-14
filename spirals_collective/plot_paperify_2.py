# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import pandas as pd
from string import atof
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
from scipy.optimize import curve_fit

##########################################################################

#
# Function definitions
#

##########################################################################
 
def loadSimData(datafile):
    """ Load initial simulation data"""
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
    dtSamp, T, box_area, nt, body_length, Pe, persistence 

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
    body_length = B*(N-1)
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)

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
    density = [0.08]
    kappa = [1.25, 2.5, 5.0, 62.5]
    fp = [0.0, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
#    density = [0.08, 0.2]
#    kappa = [2.5, 5.0, 25.0, 62.5, 125.0, 200.0]
#    fp = [0.0, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
    ## Plot properties
    ax_len = 0.9                      # Length of one subplot square box
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
        yf05 = []
        yf10 = []
        yf15 = []
        yf20 = []     
        for f in fp:
        
            ## y as a fnc. of k
            yk05 = []
            yk10 = []
            yk15 = []
            yk20 = []   
            for k in kappa:
                
                ## Load initial simulation data
                data_path = data_base + 'density_' + str(d) + '/kappa_' + \
                    str(k) + '/fp_' + str(f)
                loadSimData(data_path + '/init_info.txt')     
                    
                ## Load data to be plotted
                file_path = data_path + '/SPIRALS/spiral_evolution.data'
                order_data = pd.read_csv(file_path, sep='\t', skiprows=1, header=0) 
                y05 = order_data['s>=0.5'] 
                y10 = order_data['s>=1.0'] 
                y15 = order_data['s>=1.5'] 
                y20 = order_data['s>=2.0'] 
                
                y05 = y05.mean() 
                y10 = y10.mean()
                y15 = y15.mean()
                y20 = y20.mean()
                
                    
                ## Save as a function of kappa    
                yk05.append(y05)
                yk10.append(y10)
                yk15.append(y15)
                yk20.append(y20)
            
            ## Save as a function of propulsion force
            yf05.append(yk05)
            yf10.append(yk10)
            yf15.append(yk15)
            yf20.append(yk20)
        
        yf05 = np.transpose(yf05)
        yf10 = np.transpose(yf10)
        yf15 = np.transpose(yf15)
        yf20 = np.transpose(yf20)
        
        print np.shape(x), np.shape(yf05)
        
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
                    ax[0].loglog(x,yf20[j],'o', \
                        label=labels['xi_L']+'{0:.2f}'.format(xil),linewidth=3.,color=c,markersize=10)
#                    ax[0].loglog(x,yf10[j],'v', \
#                        label=labels['xi_L']+'{0:.1f}'.format(xil),linewidth=2.,color=c)                        
#                    ax[0].loglog(x,yf15[j],'^', \
#                        label=labels['xi_L']+'{0:.1f}'.format(xil),linewidth=2.,color=c)
#                    ax[0].loglog(x,yf20[j],'s', \
#                        label=labels['xi_L']+'{0:.1f}'.format(xil),linewidth=2.,color=c)                        
                else:
                    ax[0].loglog(x,yf20[j],'o', \
                        label='_nolegend_',linewidth=3.,color=c,markersize=10)   
#                    ax[0].loglog(x,yf10[j],'v', \
#                        label=labels['xi_L']+'{0:.1f}'.format(xil),linewidth=2.,color=c)         
#                    ax[0].loglog(x,yf15[j],'^', \
#                        label=labels['xi_L']+'{0:.1f}'.format(xil),linewidth=2.,color=c)                        
#                    ax[0].loglog(x,yf20[j],'s', \
#                        label=labels['xi_L']+'{0:.1f}'.format(xil),linewidth=2.,color=c)
                        
        ## Plot properties
        pe = convert_to_dimensionless(f, 'fp')
        #ax[0].set_title(r'Frac. of filaments w/ s > 2.0', weight='bold', size=20)
        ax[0].set_xlabel(r'$Pe$', weight='bold', size=40)
#            ax[i].set_xlim((-5,10**4))
#            ax[i].set_ylim((10**-5, 1)) 
        #ax[0].xaxis.set_ticks( [10**0, 10**2, 10**4] )
        ax[0].yaxis.set_ticks( [10**-8, 10**-6, 10**-4, 10**-2])
        ax[0].tick_params(axis='both', which='major', labelsize=40)
        ax[0].set_ylabel(r'$P(|s|)$', weight='bold', size=40)
        ax[0].legend(bbox_to_anchor=(0.02,0.,0.45,1.), loc=2, borderaxespad=0., prop={'size': 20}, mode="expand", frameon=False)
        
#        plt.figtext(0.05, 1.1, \
#            'Spiral Number, ' + labels['density'] + ' = ' + "{0:.1f}".format(d), size=35, weight='bold')        
        ## Save the previous figure to this file
        path1 = save_base + '/plots'
        path2 = save_base + '/plots/SPIRALS'
        if os.path.exists(path1) == False:  
            os.mkdir(path1)
        if os.path.exists(path2) == False: 
            os.mkdir(path2)                
        plt.savefig(path2 + '/paper_fraction_20_set_' + \
            labels['density'] + '_' + "{0:.2f}".format(d) + '.png', dpi=200, bbox_inches='tight', pad_inches=0.08)  
        plt.savefig(path2 + '/paper_fraction_20_set_' + \
            labels['density'] + '_' + "{0:.2f}".format(d) + '.eps', dpi=200, bbox_inches='tight', pad_inches=0.08)               
        plt.clf()
                


##########################################################################

if __name__ == '__main__':
    main()



