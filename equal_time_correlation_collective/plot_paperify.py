
# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import pandas as pd
from string import atof
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
    body_length = B*N
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    tau = body_length**3*2.01/4/kT
    tau = 100.
    
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
    kappa = [2.5, 5.0, 25.0, 62.5, 125.0, 200.0]
    fp = [0.24, 7.0]
#    density = [0.08, 0.2]
#    kappa = [2.5, 5.0, 25.0, 62.5, 125.0, 200.0]
#    fp = [0.0, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
    ## Plot properties
    ax_len = 0.8                      # Length of one subplot square box
    ax_b = 0.1                        # Beginning/offset of the subplot in the box
    ax_sep = 0.1                      # Separation length between two subplots
    total_subplots_in_x = 3           # Total number of subplots
    fig = plt.figure()
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x)
    subcnt = -1
    ax = []
    
    ## Labels to be used in plotting 
    labels={'Pe':r"Pe = ",
       'xi_L':r"$\xi_p/L = $",
       'density':'Density'}   
       
    data_base = '/local/duman/SIMULATIONS/many_polymers_5/'
    save_base = '/usr/users/iff_th2/duman/RolfData/many_polymers_5'
    
    for d in density:
        
        ## Load all the data
        ## y as a fnc. of fp
        x = []
        yf = []
        yf_fit = []
        for f in fp:
            
            ax.append( subp.addSubplot() )
            subcnt += 1
                       
            ## y as a fnc. of k
            yk = []
            yk_fit = []
            for k in kappa:
                
                ## Load initial simulation data
                loadSimData(data_base + 'density_' + str(d) + '/kappa_' + \
                    str(k) + '/fp_' + str(f) + '/init_info.txt')     
                    
                ## Load data to be plotted
                path = generate_path(str(d),str(k),str(f),"POLARITY_CORRELATION")
                data = np.transpose(np.loadtxt(path, dtype=float))
                x = data[0]
                y = data[1]
                
                ## Curve fitting
                popt, pcov = curve_fit(exp_law, x, y)
                yfit = exp_law(x, popt[0], popt[1])
                print popt[1]

                    
                ## Save as a function of kappa    
                yk.append(y)
                yk_fit.append(yfit)
            
            ## Save as a function of propulsion force
            yf.append(yk)
            yf_fit.append(yk_fit)
        
        print np.shape(yf)
        
        ## Plot all the data
        for i, f in enumerate(fp):
                
            ## Color code for the plots
            color = iter(plt.cm.rainbow(np.linspace(0,1,len(kappa))))
            
            for j, k in enumerate(kappa):
                
                ## Convert to dimensionless parameters
                xil = convert_to_dimensionless(k, 'kappa')
                
                ## Plot the data
                c = next(color)
                ax[i].semilogx(x,yf[i][j],'o', \
                    label=labels['xi_L']+'{0:.1f}'.format(xil),linewidth=2.,color=c)
                    
                ax[i].semilogx(x,yf_fit[i][j], \
                    label='_nolegend_',linewidth=2.,color=c)                 
                        
            ## Plot properties
            pe = convert_to_dimensionless(f, 'fp')
            #ax[i].set_title(labels['Pe'] + '{0:.0f}'.format(pe), weight='bold', size=30)
            ax[i].set_xlabel('$r/\\sigma$', weight='bold', size=30)
            #ax[i].set_xlim((0,10**4))
            #ax[i].set_ylim((10**-5, 1)) 
            ax[i].xaxis.set_ticks( [10**0, 10**1, 10**2] )
            ax[i].yaxis.set_ticks( [0, 0.5, 1])
            ax[i].tick_params(axis='both', which='major', labelsize=30)
            
        ax[0].set_ylabel('$g_{p}(r)$', weight='bold', size=30)
        
        plt.setp(ax[0].get_yticklabels(),visible=False)
        plt.setp(ax[1].get_yticklabels(),visible=False)

        ax[1].legend(bbox_to_anchor=(0.6,0.,0.45,1.), loc=2, borderaxespad=0., prop={'size': 20}, mode="expand", frameon=False)
        
#        plt.figtext(0.5, 1.1, \
#            'Polarity Correlation, ' + labels['density'] + ' = ' + "{0:.2f}".format(d), size=35, weight='bold')        
        ## Save the previous figure to this file
        path1 = save_base + '/plots'
        path2 = save_base + '/plots/POLARITY_CORRELATION'
        if os.path.exists(path1) == False:  
            os.mkdir(path1)
        if os.path.exists(path2) == False: 
            os.mkdir(path2)                
        plt.savefig(path2 + '/paper_set_' + \
            labels['density'] + '_' + "{0:.2f}".format(d) + '.png', dpi=200, bbox_inches='tight', pad_inches=0.08) 
        plt.savefig(path2 + '/paper_set_' + \
            labels['density'] + '_' + "{0:.2f}".format(d) + '.eps', dpi=200, bbox_inches='tight', pad_inches=0.08)                
        plt.clf()
                


##########################################################################

if __name__ == '__main__':
    main()



