# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os.path
import pandas as pd
from string import atof
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png


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
        return param*25.**2/1
    

##########################################################################

def generate_path(dens, kap, fpr, fn, fn2):
    """ Generate file path"""
    
    return fn + "/density_" + \
        dens + "/kappa_" + kap + "/fp_" + fpr + "/CLUSTER/" + fn2 + ".txt"
    
##########################################################################

def format_data(datafile, density, kappa, fp):
    """ Format the data"""

    # Continue indexing the data 
    folders = []
    for d in density:
        for k in kappa:
            for f in fp:
                folders.append( datafile + '/density_' + str(d) + \
                '/kappa_' + str(k) + '/fp_' + str(f) )
    
    # Format the data            
    files = []
    tmp = {}
    for folder in folders:
        
        # Just a temporary dictionary 
        tmp = {}
        
        # Formatting the data, get the density
        dname = folder.split('/')[-3]
        ddict = dict(zip(dname.split('_')[::2],map(atof,dname.split('_')[1::2])))
        
        # Get the bending rigidity
        kname = folder.split('/')[-2]
        kdict = dict(zip(kname.split('_')[::2],map(atof,kname.split('_')[1::2])))
        
        # Get the propulsion strength
        fname = folder.split('/')[-1]
        fdict = dict(zip(fname.split('_')[::2],map(atof,fname.split('_')[1::2])))
    
        tmp.update(ddict)
        tmp.update(kdict)
        tmp.update(fdict)
    
        peclet = tmp['fp']*body_length**2/kT
        persistence_over_L = tmp['kappa']/(body_length*kT)    
        
        # Index the formatted data through the dictionary
        tmp.update({'folder':folder, 'Pe':peclet, 'xi_L':persistence_over_L})
        files.append(tmp)
            

    # Load the data into Pandas to do cool stuff
    # Format is 'folder', 'Pe', 'xi_L'
    Data = pd.DataFrame(files,columns=files[0].keys())

    return Data    
    
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
    density = [0.08, 0.2, 0.4]
    kappa = [400.0]
    fp = [0.0, 0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
#    density = [0.08, 0.2]
#    kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 400.0]
#    fp = [0.0, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
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
    labels={'Pe':r"Pe",
       'xi_L':r"$\xi_p/L = $",
       'density':'Density'}   
       
    data_base = '/local/duman/SIMULATIONS/many_polymers_5/'
    save_base = '/usr/users/iff_th2/duman/RolfData/many_polymers_5'
    
    x = []
    for f in fp:
        pe = convert_to_dimensionless(f, 'fp')
        x.append(pe)
        
    for k in kappa:
        
        ## Load all the data
        ## y as a fnc. of k
        yk = []
        for d in density:
                       
            ## y as a fnc. of d
            yd = []
            
            for f in fp:
                
                ## Load initial simulation data
                loadSimData(data_base + 'density_' + str(d) + '/kappa_' + \
                    str(k) + '/fp_' + str(f) + '/init_info.txt')     
                    
                ## Load data to be plotted
                path = generate_path(str(d),str(k),str(f),data_base,"avg_size")
                #print path
                data = np.transpose(np.loadtxt(path, dtype=float))
                print data

                ## Save as a function of density
                yd.append(data)
            
            ## Save as a function of kappa
            yk.append(yd)
        
        print np.shape(x), np.shape(yk)
        print body_length
        ## Plot all the data
        for i, k in enumerate(kappa):
                
            ## Color code for the plots
            color = iter(plt.cm.rainbow(np.linspace(0,1,len(density))))
            
            for j, d in enumerate(density):
                
                ## Convert to dimensionless parameters
                pe = convert_to_dimensionless(f, 'fp')

                ## Plot the data
                c = next(color)
                print d,x,yk[i]
                ax[i].semilogx(x,yk[j],'o', \
                    label='_nolegend_',linewidth=2.,color=c)
                    
                ax[i].semilogx(x,yk[j], \
                    label='$\\phi = $'+'{0:.2f}'.format(d),linewidth=2.,color=c)    
                    
            ## Plot properties
            xil = convert_to_dimensionless(k, 'kappa')
            ax[i].set_title(labels['xi_L'] + '{0:.1f}'.format(xil), weight='bold', size=20)
            ax[i].set_xlabel(labels['Pe'], weight='bold', size=20)
            #ax[i].set_xlim((-5,10**4))
            #ax[i].set_ylim((10**-5, 1)) 
            #ax[i].xaxis.set_ticks( [10**0, 10**2, 10**4] )
            #ax[i].yaxis.set_ticks( [10**-5, 10**-3, 10**0])
            ax[i].tick_params(axis='both', which='major', labelsize=20)
            
            ax[0].set_ylabel('Size [filaments]', weight='bold', size=20)
    
            ax[0].legend(bbox_to_anchor=(0.,0.,0.45,1.), loc=2, borderaxespad=0., prop={'size': 15}, mode="expand", frameon=False)
            
#            plt.figtext(0., 1.1, \
#                'Average Cluster Size, ' + labels['xi_L'] + "{0:.1f}".format(xil), size=35, weight='bold')        
            ## Save the previous figure to this file
            path1 = save_base + '/plots'
            path2 = save_base + '/plots/phase/CLUSTER'
            if os.path.exists(path1) == False:  
                os.mkdir(path1)
            if os.path.exists(path2) == False: 
                os.mkdir(path2)                
            plt.savefig(path2 + '/SIZE_single_set_' + \
                '.png', dpi=200, bbox_inches='tight', pad_inches=0.08)                
            plt.clf()
                
        


if __name__ == '__main__':
    main()

