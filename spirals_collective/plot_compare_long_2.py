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
from matplotlib.colors import LinearSegmentedColormap
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


##########################################################################

def save_data(f, x, y, out_fname, d, k, fp):
    """ Save the data for plotting it better in the future"""
    
    ## create the folder in which saving is gonna be performed
    out_fname = '/usr/users/iff_th2/duman/RolfData/GraphData/' + out_fname
    if os.path.exists(out_fname) == False:
        os.mkdir(out_fname)
    
    ## save the data into the designated file
    save_path = out_fname + '/density_' + str(d) + '_kappa_' + str(k) + '_fp_' + str(fp) + '.data'
    print save_path
    save_file = open(save_path, 'w')
    for jj in range(len(x)):
        save_file.write(str(x[jj]) + '\t\t' + str(y[jj]) + '\n')
    save_file.close() 
    
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

#########################################################################
 
def checkMinMax(mini, minprev, maxi, maxprev):
    """ Check minimum and maximum values"""
    if mini < minprev:
        minprev = mini
    if maxi > maxprev:
        maxprev = maxi
        
    return minprev, maxprev
    
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
        self.tot = t        # Total number of subplots  in y direction
        
    # Add a subplot in the grid structure
    def addSubplot(self):
        
        # Increase the number of subplots in the figure
        self.totcnt += 1
        
        # Indices of the subplot in the figure
        self.nx = self.totcnt/(self.tot)
        self.ny = self.totcnt%(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])    
   
##########################################################################
   
def main(): 
    
    database = '/local/duman/SIMULATIONS/many_polymers_5'
    savebase = '/usr/users/iff_th2/duman/RolfData/many_polymers_5'
    path = '/density_0.08/kappa_1.25/fp_7.0/'
    loadSimData(database+path+'init_info.txt')    
    f = database+path+'SPIRALS/spiral_histogram.data'
    order_data = pd.read_csv(f, sep='\t', skiprows=1, header=0) 
    x = (order_data['upper_edge'] + order_data['lower_edge'])/2
    y = order_data['number'] 
    y = np.asarray(y)
    
    #maxx = max(x)
    #y = np.linspace(0,maxx+0.3,num=len(x))
    print x
    print y
    print max(x), max(y)
    y[:20] -= 5.95221000e+05
    y[20:] += 5.95221000e+05
    
    # Normalize distribution
    weights = np.ones_like(y)*M
    y /= weights
    ysum = np.sum(y)
    for el in range(len(y)):
        y[el] = y[el]/0.05/ysum

    savepath = savebase+'/density_0.08/kappa_1.25/fp_16.0/SPIRALS'
    os.system('mkdir -p ' + savepath)
    save_data(f, x, y, 'SPIRALS', 0.08, 1.25, 16.0)      

if __name__ == '__main__':
    main()

