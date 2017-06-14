# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os.path
import glob
import pandas as pd
from string import atof
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
import colorsys 
import sys

#
# Function definitions
#
 
# Check minimum and maximum values
def checkMinMax(mini, minprev, maxi, maxprev):
    if mini < minprev:
        minprev = mini
    if maxi > maxprev:
        maxprev = maxi
        
    return minprev, maxprev
    
    
# Load initial simulation data      
def loadSimData(datafile):
    
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


#
# Class definitions
#

# Arrange subplot grid structure (square box is assumed)
class Subplots:
    
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
   

def main(): 
    
    # Argument parsing (command line options)
    parser = argparse.ArgumentParser()
    parser.add_argument("datafile", help="Folder containing data")
    parser.add_argument("analysis", help="Type of analysis")
    parser.add_argument("savefile", help="Folder in which figures should be saved")
    args = parser.parse_args()
    
    # Load saved preliminary simulation data into relevant variables
    loadSimData(args.datafile+"/density_20.0/kappa_12.5/fp_0.24/init_info.txt")
    
    # Plot properties
    ax_len = 0.3                      # Length of one subplot square box
    ax_b = 0.1                        # Beginning/offset of the subplot in the box
    ax_sep = ax_len/1.3               # Separation length between two subplots
    total_subplots_in_y = 2           # Total number of subplots
    
    # Index the data
    density = [2.0, 5.0, 10.0, 20.0]
    kappa = [2.5, 12.5, 62.5, 200.0]
    fp = [0.0, 0.08, 0.024, 0.24, 0.8, 2.4]
    
    folders = []
    for d in density:
        for k in kappa:
            for f in fp:
                folders.append( args.datafile + '/density_' + str(d) + \
                '/kappa_' + str(k) + '/fp_' + str(f) )
    
    # Format the data            
    files = []
    tmp = {}
    for folder in folders:
        
        tmp = {}
        dname = folder.split('/')[-3]
        ddict = dict(zip(dname.split('_')[::2],map(atof,dname.split('_')[1::2])))
        
        kname = folder.split('/')[-2]
        kdict = dict(zip(kname.split('_')[::2],map(atof,kname.split('_')[1::2])))
        
        fname = folder.split('/')[-1]
        fdict = dict(zip(fname.split('_')[::2],map(atof,fname.split('_')[1::2])))
    
        tmp.update(ddict)
        tmp.update(kdict)
        tmp.update(fdict)
    
        peclet = tmp['fp']*body_length**2/kT
        persistence_over_L = tmp['kappa']/(body_length*kT)    
        
        tmp.update({'folder':folder, 'Pe':peclet, 'xi_L':persistence_over_L})
        files.append(tmp)
            

    Data = pd.DataFrame(files,columns=files[0].keys())
    
    # Choose the comparison parameter for plotting 
    comparison_array = ['density', 'Pe', 'xi_L']
    compareby = 'density'         # make this a command line option
    constants = []
    for element in comparison_array:
        if element != compareby:
            constants.append(element)

    # Group the data accordingly with the chosen comparison parameter
    group = Data.groupby(compareby)
 
    # Labels
    labels={'Pe':r"Pe",
       'xi_L':r"$\xi_p/L$",
       'density':'Density'}
       

    # Plot the data
    ax = []
    fig = plt.figure() 
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_y)
    subcnt = -1
    
    minx = Data[constants[0]].min()
    maxx = Data[constants[0]].max()
    miny = Data[constants[1]].min()
    maxy = Data[constants[1]].max()
    
    minz = 0
    maxz = 0
    minzprev = 1
    maxzprev = 1
    S_av_2 = np.zeros((len(Data['folder'])))
    cnt = -1
    for f in Data['folder']:
        loadSimData(f+'/init_info.txt')
        f += args.analysis+'.data'
        if os.path.exists(f):    
            data = np.loadtxt(f)
            #print data
            cnt += 1
            if np.size(data) != 0:
                S_av_2[cnt] = float(data)
            else:
                S_av_2[cnt] = 0
        minz = min(S_av_2)
        maxz = max(S_av_2)
        minzprev, maxzprev = checkMinMax(minz, minzprev, maxz, maxzprev)
    
    print minx, '\t', maxx, '\t', miny, '\t', maxy
    
    line = 0                
    for i, grp in group:
        # i are the key values,          [tuple]
        # grp is the grouped element     [DataFrame]
             
        #print i, '\t', grp
        
        # Comparison parameter will determine the subplots in the figure
        # (density, in this case)
        ax.append( subp.addSubplot() )
        subcnt += 1
        
        
        # Constants will determine the axes in the subplot (Pe and xi_L)
        # with order parameter being the z dimension (size or coloring)
        S_av = np.zeros((len(grp['folder'])))
        for j, f in enumerate(grp['folder']):
            loadSimData(f+'/init_info.txt')
            f += args.analysis+'.data'
            if os.path.exists(f):
                data = np.loadtxt(f)
                if np.size(data) != 0:
                    S_av[j] = float(data)
        
            #print grp[constants[0]], '\t', grp[constants[1]], '\t', S_av
            sfac = np.pi * ( 0.1 * S_av**2 )
            line = ax[subcnt].scatter(grp[constants[0]], grp[constants[1]], c=S_av, \
            vmin=minzprev, vmax=maxzprev, edgecolors='None', cmap=plt.cm.brg)
            
        ax[subcnt].set_title( labels[compareby] + " = " + "{0:.1f}".format(i), \
        weight='bold', size='large' )
        ax[subcnt].set_xlabel(labels[constants[0]], size='large')
        ax[subcnt].set_ylabel(labels[constants[1]], size='large')
        ax[subcnt].set_xlim((minx-5,maxx+5))
        ax[subcnt].set_ylim((miny-5,maxy+5))
        ax[subcnt].xaxis.set_ticks( np.arange(int(minx-5),int(maxx+5),int((maxx-minx)/4)) )
        ax[subcnt].yaxis.set_ticks( np.arange(int(miny-5),int(maxy+5),int((maxy-miny)/3)) )
        
    cax_1 = plt.axes([subp.xbeg+ax_len+0.04, subp.ybeg-ax_len-ax_sep, 0.01, 2*ax_len+ax_sep])
    plt.colorbar(line, cax=cax_1)
    cax_1.tick_params(width=0.4, labelsize=6)
    cax_1.set_title('$S_{1}$',fontsize=10)
    
    plt.figtext(subp.xbeg-ax_len, subp.ybeg+ax_len+0.1, \
        'Corr. Area of Nematic Order', size='xx-large', weight='bold')
        
        
    path1 = args.savefile + '/plots/phase'
    path2 = args.savefile + '/plots/phase/ORDER_PARAMETER'
    if os.path.exists(path1) == False:
        os.mkdir(path1)
    if os.path.exists(path2) == False:  
        os.mkdir(path2)
        
    plt.savefig(path2 + '/corr_len_nematic_set_' + \
        labels[compareby] + '.png', dpi=200, bbox_inches='tight', pad_inches=0.08)
        
    plt.clf()
        


if __name__ == '__main__':
    main()

