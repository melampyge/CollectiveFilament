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
from matplotlib.colors import LinearSegmentedColormap
import colorsys 
import sys
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

def save_data(f, x, y, ystd, out_fname):
    """ Save the data for plotting it better in the future"""
    
    ## create the folder in which saving is gonna be performed
    out_fname = '/usr/users/iff_th2/duman/RolfData/GraphData/' + out_fname
    if os.path.exists(out_fname) == False:
        os.mkdir(out_fname)
        
    ## parse the filename
    fp = f.split('/')[-3]
    k = f.split('/')[-4]
    d = f.split('/')[-5]
    
    ## save the data into the designated file
    save_path = out_fname + '/' + d + '_' + k + '_' + fp + '.data'
    save_file = open(save_path, 'w')
    for jj in range(len(x)):
        save_file.write(str(x[jj]) + '\t\t' + str(y[jj]) + '\t\t' + str(ystd[jj]) +'\n')
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

def getE2E(f):
    """ Get end-to-end distance vector"""
    
    f += "/RE2E/e2e.data"
    if os.path.exists(f):
        datafile = open(f,"r")
        datafile.readline()
        datafile.readline()
        line = datafile.readline()
        A = line.split()
        
        return float(A[-1])
        
    else:
        return 0.

##########################################################################

def setAnalytic(e2e):
    """ Set the necessary variables used in the calculation of theoretical MSD"""
    
    global Dt, v0
    
    gam = 1
    gam_l = gam*N/body_length
    Dt = kT/(gam_l*body_length)
    #Dr = 9*kT/(body_length**3*gam_l) + Fmc/(gam_l*persistence)
    v0 = Fmc*e2e/(gam_l*body_length)
    print v0, '\t', Fmc, '\t', e2e, '\n'

##########################################################################

def analyticFunc(t,a):
    """ Theoretical MSD"""
    return 4*Dt*t + 2*(v0**2)*(a*t + np.exp(-a*t) - 1)/a**2
    
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
    
    # Argument parsing (command line options)
    parser = argparse.ArgumentParser()
    parser.add_argument("datafile", help="Folder containing data")
    parser.add_argument("analysis", help="Type of analysis")
    parser.add_argument("savefile", help="Folder in which figures should be saved")
    args = parser.parse_args()
    
    # Load saved preliminary simulation data into relevant variables
    # This needs to be changed!!! It needs to be more elegant now, but hey, it works!
    loadSimData(args.datafile+'/density_0.08/kappa_5.0/fp_0.24/init_info.txt')
    
    # Plot properties
    ax_len = 0.7                      # Length of one subplot square box
    ax_b = 0.1                        # Beginning/offset of the subplot in the box
    ax_sep = 0.4                      # Separation length between two subplots
    total_subplots_in_x = 4           # Total number of subplots in the x direction
    
    # Index the data
    density = [0.08, 0.2, 0.4]
    kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 400.0]
    fp = [0.0, 0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
    # Format the data
    Data = format_data(args.datafile, density, kappa, fp)    
    
    # Choose the comparison parameter for plotting 
    comparison_array = ['density', 'Pe', 'xi_L']
    compareby = 'xi_L'         # make this a command line option!!!
    constants = []
    for element in comparison_array:
        if element != compareby:
            constants.append(element)
            
    # What is this comparison?
    # compareby is the legends of the subplots (xi_l)
    # first constant is the file name (density)
    # second constant is the subplot (Pe)

    # Group the data accordingly with the chosen comparison parameter
    group = Data.groupby(constants)
    
    # Labels to be used in plotting 
    labels={'Pe':r"Pe = ",
       'xi_L':r"$\xi_p/L = $",
       'density':'Density'}   

    # Plot the data
    # It will take some time before actually doing so though
    # These are some useful but completely unelegant holder variables which help in plotting
    # This needs to change!!! It should be much more portable
    i0prev = density[0]
    i1prev = -1.0
    i0end = density[-1]
    i0tot = 7

    ax = []
    fig = plt.figure() 
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x)
    subcnt = -1    
    
    
    for i, grp in group:
        # i are the key values,          [tuple]       (density, Pe)
        # grp is the grouped element     [DataFrame] 
                    
        # First constant will be the filename (density, in this case) 
        # To save the figure, we first need to plot stuff, hence this part 
        # will be executed later
        if i0prev != i[0]:
            
            ax[subcnt].legend(bbox_to_anchor=(1.05,0.,0.45,1.), loc=2, borderaxespad=0., prop={'size': 20}, mode="expand")
            
            # Save the previous figure to this file
            path1 = args.savefile + '/plots'
            path2 = args.savefile + '/plots/MSD'
            if os.path.exists(path1) == False: 
                os.mkdir(path1)
            if os.path.exists(path2) == False: 
                os.mkdir(path2)
                
            plt.savefig(path2 + '/set_' + \
                labels[constants[0]] + '_' + "{0:.2f}".format(i0prev) + '.png', dpi=200, bbox_inches='tight', pad_inches=0.08)
                
            plt.clf() 
            
            # Now we plotted the first figure, so we're good, we don't need this branch anymore
            i0prev = i[0]
            
            # but before, let's create a new figure
            ax = []
            fig = plt.figure()
            subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x)
            subcnt = -1
            

        
        # Second constant will determine the subplots (Pe, in this case)
        # when the second constant changes, a new subplot is added
        if i1prev != i[1]:
            i1prev = i[1]
            ax.append( subp.addSubplot() )
            subcnt += 1         # Counts the subplots
            
    
        # Comparison parameter will determine the legends (xi_L, in this case)
        x = []
        y = []
        minxprev = 0
        maxxprev = 0
        minyprev = 0
        maxyprev = 0
        slope_flag = True
        
        # To color the labels in the legend
        color = iter(plt.cm.rainbow(np.linspace(0,1,len(kappa))))

        for j, f in enumerate(grp['folder']):        
            # j is index number
            # f is folder 

            # Write the rotational diffusion for the given folder
            ofile = open(f+'/rot_dif_from_msd.data','w')
            
            # Load particular simulation data for the given folder
            loadSimData(f+'/init_info.txt')
            
            # Get the average end to end vector which will be used in MSD calculation
            e2e = getE2E(f)
            
            # To get the analysis data
            f += args.analysis+'.data'

            # Set the necessary variables for the calculation of analytic MSD
            setAnalytic(e2e)
            
            # Read in the data from the corresponding file
            if os.path.exists(f):
                order_data = pd.read_csv(f, sep='\t', skiprows=1, header=0) 
                x = order_data['Timestep']*dt
                y = order_data['MSD']
                
                # Analytic expression
                popt, pcov = curve_fit(analyticFunc, x, y)
                Dr = popt[0]
                ofile.write(str(Dr))
                ofile.write('\n')
                z = analyticFunc(x, Dr) 
                #z = 4*Dt*x + 2*v0**2*(Dr*x + np.exp(-Dr*x) - 1)/Dr**2
                 
                save_data(f, x, y, z, 'MSD')
 
            # Plot the data to the subplot
            c = next(color)
            ax[subcnt].loglog(x,z/body_length**2, label="_nolegend_", color=c)
            ax[subcnt].loglog(x,y/body_length**2, 'o', label=labels[compareby] + "{0:.1f}".format(grp['xi_L'].iloc[j]), linewidth=2., color=c)
            
            # Plot slopes of 1 and 2 for each subplot once
            if slope_flag:
                slope_flag = False
                popt, pcov = curve_fit(linearFunc, x[30:40]/5, y[30:40]/5/body_length**2)
                ylinfit = linearFunc(x[30:40], popt[0])
                popt, pcov = curve_fit(quadFunc, x[:10], y[:10]/body_length**2)
                yquadfit = quadFunc(x[:10], popt[0])
                ax[subcnt].loglog(x[30:40]/5,ylinfit, '--', label='_nolegend_', linewidth=2., color='gray')
                ax[subcnt].loglog(x[:10],yquadfit, '--', label='_nolegend_', linewidth=2., color='gray')

            if len(x) > 0:
                minx = min(x)
                maxx = max(x)     
            else:
                minx = 0
                maxx = 0
                
            if minx < minxprev:
                minxprev = minx
            if maxx > maxxprev:
                maxxprev = maxx
    
        minx = minxprev
        maxx = maxxprev
        
        slope_flag = True
        
        if maxx != 0:
            #ax[subcnt].set_xlim((minx,maxx))
            #ax[subcnt].set_ylim((0,1))
            #ax[subcnt].xaxis.set_ticks( np.arange(int(minx),int(maxx),int((maxx-minx)/5)) )
            #ax[subcnt].yaxis.set_ticks( np.arange(0.,1.,0.3) )
            ax[subcnt].tick_params(axis='both', which='major', labelsize=25)
        else:
            plt.setp(ax[subcnt].get_xticklabels(),visible=False)
            plt.setp(ax[subcnt].get_yticklabels(),visible=False)
            
            
        # Figure texts
        ax[subcnt].set_title(constants[1] + " = " + "{0:.1f}".format(i[1]), weight='bold', size=30)
        #ax[subcnt].legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0., prop={'size': 12})
        ax[subcnt].set_xlabel('$\\Delta t [t. u.]$', weight='bold', size='large')
        ax[subcnt].set_ylabel('$\\Delta r^2/L^2$', weight='bold', size='large')
            
        if subcnt == i0tot:
            plt.figtext((subp.xbeg+ax_len+ax_sep)/3, subp.ybeg+ax_len+0.2, \
            'Mean Square Displacement, ' + labels[constants[0]] + ' = ' + "{0:.2f}".format(i0prev), size=35, weight='bold')
            
        if i[0] == i0end and subcnt == i0tot:
            
            ax[subcnt].legend(bbox_to_anchor=(1.05,0.,0.45,1.), loc=2, borderaxespad=0., prop={'size': 20}, mode="expand")
            
            # Save the previous figure to this file
            path1 = args.savefile + '/plots'
            path2 = args.savefile + '/plots/MSD'
            if os.path.exists(path1) == False:  
                os.mkdir(path1)
            if os.path.exists(path2) == False: 
                os.mkdir(path2)
                
            plt.savefig(path2 + '/set_' + \
                labels[constants[0]] + '_' + "{0:.2f}".format(i0prev) + '.png', dpi=200, bbox_inches='tight', pad_inches=0.08)
                
            plt.clf()
            
        #print i, i[0], i[1], i0prev, i1prev



if __name__ == '__main__':
    main()
