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
from scipy.special import gammainc
from scipy.integrate import quad, dblquad, nquad

#
# Function definitions
#
 
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
    

def integrand(x,xp,y,yp,a):
    return np.exp( -np.sqrt( (x-xp)**2 + (y-yp)**2 ) / a )

def bounds(w):
    return [0, w]
    
def func(x, a):
    return nquad(integrand, [bounds(x), bounds(x), bounds(x), bounds(x)], args=(a,)) / x**4
    #x_a = x/a
    #return 4 ( 2 - x_a**2 - np.exp(x_a)*(1-x_a) - np.exp(-x_a)*(1+x_a) ) / x_a**4
    #return 2 * a**2 * (1 - (np.exp(-x/a)*(1+x/a))) / x**2
    #return 2 * x**(-2) * gammainc(2/a, x**a) / a

def powerlaw(x, a, b):
    #return (x**(-a))*np.exp(-x/b)
    return a * x**(-b)

def explaw(x, a, b):
    return a * np.exp(-x/b)

def powerwithexp(x, a, b):
    return (x**(-a))*np.exp(-x/b)

def powerwithupexp(x, a, b, c, d, f):
    return (x**(-a))*np.exp(-x/b) + c*(x**(d))*np.exp(-x/f)

def stretchexp(x, a, b, c):
    return a * np.exp(-(x**b)/c)

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
    # This needs to be changed!!! It needs to be more elegant now, but hey, it works!
    loadSimData(args.datafile+'/density_0.08/kappa_5.0/fp_0.24/init_info.txt')
    
    # Plot properties
    ax_len = 0.6                      # Length of one subplot square box
    ax_b = 0.1                        # Beginning/offset of the subplot in the box
    ax_sep = ax_len/1.3               # Separation length between two subplots
    total_subplots_in_y = 2           # Total number of subplots
    
    # Index the data
#    density = [2.0, 5.0, 10.0, 20.0]
#    kappa = [2.5, 12.5, 62.5, 200.0]
#    fp = [0.0, 0.08, 0.024, 0.24, 0.8, 2.4]
    #density = [0.08, 0.2, 0.4, 0.6, 0.8]
    density = [0.08, 0.2, 0.4]
    kappa = [2.5, 5.0, 25.0, 62.5, 125.0]
    fp = [0.0, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
     
    # Continue indexing the data 
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
    # It will take some time before actually doing though
    # These are some useful but completely unelegant holder variables which help in plotting
    # This needs to change!!! It should be much more portable
    i0prev = density[0]
    i1prev = -1.0
    i0end = density[-1]
    i0cnt = 0
    i0tot = 7

    ax = []
    fig = plt.figure() 
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_y)
    subcnt = -1
    slope_flag = True
    
    for i, grp in group:
        # i are the key values,          [tuple]
        # grp is the grouped element     [DataFrame]
        
        #print i, i[0], i[1], i0prev, i1prev
            
            
        # First constant will be the filename (density in this case) 
        # To save the figure, we first need to plot stuff
        if i0prev != i[0]:
                        
            # Save the previous figure to this file
            path1 = args.savefile + '/plots'
            path2 = args.savefile + '/plots/ORDER_PARAMETER'
            if os.path.exists(path1) == False:
                os.mkdir(path1)
            if os.path.exists(path2) == False:  
                os.mkdir(path2)
                
            plt.savefig(path2 + '/nematic_order_set_2_' + \
                labels[constants[0]] + '_' + "{0:.2f}".format(i0prev) + '.png', dpi=200, bbox_inches='tight', pad_inches=0.08)
                
            plt.clf() 
            
            # Now we plotted the first figure, so we're good, we don't need this branch anymore
            i0prev = i[0]
            
            # but before, let's create a new figure
            ax = []
            fig = plt.figure()
            subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_y)
            subcnt = -1
            

        
        # Second constant will determine the subplots
        if i1prev != i[1]:
            i1prev = i[1]
            ax.append( subp.addSubplot() )
            subcnt += 1
            
    
        # Comparison parameter will determine the legends
        x = []
        y = []
        ystd = []
        yfit = []
        minxprev = 0
        maxxprev = 0
        
        # To color the labels in the legend
        color = iter(plt.cm.rainbow(np.linspace(0,1,len(kappa))))
        
        
        for j, f in enumerate(grp['folder']):
            #pr( j , i
            #print f
        
            #ofile = open(f+'/nematic_corr_len.data','w')
        
            # Load particular simulation data
            loadSimData(f+'/init_info.txt')
                        
            # To get the analysis data
            f += args.analysis+'.data'
            
            if os.path.exists(f):
                order_data = pd.read_csv(f, sep='\t', skiprows=4, header=None) 
                cols = np.arange(1,len(order_data.columns))
                t = order_data[0]
                y = order_data[cols].mean(axis=0)
                ystd = y/np.sqrt(len(t))
                ystd.index = cols-1
                order_data = pd.read_csv(f, sep='\t', skiprows=3, nrows=1, header=None)
                x = np.asarray(order_data[cols].transpose())
                lnx = np.log(x)
                
                # Curve fitting
                xdata = np.ones((16,1))
                ydata = np.ones((16,1))
                xdata2 = np.transpose(np.array(x, dtype=float))
                ydata2 = np.array(y, dtype=float)
                for iii in np.arange(16):
                    xdata[iii] = xdata2[0][iii]
                    ydata[iii] = ydata2[iii]
                xdata = xdata[:,0]
                ydata = ydata[:,0]
                
                popt, pcov = curve_fit(stretchexp, xdata, ydata)
                yfit = stretchexp(xdata, popt[0], popt[1], popt[2])
                    
#                if float(grp['xi_L'].iloc[j]) > 2.0:
#                    if float(i[1]) > 10 and float(i[1]) < 500:
#                        print i[1], '\t', f
#                        popt, pcov = curve_fit(powerwithupexp, xdata, ydata)
#                        #yfit = powerwithexp(xdata, popt[0], popt[1])
#                        yfit = powerwithupexp(xdata, popt[0], popt[1], popt[2], popt[3], popt[4])
#                    else:
#                        popt, pcov = curve_fit(stretchexp, xdata, ydata)
#                        yfit = stretchexp(xdata, popt[0], popt[1], popt[2])
#                else:
#                    popt, pcov = curve_fit(stretchexp, xdata, ydata)
#                    yfit = stretchexp(xdata, popt[0], popt[1], popt[2])
                

                        
#                ofile.write(str(popt[0]))
#                ofile.write('\n')

            #print subcnt

            # Plot the data to the subplot
            c = next(color)
            
            ax[subcnt].errorbar(x, y, yerr = ystd, \
            label=labels[compareby] + "{0:.1f}".format(grp['xi_L'].iloc[j]), linewidth=2., color=c, marker ='o')
            ax[subcnt].plot(x, yfit, '--', \
            label=labels[compareby] + "{0:.1f}".format(grp['xi_L'].iloc[j]) + "_f", linewidth=2., color=c, marker ='o')
            
#            ax[subcnt].plot(x, yfit, '--', \
#            label=labels[compareby] + "{0:.1f}".format(grp['xi_L'].iloc[j]) + "_f", linewidth=1., color=c)
            
            
            if len(x) > 0:
                minx = min(lnx)
                maxx = max(lnx)     
            else:
                minx = 0
                maxx = 0
                
            if minx < minxprev:
                minxprev = minx
            if maxx > maxxprev:
                maxxprev = maxx
    
        minx = minxprev
        maxx = maxxprev
        
        if maxx != 0:
            #ax[subcnt].set_xlim((minx,maxx))
            ax[subcnt].set_ylim((0,1))
            #ax[subcnt].xaxis.set_ticks( np.arange(int(minx),int(maxx),int((maxx-minx)/5)) )
            ax[subcnt].yaxis.set_ticks( np.arange(0.,1.,0.3) )
            ax[subcnt].tick_params(axis='both', which='major', labelsize=12)
        else:
            plt.setp(ax[subcnt].get_xticklabels(),visible=False)
            plt.setp(ax[subcnt].get_yticklabels(),visible=False)
            
            
        # Figure texts
        ax[subcnt].set_title(constants[1] + " = " + "{0:.1f}".format(i[1]), weight='bold', size='large')
        ax[subcnt].set_xscale('log')
        if maxx != 0:
            ax[subcnt].legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0., prop={'size': 12})
        ax[subcnt].set_xlabel('$w [\\sigma]$', weight='bold', size='large')
        ax[subcnt].set_ylabel('$S_{2}^2$', weight='bold', size='large')
            
        if subcnt == i0tot:
            plt.figtext((subp.xbeg+ax_len+ax_sep)/3, subp.ybeg+2*ax_sep, \
            'Nematic Order Parameter, ' + labels[constants[0]] + ' = ' + "{0:.2f}".format(i0prev), size='xx-large', weight='bold')
            
        if i[0] == i0end and subcnt == i0tot:
            # Save the previous figure to this file
            path1 = args.savefile + '/plots'
            path2 = args.savefile + '/plots/ORDER_PARAMETER'
            if os.path.exists(path1) == False:
                os.mkdir(path1)
            if os.path.exists(path2) == False:  
                os.mkdir(path2)
                
            plt.savefig(path2 + '/nematic_order_set_2_' + \
                labels[constants[0]] + '_' + "{0:.2f}".format(i0prev) + '.png', dpi=200, bbox_inches='tight', pad_inches=0.08)
                
            plt.clf()
            
        #print i, i[0], i[1], i0prev, i1prev



if __name__ == '__main__':
    main()



