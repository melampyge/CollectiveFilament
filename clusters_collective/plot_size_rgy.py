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
    body_length = B*N
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)

##########################################################################

def save_data(f, x, y, out_fname):
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

def powerlaw(x, a, b):
    return a * x**b

##########################################################################

def powerlaw_less(x, a):
    return x**a

##########################################################################

def powerwithexp(x, a, b):
    return (x**(-a))*np.exp(-x/b)

##########################################################################

def powerwithupexp(x, a, b, c, d, f):
    return (x**(-a))*np.exp(-x/b) + c*(x**(d))*np.exp(-x/f)

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
    ## - row major -
    
    totcnt = -1             # Total number of subplots 
    
    # Constructor
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction
        
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
    loadSimData(args.datafile+'/density_0.4/kappa_5.0/fp_0.24/init_info.txt')
    
    # Plot properties
    ax_len = 0.7                      # Length of one subplot square box
    ax_b = 0.1                        # Beginning/offset of the subplot in the box
    ax_sep = 0.4                      # Separation length between two subplots
    total_subplots_in_x = 4           # Total number of subplots
    
    # Index the data
    density = [0.2, 0.4]
    kappa = [2.5, 5.0, 25.0, 62.5, 125.0, 200.0]
    fp = [0.0, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
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
    
    y = np.linspace(10, M+1, num=M+1-10, endpoint=True)
    for i, grp in group:
        # i are the key values,          [tuple]       (density, Pe)
        # grp is the grouped element     [DataFrame] 
                    
        # First constant will be the filename (density, in this case) 
        # To save the figure, we first need to plot stuff, hence this part 
        # will be executed later
        if i0prev != i[0]:
            
            ax[subcnt].legend(bbox_to_anchor=(1.05,0.,0.53,1.), loc=2, borderaxespad=0., prop={'size': 25}, mode="expand")
            
            # Save the previous figure to this file
            path1 = args.savefile + '/plots'
            path2 = args.savefile + '/plots/CLUSTER'
            if os.path.exists(path1) == False: 
                os.mkdir(path1)
            if os.path.exists(path2) == False: 
                os.mkdir(path2)
                
            plt.savefig(path2 + '/size_rgy_set_' + \
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
        yfit = []
        minxprev = 0
        maxxprev = 0
        minyprev = 0
        maxyprev = 0
        slope_flag = True
        
        # To color the labels in the legend
        color = iter(plt.cm.rainbow(np.linspace(0,1,len(kappa))))
        
        # j is index
        # f is folder
        for j, f in enumerate(grp['folder']):
            
            # Load the preliminary data for the given folder
            loadSimData(f+'/init_info.txt')
            
            # Path for the data
            f += args.analysis+'.txt'
            
            if os.path.exists(f):
                
                x = np.loadtxt(f, dtype=float)
                x = x[10:]
                
                # Curve fitting
                popt, pcov = curve_fit(powerlaw, x, y)
                yfit = powerlaw(x, popt[0], popt[1])
                print popt[1]
                
                #save_data(f, x, y, 'CLUSTER')
                
#                # Curve fitting
#                print f
#                print i[0], float(i[0])
#                if float(i[0]) == 0.08:
#                    # First two kappa values
#                    if float(grp['xi_L'].iloc[j]) < 0.5:
#                        
#                        popt, pcov = curve_fit(explaw, x, y)
#                        yfit = explaw(x, popt[0], popt[1])
#                    
#                    # Middle kappa value
#                    elif float(grp['xi_L'].iloc[j]) > 0.5 and float(grp['xi_L'].iloc[j]) < 2.0:
#                        if float(i[1]) != 52.02:
#                            if float(i[1]) == 0.0:
#                                popt, pcov = curve_fit(explaw, x, y)
#                                yfit = explaw(x, popt[0], popt[1])
#                            elif float(i[1]) == 15.606:
#                                popt, pcov = curve_fit(powerwithexp, x, y)
#                                yfit = powerwithexp(x, popt[0], popt[1])
#                            elif float(i[1]) > 1562:
#                                popt, pcov = curve_fit(explaw, x, y)
#                                yfit = explaw(x, popt[0], popt[1])                            
#                            else:
#                                popt, pcov = curve_fit(powerwithupexp, x, y)
#                                yfit = powerwithupexp(x, popt[0], popt[1], popt[2], popt[3], popt[4])   
#                            
#                    # Last three kappa values    
#                    else:
#                        if float(i[1]) == 0.0:
#                            popt, pcov = curve_fit(explaw, x, y)
#                            yfit = explaw(x, popt[0], popt[1])   
##                        elif float(i[1]) > 1562:
##                            popt, pcov = curve_fit(powerwithexp, x, y)
##                            yfit = powerwithexp(x, popt[0], popt[1])
#                        else:
#                            popt, pcov = curve_fit(powerwithupexp, x, y)
#                            yfit = powerwithupexp(x, popt[0], popt[1], popt[2], popt[3], popt[4]) 
#                                        
#                elif float(i[0]) == 0.2:
#                    print 'DENSITY CHANGE to 0.2'
#                    # First two kappa values
#                    if float(grp['xi_L'].iloc[j]) < 0.5:
#                        if float(i[1]) == 0.0:
#                            popt, pcov = curve_fit(explaw, x, y)
#                            yfit = explaw(x, popt[0], popt[1])                            
#                        else:
#                            popt, pcov = curve_fit(powerwithupexp, x, y)
#                            yfit = powerwithupexp(x, popt[0], popt[1], popt[2], popt[3], popt[4])        
#                            
#                    else:
#                        if float(i[1]) == 0.0:
#                            popt, pcov = curve_fit(explaw, x, y)
#                            yfit = explaw(x, popt[0], popt[1])  
#                            
#                        elif float(i[1]) > 1562:
#                            if float(grp['xi_L'].iloc[j]) > 4:
#                                popt, pcov = curve_fit(powerlaw, x, y)
#                                yfit = powerlaw(x, popt[0], popt[1])
#                            else:
#                                popt, pcov = curve_fit(powerwithupexp, x, y)
#                                yfit = powerwithupexp(x, popt[0], popt[1], popt[2], popt[3], popt[4])                     
#                        else:
#                            popt, pcov = curve_fit(powerlaw, x, y)
#                            yfit = powerlaw(x, popt[0], popt[1])
                    
   
 
            # Plot the data to the subplot
            c = next(color)
            print np.shape(x), np.shape(y)
            ax[subcnt].loglog(x, y, 'o', \
                label=labels[compareby] + "{0:.1f}".format(grp['xi_L'].iloc[j]), linewidth=1.0, color=c)
            
            ax[subcnt].loglog(x, yfit, \
                label='_nolegend_', linewidth=1.0, color=c)
            
#            # Plot slopes of 1 and 2 for each subplot once
#            if slope_flag:
#                slope_flag = False
#                popt, pcov = curve_fit(linearFunc, x[5:15], y[5:15])
#                ylinfit = linearFunc(x[5:15], popt[0])
#                popt, pcov = curve_fit(quadFunc, x[300:600], y[300:600]/body_length**2)
#                yquadfit = quadFunc(x[300:600], popt[0])
#                ax[subcnt].loglog(x[5:15]/5,ylinfit, '--', label='_nolegend_', linewidth=2., color='gray')
#                ax[subcnt].loglog(x[300:600],yquadfit, '--', label='_nolegend_', linewidth=2., color='gray')
                
            if len(x) > 0:
                minx = min(x)
                maxx = max(x) 
                miny = min(y)
                maxy = max(y)
            else:
                minx = 0
                maxx = 0
            
            minxprev, maxxprev = checkMinMax(minx, minxprev, maxx, maxxprev)
            minyprev, maxyprev = checkMinMax(miny, minyprev, maxy, maxyprev)
    
        minx = minxprev
        maxx = maxxprev
        miny = minyprev
        maxy = maxyprev
        
        if abs(maxx) != 0:
            #ax[subcnt].set_xlim((-5,M+2))
            #ax[subcnt].set_ylim((10**-5, 1))
            #ax[subcnt].xaxis.set_ticks( np.arange(int(minx),int(maxx),int((maxx-minx)/5)) )
            #ax[subcnt].yaxis.set_ticks( np.arange(0.,1.,0.3) )
            ax[subcnt].tick_params(axis='both', which='major', labelsize=25)
        else:
            plt.setp(ax[subcnt].get_xticklabels(),visible=False)
            plt.setp(ax[subcnt].get_yticklabels(),visible=False)
            
            
        # Figure texts
        ax[subcnt].set_title(constants[1] + " = " + "{0:.0f}".format(i[1]), weight='bold', size=30)
        ax[subcnt].set_xlabel(r'$R_{g}$', weight='bold', size=30)
        ax[subcnt].set_ylabel(r'm', weight='bold', size=30)
            
        if subcnt == i0tot:
            plt.figtext((subp.xbeg+ax_len+ax_sep)/3, subp.ybeg+ax_len+0.2, \
            'Cluster Mass - Radius of Gyration, ' + labels[constants[0]] + ' = ' + "{0:.2f}".format(i0prev), size=35, weight='bold')
            
        # Special handling of the last plot to save
        if i[0] == i0end and subcnt == i0tot:
            
            ax[subcnt].legend(bbox_to_anchor=(1.05,0.,0.53,1.), loc=2, borderaxespad=0., prop={'size': 25}, mode="expand")
            
            # Save the previous figure to this file
            path1 = args.savefile + '/plots'
            path2 = args.savefile + '/plots/CLUSTER'
            if os.path.exists(path1) == False:  
                os.mkdir(path1)
            if os.path.exists(path2) == False: 
                os.mkdir(path2)
                
            plt.savefig(path2 + '/size_rgy_set_' + \
                labels[constants[0]] + '_' + "{0:.2f}".format(i0prev) + '.png', dpi=200, bbox_inches='tight', pad_inches=0.08)
                
            plt.clf()
            



if __name__ == '__main__':
    main()



