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
 
# Load initial simulation data      
def loadSimData(datafile):
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
    dtSamp, T, box_area, nt, body_length, Pe, persistence 

    datafile = open(datafile,"r")
    for line in datafile:
        A = line.split()
        if A[0] == "dt":                    # Time interval between MD steps
            dt = float(A[1])
        elif A[0] == "ti":                  # Beginning time for data acquisition
            ti = float(A[1])
        elif A[0] == "Lx":                  # Box size in x
            Lx = float(A[1])            
        elif A[0] == "Ly":                  # Box size in y
            Ly = float(A[1])
        elif A[0] == "totalStep":           # Total MD steps
            totalStep = float(A[1])
        elif A[0] == "nsamp":               # Data sampling frequency
            nsamp = float(A[1])
        elif A[0] == "nfil":                # Number of particles per polymer
            N = float(A[1])
        elif A[0] == "L":                   # Number of particles
            L = float(A[1])
        elif A[0] == "B":                   # Bond length between particles of a body
            B = float(A[1])
        elif A[0] == "kT":                  # Boltzmann constant*Temperature
            kT = float(A[1])
        elif A[0] == "Fmc":                 # Self propulsion force constant
            Fmc = float(A[1])     
        elif A[0] == "Kbend":               # Bending constant
            Kbend = float(A[1])
    
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
    parser.add_argument("initfile", help="Folder containing initial simulation data")
    parser.add_argument("datafile", help="Folder containing data")
    parser.add_argument("analysis", help="Type of analysis")
    parser.add_argument("savefile", help="Folder in which figures should be saved")
    args = parser.parse_args()
    
    # Load saved preliminary simulation data into relevant variables
    # This needs to be changed!!!
    loadSimData(args.initfile)
    
    # Plot properties
    ax_len = 0.4                      # Length of one subplot square box
    ax_b = 0.1                        # Beginning/offset of the subplot in the box
    ax_sep = ax_len/1.3                 # Separation length between two subplots
    total_subplots_in_y = 2           # Total number of subplots
    
    # Index the data
    density = [2.0, 5.0, 10.0, 20.0]
    kappa = [2.5, 12.5, 62.5, 200.0]
    fp = [0.0, 0.08, 0.024, 0.24, 0.8, 2.4, 8.0]
    
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
        
        tmp.update({'folder':folder+args.analysis+'.data', 'Pe':peclet, 'xi_L':persistence_over_L})
        files.append(tmp)
            

    Data = pd.DataFrame(files,columns=files[0].keys())
    
    # Choose the comparison parameter for plotting 
    comparison_array = ['density', 'Pe', 'xi_L']
    compareby = 'xi_L'         # make this a command line option
    constants = []
    for element in comparison_array:
        if element != compareby:
            constants.append(element)

    # Group the data accordingly with the chosen comparison parameter
    group = Data.groupby(constants)
    
    # Labels
    labels={'Pe':r"Pe = ",
       'xi_L':r"$\xi_p/L = $",
       'density':'Density'}
       

    # Plot the data
    i0prev = density[0]
    i1prev = -1.0
    i0end = density[-1]
    i0cnt = 0
    i0tot = 6

    ax = []
    fig = plt.figure() 
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_y)
    subcnt = -1
    
    for i, grp in group:
        # i are the key values,          [tuple]
        # grp is the grouped element     [DataFrame]
        
        #print i, i[0], i[1], i0prev, i1prev
            
            
        # First constant will be the filename (density in this case)    
        if i0prev != i[0]:
                        
            # Save the previous figure to this file
            path1 = args.savefile + '/plots'
            path2 = args.savefile + '/plots/RDF_SITE'
            if os.path.exists(path1) == False:  
                os.mkdir(path1)
            if os.path.exists(path2) == False:
                os.mkdir(path2)
                
            plt.savefig(path2 + '/site_set_' + \
                labels[constants[0]] + '_' + str(i0prev) + '.png', dpi=200, bbox_inches='tight', pad_inches=0.08)
                
            plt.clf() 
            
            i0prev = i[0]
            
            # Create a new figure
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
        minxprev = 0
        maxxprev = 0
        for j, f in enumerate(grp['folder']):
            #print j , i
            #print f
            if os.path.exists(f):
                order_data = pd.read_csv(f, sep='\t', skiprows=1, header=0) 
                x = (order_data['r_max'] + order_data['r_min'])/2
                y = order_data['rdf'] 
            #print subcnt
            ax[subcnt].plot(x,y,label=labels[compareby] + str( grp['xi_L'].iloc[j] ))
            
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
        
        if maxx != 0:
            ax[subcnt].set_xlim((minx,maxx))
            #ax[subcnt].set_ylim((0,1))
            ax[subcnt].xaxis.set_ticks( np.arange(int(minx),int(maxx),int((maxx-minx)/5)) )
            #ax[subcnt].yaxis.set_ticks( np.arange(0.,1.,0.3) )
            ax[subcnt].tick_params(axis='both', which='major', labelsize=12)
        else:
            plt.setp(ax[subcnt].get_xticklabels(),visible=False)
            plt.setp(ax[subcnt].get_yticklabels(),visible=False)
            
            
        # Figure texts
        ax[subcnt].set_title(labels[constants[1]] + str(i[1]), weight='bold', size='large')
        ax[subcnt].legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0., prop={'size': 12})
        ax[subcnt].set_xlabel('$r [\\sigma]$', weight='bold', size='large')
        ax[subcnt].set_ylabel('$g(r)$', weight='bold', size='large')
            
        if subcnt == i0tot:
            plt.figtext((subp.xbeg+ax_len+ax_sep)/2.3, subp.ybeg+2*ax_len+1.3*ax_sep, \
            'Radial Distribution Function', size='xx-large', weight='bold')
            
        if i[0] == i0end and subcnt == i0tot:
            # Save the previous figure to this file
            path1 = args.savefile + '/plots'
            path2 = args.savefile + '/plots/RDF_SITE'
            if os.path.exists(path1) == False:  
                os.mkdir(path1)
            if os.path.exists(path2) == False:
                os.mkdir(path2)
                
            plt.savefig(path2 + '/site_set_' + \
                labels[constants[0]] + '_' + str(i0prev) + '.png', dpi=200, bbox_inches='tight', pad_inches=0.08)
                
            plt.clf()
            
        #print i, i[0], i[1], i0prev, i1prev



if __name__ == '__main__':
    main()



#
#
## Function for reading data from file
#def readData( indata ):
#
#	#d = indata.shape[1] # dimension information
#	N = indata.shape[0] # number of data points
#
#	x = np.zeros(N)
#	y = np.zeros(N)
#	for i in np.arange(N):
#		x[i] = indata[i][0]
#		y[i] = indata[i][1]
#	return x, y, N
#
## Functions for plotting specific analysis modules
#def linlog_plotter( x, y, ax, label_name, color_idx ):
#
#	ax[0].plot(x,y,label=label_name,linewidth=1.0,color=color_idx)
#	ax[1].loglog(x,y,label=label_name,linewidth=1.0,color=color_idx)
#
#	return ax
#
#def semilogx_plotter( x, y, ax, label_name, color_idx ):
#
#	ax.semilogx(x,y,label=label_name,linewidth=1.0,color=color_idx)
#
#	return ax
#
#def lin_plotter( x, y, ax, label_name, color_idx ):
#
#	ax.plot(x,y,label=label_name,linewidth=1.0,color=color_idx)
#
#	return ax
#
#
##
## Read the data
##
#
## Folder names
##phi = [0.3, 0.4, 0.5, 0.6, 0.7]
##eps = [10, 12, 14]
##Fm = [1, 2]
#
#phi = [0.3]
#eps = [10, 12, 14]
#Fm = [1, 2]
#
#
##dir_name = 'eps' + str(eps[0]) + 'Fm' + str(Fm[0])
#dir_name = args.cell_file
#os.mkdir(dir_name)
#
#folders = []
#for P in phi:
#    for E in eps:
#        for F in Fm:
#            folders.append('phi'+str(P)+'eps'+str(E)+'fm'+str(F))
#
#
## Analysis modules
#analysis = ['msd_post','inter_scattering','dynamic_struct','static_struct', \
#'pair_corr','vacf_cm','sp_vacf','bond_corr']
#
#
#for a in analysis:
#    
#    x = []
#    y = []
#    cnt = 0
#    
#    # Load the data from the files
#    for folder in folders:
#
#        path = '../' + folder + '/' + a + '.txt'
#        indata = np.loadtxt(path,dtype=float)
#        outdata = readData( indata )
#        x.append( outdata[0] )
#        y.append( outdata[1] )
#        N = outdata[2]
#        cnt += 1
#        print path, cnt
#        
#    # Set general graph properties
#    plt.figure(figsize=(20,17))
#    colors = plt.cm.spectral(np.linspace(0,1,cnt))
#    if a == 'msd_post':
#        fig, ax = plt.subplots(1,2)
#    else:
#        fig, ax = plt.subplots(1,1)
#    
#    # Plot the data
#    for i in np.arange(cnt):
#        if a == 'msd_post':
#            ax = linlog_plotter(x[i]/tr,y[i],ax,folders[i],colors[i])
#        elif a == 'inter_scattering':
#            ax = semilogx_plotter(x[i]/tr,y[i],ax,folders[i],colors[i])
#        elif a == 'vacf_cm':
#            ax = lin_plotter(x[i]/tr,y[i],ax,folders[i],colors[i])
#        else:
#            ax = lin_plotter(x[i],y[i],ax,folders[i],colors[i])
#            
#    # Customize the plots
#    if a == 'msd_post':
#        ax[1].legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0., prop={'size': 9})
#        ax[0].tick_params(axis='both',which='major',labelsize=11)
#        ax[1].tick_params(axis='both',which='major',labelsize=11)
#    else:
#        ax.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0., prop={'size': 9})
#        ax.tick_params(axis='both',which='major',labelsize=15)
#        
#    if a == 'static_struct':
#        ax.set_title('Static structure factor')
#        ax.set_xlabel('k*R')
#        ax.set_ylabel('S(k)')
#    elif a == 'inter_scattering':
#        ax.set_title('Intermediate scattering function')
#        ax.set_xlabel('$t/\\tau_{R}$')
#        ax.set_ylabel('F(k,t)') 
#    elif a == 'dynamic_struct':
#        ax.set_title('Dynamic structure factor')
#        ax.set_xlabel('w')
#        ax.set_ylabel('DS(k,w)')
#    elif a == 'msd_post':
#        ax[0].set_xlabel('$t/\\tau_{R}$')
#        ax[0].set_ylabel('$\Delta r^2/4R^2$')
#        ax[1].set_xlabel('$t/\\tau_{R}$')
#    elif a == 'pair_corr':
#        ax.set_title('Pair correlation function')
#        ax.set_xlabel('r/2R')
#        ax.set_ylabel('g(r)') 
#        ax.set_xlim((0,9)) 
#    elif a == 'sp_vacf':
#        ax.set_title('Spatial velocity correlation')
#        ax.set_xlabel('r/2R')
#        ax.set_ylabel('$C_{vv}(r)$')
#    elif a == 'vacf_cm':
#        ax.set_title('Velocity autocorrelation')
#        ax.set_xlabel('$t/\\tau_{R}$')
#        ax.set_ylabel('$C_{vv}(t)$') 
#    elif a == 'bond_corr':
#        ax.set_title('Spatial orientational order correlation')
#        ax.set_xlabel('r/2R')
#        ax.set_ylabel('$g_6(r)$') 
#        
#    plt.savefig(a+'.eps', bbox_inches='tight')
#    
#    shu.move(a+'.eps',dir_name)


