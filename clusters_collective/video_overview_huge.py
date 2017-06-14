
##############################################################################

## 
## Generate a huge overview video of filaments with color code as the bead orientations
##

##############################################################################

## load needed libraries and necessary files

import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import h5py
import os
from matplotlib.patches import Rectangle

##############################################################################

def loadSimData(datafile):
    """ load initial simulation data"""
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
        dtSamp, T, box_area, nt, body_length, pe, xil, flexure, tau_D, tau_A, \
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
    pe = Fmc*body_length**2/kT
    xil = Kbend/(kT*body_length)
    flexure = pe/xil
    T = totalStep - ti
    nt = T/nsamp
    nsamp = int(nsamp)
    tau_D = body_length**2*(N+1)/4/kT
    if Fmc != 0:
        tau_A = (N+1)/Fmc
    else:
        tau_A = 0.000001
    nmol = int(M)
    nfil = int(N)
    nbeads = int(L)

    print '\n\n*** SIMULATION PARAMETERS FOR ', datafile
    print 'dt = ', dt
    print 'L = ', body_length
    print 'Pe = ', pe
    print 'xi_p/L = ', xil
    print 'T = ', T
    print 'nfil = ', nfil
    print 'nmol = ', nmol
    print 'nbeads = ', nbeads
    print 'lx = ly = ', Lx
    print "Diffusive time scale is ", tau_D
    print "Advective time scale is ", tau_A
    print '***\n\n'
    
    return


##############################################################################

def gen_mol_info(natoms, nfil):
    """ generate information about molecular ID"""
    
    mol = np.zeros((natoms), dtype = np.int32)
    nmol = natoms/nfil
    k = 0
    for i in range(nmol):
        for j in range(nfil):
            mol[k] = i
            k = k + 1
            
    return mol, nmol
    
##############################################################################    

class Particles:
    """ data structure for storing particle information"""
    
    def __init__(self, path):
        f = h5py.File(path, 'r')
        beads = f['beads']
        self.xi = np.asarray(beads['x'], dtype=float)/B                 # Image particle positions in x
        self.yi = np.asarray(beads['y'], dtype=float)/B                 # Image particle positions in y 
        self.phi = np.asarray(beads['phi'], dtype=float)                # Bead orientation 
        self.cidx = np.asarray(beads['idx'], dtype=int)                 # Cluster index
        f.close()
        
        return
        
    def assign_index(self):
        nbeads = len(self.xi)
        self.idx = np.zeros((nbeads))
        nfil = 51
        for j in range(nbeads):
            self.idx[j] = int(j/nfil)
            
        return
        
##############################################################################
        
class ClusterSize:
    """ data structure for storing cluster sizes"""
    
    def __init__(self, path):
        f = h5py.File(path, 'r')
        size_grp = f['size']
        self.cs = np.asarray(size_grp['size'], dtype=int)       # Cluster sizes
        self.totcs = len(self.cs)                               # Total number of clusters
        f.close()
        
        return

##############################################################################

class Subplots:
    """ plot structure"""
    
    totcnt = -1             # Total number of subplots 
    
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction
        
        return
        
    def addSubplot(self):
        """ add a subplot in the grid structure"""
        
        ## increase the number of subplots in the figure
        
        self.totcnt += 1
        
        ## get indices of the subplot in the figure
        
        self.nx = self.totcnt%(self.tot)
        self.ny = self.totcnt/(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])
        
##############################################################################

def plot_frames(ti, tf, database, savebase, add_patch, change_limits, save, color_switch):
    """ plot some of the time frames"""

    ## set plot properties

    ax_len = 1.0                          # Length of one subplot square box
    ax_b = 0.0                            # Beginning/offset of the subplot in the box
    ax_sep = 0.0                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots    
    fig = plt.figure()
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
    
    ## set more plot properties
    
    if change_limits:
        downlim = 0
        uplim = 500
    else:
        downlim = 0
        uplim = max(Lx,Ly)
    quant_steps = 2056
    if color_switch == 'orient':
        norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)  
    else:
        norm = mpl.colors.Normalize(vmin=0, vmax=2000)
    num_ticks = 5
    
    ## plot the frames
    
    for frame in np.arange(ti, tf, nsamp):
        
        print frame, " of ", tf, " with ", nsamp, " intervals" 
        time = frame*dt
        
        ## data
        
        path = database + "cluster_" + str(frame) + ".hdf5"
        savepath1 = savebase + "frame-" + "{0:05d}".format(int(frame)) + ".png"
        savepath2 = savebase + "frame-" + "{0:05d}".format(int(frame)) + ".eps"
        if add_patch:
            savepath1 = savebase + "frame-" + "{0:05d}".format(int(frame)) + "_z.png"
            savepath2 = savebase + "frame-" + "{0:05d}".format(int(frame)) + "_z.eps"
        
        p = Particles(path)
        p.assign_index()
        
        ## plot
        
        subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
        ax0 = subp.addSubplot()
        
        if color_switch == 'orient':
            line0 = ax0.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                        edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, norm=norm, rasterized=True)
        else:
            line0 = ax0.scatter(p.xi, p.yi, s=1, c=p.idx, cmap=plt.cm.get_cmap('jet',quant_steps), 
                        edgecolors='None', alpha=0.7, vmin=0, vmax=2000, norm=norm, rasterized=True)    
            
        ax0.axis('scaled')
        
        ## title
        
        ax0.set_title("$t/\\tau_{D}$ = " + "{0:.2f}".format(time/tau_D) + \
            ", $t/\\tau_{A}$ = " + "{0:.2f}".format(time/tau_A), fontsize=30)
        
        ## labels
            
        ax0.set_xlabel("$x/r_{0}$", fontsize=30)
        ax0.set_ylabel("$y/r_{0}$", fontsize=30)

        ## limits

        ax0.set_xlim((downlim, uplim))
        ax0.set_ylim((downlim, uplim))
        
        ## ticks
        
        #ax0.xaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
        #ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
        ax0.xaxis.set_ticks([0, 500, 1000, 1500])
        ax0.yaxis.set_ticks([0, 500, 1000, 1500])
        ax0.tick_params(axis='both', which='major', labelsize=30)
        
#        if add_patch:
#            ax0.add_patch( Rectangle( (0,0),500,500,edgecolor='gray',facecolor='gray',alpha=0.3 ) )
            
        
        ## colorbar
        
#        if color_switch == 'orient':
#        
#            cax0 = plt.axes([subp.xbeg+ax_len+0.01, subp.ybeg+ax_len/3, ax_len/4.6, ax_len/4.6], projection='polar')
#            xval = np.arange(-np.pi, np.pi, 0.01)
#            yval = np.ones_like(xval)
#            cax0.scatter(xval, yval, c=xval, s=300, cmap=plt.cm.get_cmap('hsv',quant_steps), norm=norm, linewidths=0)
#            cax0.set_xticks([])
#            cax0.set_yticks([])
#            cax0.set_title('$\\phi$',fontsize=20)
#            cax0.set_rlim([-1,1])
#            cax0.set_axis_off()
                    
        
        ## save
        
        plt.savefig(savepath1, dpi=200, bbox_inches='tight', pad_inches=0.08)
        if save:
            plt.savefig(savepath2, dpi=200, bbox_inches='tight', pad_inches=0.08)        
        fig.clf()

    
    return

##############################################################################        

def main():
    
    ## do argument parsing
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--density", help="Density of the system")
    parser.add_argument("-k", "--kappa", help="Bending rigidity")
    parser.add_argument("-f", "--fp", help="Propulsion force")
    parser.add_argument("-l", "--length", help="Filament length, like many_polymers_5/ or long_filaments/")
    parser.add_argument("-z", "--zoom", help="Decide whether to zoom in or not", action="store_true")
    parser.add_argument("-ti","--init_time", nargs="?", const=10000000, type=int, help="First time frame of the video, in timesteps, you can also leave it empty")
    parser.add_argument("-tf","--fin_time", nargs="?", const=100000000, type=int, help="Last time frame of the video, in timesteps, you can also leave it empty")
    parser.add_argument("-c","--colorbar", help="Decide on the colorbar, id or orient")    
    parser.add_argument("-s","--save", action="store_true", help="Decide whether to save in eps or not")        
    args = parser.parse_args()
    
    ## time frames
    
    ti = args.init_time
    tf = args.fin_time
    
    ## load saved preliminary simulation data into relevant variables
    
    folderbase = "/local/duman/SIMULATIONS/"     # the data is in iff416 (assumption!)
    datafolder = folderbase + args.length + "density_" + args.density + \
        "/kappa_" + args.kappa + "/fp_" + args.fp + "/"
    initfilepath = datafolder + "init_info.txt"
    assert os.path.exists(initfilepath), "init_info.txt does NOT exist for the following file " + initfilepath
    loadSimData(initfilepath)
        
    ## set data folder paths 
        
    clusterpath = datafolder + "CLUSTER/"
    assert os.path.exists(clusterpath), "Cluster analysis has NOT been conducted for : " + clusterpath

    ## set save folder paths 

    savebase = "/usr/users/iff_th2/duman/RolfData/MOVIES/" + args.length
    os.system("mkdir -p " + savebase)
    savepath = savebase + "density_" + args.density + "_kappa_" + args.kappa + \
        "_fp_" + args.fp 
    os.system("mkdir -p " + savepath)
    savepath = savepath + '/'
    
    ## plot certain time frames to generate a movie later on
    
    print "Generating images for the following file : " + datafolder
    print "From frame : " + str(ti) + " to frame : " + str(tf) + " with nsamp : " + str(nsamp)  
    if args.zoom:
        plot_frames(ti, ti+nsamp, clusterpath, savepath, args.zoom, False, args.save, args.colorbar)
        plot_frames(ti+nsamp, tf, clusterpath, savepath, False, True, args.save, args.colorbar)        
    else:              
        plot_frames(ti, tf, clusterpath, savepath, args.zoom, False, args.save, args.colorbar)
    
    return
        
##############################################################################
        
if __name__ == "__main__":
    main()
    
##############################################################################    
    

