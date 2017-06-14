
##############################################################################

## 
## Generate a video with color code of spiral number
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
import math

##############################################################################

def loadSimData(datafile):
    """ load initial simulation data"""
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
    dtSamp, T, box_area, nt, body_length, Pe, persistence, flexure, tau_D, tau_A 

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
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    flexure = Pe/persistence
    T = totalStep - ti
    nt = T/nsamp
    nsamp = int(nsamp)
    tau_D = body_length**2*(N+1)/4/kT
    tau_A = (N+1)/Fmc
    print "Diffusive time scale is ", tau_D
    print "Advective time scale is ", tau_A
    
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
        self.xi = np.asarray(beads['x'], dtype=float)                   # Image particle positions in x
        self.yi = np.asarray(beads['y'], dtype=float)                   # Image particle positions in y 
        self.phi = np.asarray(beads['phi'], dtype=float)                # Bead orientation 
        self.cidx = np.asarray(beads['idx'], dtype=int)                 # Cluster index
        f.close()
        
        return
        
    ###
        
    def assign_index(self):
        nbeads = len(self.xi)
        self.idx = np.zeros((nbeads))
        nfil = 201
        for j in range(nbeads):
            self.idx[j] = int(j/nfil)
            
        return
        
    ###
        
    def neigh_min(self, dx, lx):
        
        dx1 = dx + lx
        dx2 = dx - lx
        if dx**2 < dx1**2 and dx**2 < dx2**2:
            return dx
        if dx1**2 < dx2**2:
            return dx1
        return dx2
        
    ###
        
    def compute_spiral_number(self, xt, yt, lx, ly):
        
        x = np.copy(xt)
        y = np.copy(yt)
        nbeads = len(x)
                    
        # correct pbcs
        for i in range(1,nbeads):
            dx = x[i] - x[i-1]
            dy = y[i] - y[i-1]
            x[i] = x[i-1] + self.neigh_min(dx,lx)
            y[i] = y[i-1] + self.neigh_min(dy,ly)
            
        # compute all bond orientations and center of mass
        phi = np.zeros((nbeads - 1))
        for i in range(1,nbeads):
            dx = x[i] - x[i-1]
            dy = y[i] - y[i-1]
            dphi = math.atan2(dy,dx)
            phi[i-1] = dphi
                            
        # correct for 2*pi periodicity
        phi2 = np.copy(phi)
        nbonds = len(phi)
        for i in range(1,nbonds):
            dphi = phi[i] - phi[i-1]
            if dphi < -np.pi:
                dphi += 2*np.pi
            elif dphi > np.pi:
                dphi -= 2*np.pi
            phi2[i] = phi2[i-1] + dphi
            
        # compute the spiral number
        s = (phi2[-1] - phi2[0])/2/np.pi
        
        return s    
        
    ###
                
    def assign_spiral_number(self, nbeads, nmol, nfil, lx, ly):
        
        self.spiral_number = np.zeros((nbeads), dtype=np.float64)
        self.color_num = np.zeros((nbeads), dtype=np.int32)
        self.color_idx = []
        
        for i in range(nmol):
        
            s = self.compute_spiral_number(self.xi[i*nfil:(i+1)*nfil-1], \
                self.yi[i*nfil:(i+1)*nfil-1], lx, ly)     
            s = math.fabs(s)
            
            if s < 1.6:
                c = 1
                ci = 'grey'
            elif s >= 1.6 and s < 2.2:
                c = 2
                ci = 'k'
            elif s >= 2.2 and s < 3.5:
                c = 3
                ci = 'r'
            elif s >= 3.5 and s < 6.0:
                c = 4
                ci = 'b'
            elif s >= 6.0:
                c = 5
                ci = 'g'
                
            for j in range(nfil):
                k = i*nfil + j
                self.spiral_number[k] = s
                self.color_num[k] = c
                self.color_idx.append(ci)
                
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

def plot_frames(ti, tf, database, savebase, save):
    """ plot some of the time frames"""

    ## set plot properties

    ax_len = 0.9                          # Length of one subplot square box
    ax_b = 0.1                            # Beginning/offset of the subplot in the box
    ax_sep = 0.1                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots    
    fig = plt.figure()
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
    
    ## set more plot properties

    downlim = -5
    uplim = max(Lx,Ly)+5
    quant_steps = 2056
    norm = mpl.colors.Normalize(vmin=0, vmax=5)  
    num_ticks = 5
    
    ## plot the frames
    
    for frame in np.arange(ti, tf, nsamp):
        
        print frame, " of ", tf, " with ", nsamp, " intervals" 
        time = frame*dt
        
        ## data
        
        path = database + "cluster_" + str(frame) + ".hdf5"
        savepath1 = savebase + "frame-" + "{0:05d}".format(int(frame)) + ".png"
        savepath2 = savebase + "frame-" + "{0:05d}".format(int(frame)) + ".eps"
        
        p = Particles(path)
        p.assign_spiral_number(402000, 2000, 201, 1000.0, 1000.0)
        
        ## plot
        
        subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
        ax0 = subp.addSubplot()
        
        line0 = ax0.scatter(p.xi/B, p.yi/B, s=1, c=p.color_idx, \
            edgecolors='None', alpha=0.7, rasterized=True)

        ## title
        
        ax0.set_title("$t/\\tau_{D}$ = " + "{0:.4f}".format(time/tau_D) + \
            ", $t/\\tau_{A}$ = " + "{0:.4f}".format(time/tau_A), fontsize=20)
        
        ## labels
            
        ax0.set_xlabel("$x/r_{0}$", fontsize=30)
        ax0.set_ylabel("$y/r_{0}$", fontsize=30)

        ## limits

        ax0.set_xlim((downlim, uplim))
        ax0.set_ylim((downlim, uplim))
        
        ## ticks
        
        ax0.xaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
        ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
        ax0.tick_params(axis='both', which='major', labelsize=20)

                    
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
    parser.add_argument("-ti","--init_time", nargs="?", const=10000000, type=int, help="First time frame of the video, in timesteps, you can also leave it empty")
    parser.add_argument("-tf","--fin_time", nargs="?", const=100000000, type=int, help="Last time frame of the video, in timesteps, you can also leave it empty")
    parser.add_argument("-s","--save", action="store_true", help="Decide whether to save in eps or not")        
    args = parser.parse_args()
    
    ## time frames
    ti = args.init_time
    tf = args.fin_time
    
    ## load saved preliminary simulation data into relevant variables
    
    folderbase = "/local/duman/SIMULATIONS/"     # the data is in iff416 (assumption!)
    datafolder = folderbase + "long_filaments/density_0.2/kappa_5.0/fp_1.0/"
    initfilepath = datafolder + "init_info.txt"
    assert os.path.exists(initfilepath), "init_info.txt does NOT exist for the following file " + initfilepath
    loadSimData(initfilepath)
        
    ## set data folder paths 
        
    clusterpath = datafolder + "CLUSTER/"
    assert os.path.exists(clusterpath), "Cluster analysis has NOT been conducted for : " + clusterpath

    ## set save folder paths 

    savebase = "/usr/users/iff_th2/duman/RolfData/MOVIES/long_filaments/"
    os.system("mkdir -p " + savebase)
    savepath = savebase + "/density_0.2_kappa_5.0_fp_1.0"
    os.system("mkdir -p " + savepath)
    savepath = savepath + '/'
    
    ## plot certain time frames to generate a movie later on
    
    print "Generating images for the following file : " + datafolder
    print "From frame : " + str(ti) + " to frame : " + str(tf) + " with nsamp : " + str(nsamp)  
    plot_frames(ti, tf, clusterpath, savepath, args.save)
    
    return
        
##############################################################################
        
if __name__ == "__main__":
    main()
    
##############################################################################    
    

