
# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import h5py
import bottleneck as bot
from matplotlib.patches import Rectangle

################################################################################

#
# Function definitions
#
 
################################################################################
 
def loadSimData(datafile):
    """ Load initial simulation data"""
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
    dtSamp, T, box_area, nt, body_length, Pe, persistence, flexure 

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
    
    L = int(L)
    N = int(N)
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
    flexure = Pe/persistence
    
    return

################################################################################

def findIndex(cluster_size, cluster_idx, n):
    """ Find the indices of the largest n clusters in the system"""
    
    ## Cluster identity and index
    cl_id = np.zeros(n, dtype=int)
    cl_idx = np.zeros(n, dtype=int)
    
    ## Calculate the largest n clusters
    size_arr = np.asarray(cluster_size)
    maxs = -bot.partsort(-size_arr, n)[:n]         # n largest clusters
    maxs = sorted(maxs, reverse=True)
    
    ## Index the largest n clusters
    cnt = 0
    for i, cs in enumerate(cluster_size):
        if cs in maxs:
            cl_idx[cnt] = i
            cl_id[cnt] = cluster_idx[i]
            cnt += 1
            if cnt == n:
                break
            else:
                continue
            
    
#    ## Index the largest n clusters
#    for i, cs in enumerate(cluster_size):
#        if cs == maxs[0]:
#            cl_idx[0] = i             
#            cl_id[0] = cluster_idx[i]
#            continue
#        elif cs == maxs[1]:
#            cl_idx[1] = i
#            cl_id[1] = cluster_idx[i]
#            continue
#        elif cs == maxs[2]:
#            cl_idx[2] = i
#            cl_id[2] = cluster_idx[i]
#            continue
 
    return cl_id, cl_idx
            
################################################################################

class Particles:
    """ Data structure for particle data"""
    
    def __init__(self, path):
        f = h5py.File(path, 'r')
        beads = f['beads']
        self.xi = np.asarray(beads['x'], dtype=float)/B                   # Particle positions in x
        self.yi = np.asarray(beads['y'], dtype=float)/B                   # Particle positions in y 
        self.phi = np.asarray(beads['phi'], dtype=float)                  # Bead orientation 
        self.cidx = np.asarray(beads['idx'], dtype=int)                   # Cluster index
        f.close()
        
################################################################################         
        
class ClusterInfo:
    """ Data structure for cluster information"""
    
    def __init__(self, path):
        f = h5py.File(path, 'r')
        props = f['props']    
        self.idx = np.asarray(props['idx'], dtype=int)
        self.size = np.asarray(props['size'], dtype=int)
        self.comx = np.asarray(props['comx'], dtype=float)
        self.comy = np.asarray(props['comy'], dtype=float)
        self.rgysq = np.asarray(props['rgysq'], dtype=float)
        self.pl = np.asarray(props['planarity'], dtype=float)
        self.st = np.asarray(props['straightness'], dtype=float)
        self.sw = np.asarray(props['swirliness'], dtype=float)
        self.ens = np.asarray(props['enstrophy'], dtype=float) 
        fil_grp = props['filament_list']
        self.fils = []
        for el in np.arange(len(self.size)):
            fil_list = np.asarray(fil_grp[str(el)], dtype=int)
            self.fils.append(fil_list)
       
        f.close()
         
        
################################################################################      

class Subplots:
    """ Arrange subplot grid structure (square box is assumed)"""
    
    totcnt = -1             # Total number of subplots 
    
    ## constructor
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction
        
    ## add a subplot in the grid structure
    def addSubplot(self):
        
        # Increase the number of subplots in the figure
        self.totcnt += 1
        
        # Indices of the subplot in the figure
        self.nx = self.totcnt%(self.tot)
        self.ny = self.totcnt/(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])
        
##################################################################

def plot_data(dfolder_1, dfolder_2, sfolder):
    """ plot data"""   
    
    fig = plt.figure()
    downlim = 0
    uplim = Lx
    ax_len = 0.4                          # Length of one subplot square box
    ax_b = 0.1                            # Beginning/offset of the subplot in the box
    ax_sep = 0.1                         # Separation length between two subplots
    total_subplots_in_x = 2               # Total number of subplots in the x direction  
    tick_num = 4                          # Number of ticks
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x)  
    
    comx_curr_1 = []
    comy_curr_1 = []
    comx_trace_1 = []
    comy_trace_1 = []
    comx_curr_2 = []
    comy_curr_2 = []
    comx_trace_2 = []
    comy_trace_2 = []
    
    first_time = True
    cl_id_1 = []
    cl_idx_1 = []
    cl_id_2 = []
    cl_idx_2 = []
    
    n = 2           ## Number of largest clusters to track
    
    ## Time frame loop
    for frame in np.arange(int(ti),int(totalStep),int(nsamp)):
              
        print frame, ' of ', totalStep
        
        ## Load particle data  
        path_1 = dfolder_1 + '/cluster_' + str(frame) + '.hdf5'
        p_1 = Particles(path_1)                    
        path_2 = dfolder_2 + '/cluster_' + str(frame) + '.hdf5'
        p_2 = Particles(path_2)
        
        ## Load cluster information
        cl_1 = ClusterInfo(path_1) 
        cl_2 = ClusterInfo(path_2) 
        
        ## Find the index of the largest n clusters to track them for the first time
        if first_time:
            
            first_time = False
            cl_id_1, cl_idx_1 = findIndex(cl_1.size, cl_1.idx, n)           
            cl_id_2, cl_idx_2 = findIndex(cl_2.size, cl_2.idx, n)           
        
        ## Track the largest n clusters at the beginning
        comx_curr_1 = []
        comy_curr_1 = []
        tracked_idx_1 = []
        comx_curr_2 = []
        comy_curr_2 = []
        tracked_idx_2 = []        
        for j in np.arange(len(cl_1.comx)):
            if cl_1.idx[j] in cl_id_1:
                comx_curr_1.append(cl_1.comx[j])
                comy_curr_1.append(cl_1.comy[j])
                comx_trace_1.append(cl_1.comx[j])
                comy_trace_1.append(cl_1.comy[j])   
                tracked_idx_1.append(cl_1.fils[j])   
            if cl_2.idx[j] in cl_id_2:
                comx_curr_2.append(cl_2.comx[j])
                comy_curr_2.append(cl_2.comy[j])
                comx_trace_2.append(cl_2.comx[j])
                comy_trace_2.append(cl_2.comy[j])   
                tracked_idx_2.append(cl_2.fils[j])               
        
#        ## Save clusters larger than 50 filaments
#        comx_curr = []
#        comy_curr = []
#        tracked_idx = []
#        for j in np.arange(len(cl.comx)):
#            if cl.size[j] > 50:
#                comx_curr.append(cl.comx[j])
#                comy_curr.append(cl.comy[j])
#                comx_trace.append(cl.comx[j])
#                comy_trace.append(cl.comy[j])   
#                tracked_idx.append(cl.fils[j])
         
#        ## Get the indices of the overlaid trajectories of beads
#        bead_pos_x = []
#        bead_pos_y = []
#        for j, fil_list in enumerate(tracked_idx):
#            for fil_idx in fil_list:
#                for el in np.arange(N):
#                    bead_pos_x.append(p.xi[int(fil_idx*N+el)])
#                    bead_pos_y.append(p.yi[int(fil_idx*N+el)])
        
        ## Highlight the beads of tracked clusters
        clidx_1 = []
        clidx_2 = []
        for i in range(L):
            clidx_1.append('gray')
            clidx_2.append('gray')
        
        ## Get the indices of the overlaid trajectories of beads
        for j, fil_list in enumerate(tracked_idx_1):
            for fil_idx in fil_list:
                for el in np.arange(N):
                    clidx_1[int(fil_idx*N+el)] = 'blue'
        for j, fil_list in enumerate(tracked_idx_2):
            for fil_idx in fil_list:
                for el in np.arange(N):
                    clidx_2[int(fil_idx*N+el)] = 'blue'

        ## Plot the data
        ax1 = subp.addSubplot()
        
        #size = np.pi * ( (0.001*cl.size)**2 )
        #ax1.scatter(cl.comx, cl.comy, s=size, edgecolors='None', alpha=1)
        ax1.scatter(comx_curr_2, comy_curr_2, s=30, color='red', alpha=1)
        ax1.scatter(comx_trace_2, comy_trace_2, s=5, color='red', edgecolors='None', alpha=0.7)
        ax1.scatter(p_2.xi, p_2.yi, s=0.1, c=clidx_2, edgecolors='None', alpha=0.4)
        #ax1.scatter(bead_pos_x, bead_pos_y, s=0.5, color='gray', edgecolors='None', alpha=0.1)
        ax1.set_xlabel('x/r',fontsize=20)
        ax1.set_ylabel('y/r',fontsize=20)
        ax1.set_xlim((downlim,uplim))
        ax1.set_ylim((downlim,uplim))
        ax1.xaxis.set_ticks( np.linspace(downlim, uplim, num=tick_num, endpoint=True) )
        ax1.yaxis.set_ticks( np.linspace(downlim, uplim, num=tick_num, endpoint=True) )
        ax1.tick_params(axis='both', which='major', labelsize=12)         
               
        ax2 = subp.addSubplot()
        
        #size = np.pi * ( (0.001*cl.size)**2 )
        #ax2.scatter(cl.comx, cl.comy, s=size, edgecolors='None', alpha=1)
        ax2.scatter(comx_curr_1, comy_curr_1, s=30, color='red', alpha=1)
        ax2.scatter(comx_trace_1, comy_trace_1, s=5, color='red', edgecolors='None', alpha=0.7)
        ax2.scatter(p_1.xi, p_1.yi, s=0.1, c=clidx_1, edgecolors='None', alpha=0.4)
        #ax2.scatter(bead_pos_x, bead_pos_y, s=0.5, color='gray', edgecolors='None', alpha=0.1)
        ax2.set_xlabel('x/r',fontsize=20)
        #ax2.set_ylabel('y/r',fontsize=20)
        ax2.set_xlim((downlim,uplim))
        ax2.set_ylim((downlim,uplim))
        ax2.xaxis.set_ticks( np.linspace(downlim, uplim, num=tick_num, endpoint=True) )
        ax2.yaxis.set_ticks( np.linspace(downlim, uplim, num=tick_num, endpoint=True) )
        ax2.tick_params(axis='both', which='major', labelsize=12)         
        plt.setp(ax2.get_yticklabels(),visible=False)
        
        if frame == 16350000:
            ax1.set_rasterized(True)
            ax2.set_rasterized(True)
            plt.savefig(sfolder + '/frame-'+'{0:05d}'.format(int(frame))+'.eps',dpi=200,bbox_inches='tight',pad_inches=0.08)
            plt.savefig(sfolder + '/frame-'+'{0:05d}'.format(int(frame))+'.pdf',dpi=200,bbox_inches='tight',pad_inches=0.08)            
        plt.savefig(sfolder + '/frame-'+'{0:05d}'.format(int(frame))+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)

        plt.clf()
        
    
    return

##################################################################

def main():
    """ main function, called when the script is started"""
    
    ## Argument parsing (command line options)
    parser = argparse.ArgumentParser()
    parser.add_argument("datafolder_1", help="File containing simulation data")    
    parser.add_argument("datafolder_2", help="File containing simulation data")        
    parser.add_argument("savefolder", help="File in which figures should be saved")
    args = parser.parse_args()
    
    ## Load saved preliminary simulation data into relevant variables
    #loadSimData("/local/duman/SIMULATIONS/many_polymers_5/density_0.08/kappa_125.0/fp_0.8/init_info.txt")
    loadSimData(args.datafolder_1+"/init_info.txt")
    datafolder_1 = args.datafolder_1+"/CLUSTER"
    datafolder_2 = args.datafolder_2+"/CLUSTER"
    
    savefolder = args.savefolder+"/traj"
    os.system("mkdir -p " + savefolder)
    
    plot_data(datafolder_1, datafolder_2, savefolder)


        
if __name__ == '__main__':
    main()        