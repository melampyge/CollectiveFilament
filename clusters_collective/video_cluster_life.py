#!/usr/bin/python


# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import bottleneck as bot
from matplotlib.patches import Rectangle
import os
import h5py

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
            N = int(float(A[-1]))
        elif A[0] == "L":                   # Number of particles
            L = int(float(A[-1]))
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
    flexure = Pe/persistence
    
################################################################################
    
def colorToBeads(clidx, filament_list, color, part_list_x, part_list_y, part_list_phi):
    """ Assign colors to beads indices based on the cluster they belong"""
    cl_beads_x = []
    cl_beads_y = []
    cl_beads_p = []
    for i in filament_list:
        for j in range(N):
            ind = int(i*N+j)
            clidx[ind] = color
            cl_beads_x.append(part_list_x[ind])
            cl_beads_y.append(part_list_y[ind])
            cl_beads_p.append(part_list_phi[ind])
            
    xmin = min(cl_beads_x)
    ymin = min(cl_beads_y)
    xmax = max(cl_beads_x)
    ymax = max(cl_beads_y)
    
    flag = False
    x_span = xmax-xmin
    y_span = ymax-ymin
    if x_span > Lx/2:
        flag = True
    elif y_span > Ly/2:
        flag = True
    
    return cl_beads_x, cl_beads_y, cl_beads_p, xmin, xmax, ymin, ymax, flag        
    

################################################################################

def findIndex(cluster_size, cluster_idx, tracked_clusters, index_tracked_clusters):
    """ Find the indices of the largest 5 clusters in the system"""
    
    print 'Update is conducting'
    # Calculate the largest 5 clusters
    size_arr = np.asarray(cluster_size)
    maxs = -bot.partsort(-size_arr, 5)[:5]         # 5 largest clusters
    maxs = sorted(maxs, reverse=True)
    
    # Index the largest 5 clusters
    for i, cs in enumerate(cluster_size):
        if cs == maxs[0]:
            index_tracked_clusters[0] = i             
            tracked_clusters[0] = cluster_idx[i]
            continue
        elif cs == maxs[1]:
            index_tracked_clusters[1] = i
            tracked_clusters[1] = cluster_idx[i]
            continue
        elif cs == maxs[2]:
            index_tracked_clusters[2] = i
            tracked_clusters[2] = cluster_idx[i]
            continue
        elif cs == maxs[3]:
            index_tracked_clusters[3] = i
            tracked_clusters[3] = cluster_idx[i]
            continue
        elif cs == maxs[4]:
            index_tracked_clusters[4] = i
            tracked_clusters[4] = cluster_idx[i] 
            continue
            
################################################################################

def findSpecificIndex(cluster_size, cluster_idx, tracked_clusters, index_tracked_clusters, indi):
    """ Find the indices of the largest cluster in the system"""
      
    # Calculate the largest cluster
    maxi = max(cluster_size)
    
    # Index the largest cluster
    for i, cs in enumerate(cluster_size):
        if cs == maxi:
            index_tracked_clusters[indi] = i             
            tracked_clusters[indi] = cluster_idx[i]   
            break
                        
################################################################################

def findArrayIndex(cluster_idx, tracked_clusters, index_tracked_clusters):
    """ Find the array index of cluster indices to access elements of the cluster information"""
    
    found = np.zeros(5, dtype=int)
    
    # Index the clusters again since the unique cluster ids may be changed
    for i, indi in enumerate(cluster_idx):
        if indi == tracked_clusters[0]:
            index_tracked_clusters[0] = i
            found[0] = 1
            continue
        elif indi == tracked_clusters[1]:
            index_tracked_clusters[1] = i
            found[1] = 1
            continue
        elif indi == tracked_clusters[2]:
            index_tracked_clusters[2] = i
            found[2] = 1
            continue
        elif indi == tracked_clusters[3]:
            index_tracked_clusters[3] = i
            found[3] = 1
            continue
        elif indi == tracked_clusters[4]:
            index_tracked_clusters[4] = i
            found[4] = 1
            continue
            
    return found
                
################################################################################

#
# Class definitions
#

################################################################################

class Particles:
    """ Data structure for particle data"""
    
    def __init__(self, path):
        f = h5py.File(path, 'r')
        beads = f['beads']
        self.xi = np.asarray(beads['x'], dtype=float)/B                 # Image particle positions in x
        self.yi = np.asarray(beads['y'], dtype=float)/B                 # Image particle positions in y 
        self.phi = np.asarray(beads['phi'], dtype=float)                # Bead orientation 
        self.cidx = np.asarray(beads['idx'], dtype=int)                 # Cluster index
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
        
################################################################################
      
def main():
      
    # Argument parsing (command line options)
    parser = argparse.ArgumentParser()
    parser.add_argument("initfile", help="File containing initial simulation data")
    parser.add_argument("savefile", help="File in which figures should be saved")
    args = parser.parse_args()
    
    # Load saved preliminary simulation data into relevant variables
    loadSimData(args.initfile+"/init_info.txt")
    datafile = args.initfile+"/CLUSTER"
    if os.path.exists(args.savefile+'/climg') == False:
        os.mkdir(args.savefile+'/climg')
        
    # Plot properties
    downlim = -5
    uplim = max(Lx,Ly)+5
    ax_len = 0.38                         # Length of one subplot square box
    ax_b = 0.1                            # Beginning/offset of the subplot in the box
    ax_sep = 0.2                         # Separation length between two subplots
    total_subplots_in_y = 2               # Total number of subplots
    tick_interval = int(uplim/5)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    
    first_time = True
    update_limits = True
    limit_fac = 150
    tracked_clusters = np.zeros(5, dtype=int)
    index_tracked_clusters = np.zeros(5, dtype=int)
    fig = plt.figure()
    
    ################################################################################
    
    # Time frame loop
    for frame in np.arange(int(ti),int(totalStep),int(nsamp)):
        
        
        print frame 
    
        #
        # Load and set the data
        #  
        
        # Load particle data  -- 1xL --
        path = datafile + '/cluster_' + str(frame) + '.hdf5'
        p = Particles(path)                    
    
        
    #    # Load cluster sizes
    #    path = datafile + '/cluster_sizes_' + str(frame) + '.txt'
    #    clsizes = ClusterSize(path)
        
        
        # Load cluster information
        #path = datafile + '/cluster_' + str(frame) + '.hdf5'
        clinfo = ClusterInfo(path)
    
    
        # Find the index of the largest 5 clusters to track them for the first time
        if first_time:
            
            first_time = False
            findIndex(clinfo.size, clinfo.idx, tracked_clusters, index_tracked_clusters)   
        
        # If it is not the first time, then check if update is necessary to get new tracked clusters    
        else:
            
            # Update the array index of tracked clusters to access the cluster properties 
            # if the tracked cluster cannot be found in the new cluster list, update all the clusters
            found = findArrayIndex(clinfo.idx, tracked_clusters, index_tracked_clusters)  
            
            print found, index_tracked_clusters, tracked_clusters
            
            # If one of the clusters is not found in the new time frame, update all the clusters and initiate the break        
            for i in found:
                if i == 0:
                    findIndex(clinfo.size, clinfo.idx, tracked_clusters, index_tracked_clusters)
                    update_limits = True
                    break
                
            # If the break is not initiated, that is, if the the tracked clusters are here in this time frame as well
            # then this else condition will be initiated
            else:
                
                # Check whether an update to the tracked clusters is necessary or not
                for i in range(len(index_tracked_clusters)):
                        
                    # Update the index for ALL the tracked clusters, if one of them has a size smaller than 15 filaments
                    if clinfo.size[index_tracked_clusters[i]] < 15:
                        findIndex(clinfo.size, clinfo.idx, tracked_clusters, index_tracked_clusters)
                        update_limits = True
                        break    # No need to run after the update is done
                    
            print found, index_tracked_clusters, tracked_clusters
                
         
        # Highlight the indexed clusters
        clidx = []
        for i in range(L):
            clidx.append('gray')
            
        flags = ['False', 'False', 'False', 'False', 'False']
        c0_x, c0_y, c0_p, min0_x, max0_x, min0_y, max0_y, flags[0] = colorToBeads(clidx, clinfo.fils[index_tracked_clusters[0]], 'maroon', p.xi, p.yi, p.phi)
        c1_x, c1_y, c1_p, min1_x, max1_x, min1_y, max1_y, flags[1] = colorToBeads(clidx, clinfo.fils[index_tracked_clusters[1]], 'darkviolet', p.xi, p.yi, p.phi)
        c2_x, c2_y, c2_p, min2_x, max2_x, min2_y, max2_y, flags[2] = colorToBeads(clidx, clinfo.fils[index_tracked_clusters[2]], 'red', p.xi, p.yi, p.phi)
        c3_x, c3_y, c3_p, min3_x, max3_x, min3_y, max3_y, flags[3] = colorToBeads(clidx, clinfo.fils[index_tracked_clusters[3]], 'violet', p.xi, p.yi, p.phi)
        c4_x, c4_y, c4_p, min4_x, max4_x, min4_y, max4_y, flags[4] = colorToBeads(clidx, clinfo.fils[index_tracked_clusters[4]], 'lightsalmon', p.xi, p.yi, p.phi)
    
        #first_time = [cond for cond in flags if cond == True]
        
        if update_limits:
            update_limits = False
            min0_x_beg = min0_x-limit_fac
            max0_x_beg = max0_x+limit_fac
            min0_y_beg = min0_y-limit_fac
            max0_y_beg = max0_y+limit_fac        
            
            min1_x_beg = min1_x-limit_fac
            max1_x_beg = max1_x+limit_fac
            min1_y_beg = min1_y-limit_fac
            max1_y_beg = max1_y+limit_fac   
            
            min2_x_beg = min2_x-limit_fac
            max2_x_beg = max2_x+limit_fac
            min2_y_beg = min2_y-limit_fac
            max2_y_beg = max2_y+limit_fac
            
            min3_x_beg = min3_x-limit_fac
            max3_x_beg = max3_x+limit_fac
            min3_y_beg = min3_y-limit_fac
            max3_y_beg = max3_y+limit_fac
        
            min4_x_beg = min4_x-limit_fac
            max4_x_beg = max4_x+limit_fac
            min4_y_beg = min4_y-limit_fac
            max4_y_beg = max4_y+limit_fac
        
    
    ################################################################################
    
        
        #
        # Plot the data
        #   
       
        # Set up the plot and the subplots structure
        #fig = plt.figure()
        
        subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_y)  
        
        
        # Cluster id 
        ax0 = subp.addSubplot()
        ax0.scatter(p.xi, p.yi, s=1, c=clidx, edgecolors='None', alpha=0.7)
        ax0.set_xlabel('x/b',fontsize=8)
        ax0.set_ylabel('y/b',fontsize=8)
        ax0.set_title('Cluster id')
        ax0.set_xlim((downlim,uplim))
        ax0.set_ylim((downlim,uplim))
        ax0.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
        ax0.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
        ax0.tick_params(axis='both', which='major', labelsize=8)    
     
        
        # Bead orientations
        ax1 = subp.addSubplot()  
        quant_steps = 2056
        norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)
        
        line1 = ax1.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, norm=norm)
        ax1.set_ylabel('y/b',fontsize=8)
        ax1.set_title('Bead Orientations')
        ax1.set_xlim((downlim,uplim))
        ax1.set_ylim((downlim,uplim))
        ax1.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
        ax1.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
        ax1.tick_params(axis='both', which='major', labelsize=8)
        
        cax1 = plt.axes([subp.xbeg+ax_len+0.01, subp.ybeg+ax_len/3, ax_len/4.6, ax_len/4.6], projection='polar')
        xval = np.arange(-np.pi, np.pi, 0.01)
        yval = np.ones_like(xval)
        cax1.scatter(xval, yval, c=xval, s=300, cmap=plt.cm.get_cmap('hsv',quant_steps), norm=norm, linewidths=0)
        cax1.set_yticks([])
        cax1.set_xticks([])
        cax1.set_title('$\\phi$',fontsize=10)
        cax1.set_rlim([-1,1])
        cax1.set_axis_off()
         
        ax1.add_patch( Rectangle( (min0_x,min0_y),max0_x-min0_x,max0_y-min0_y,edgecolor='gray',facecolor='gray',alpha=0.3 ) )
        ax1.add_patch( Rectangle( (min1_x,min1_y),max1_x-min1_x,max1_y-min1_y,edgecolor='gray',facecolor='gray',alpha=0.3 ) )
        ax1.add_patch( Rectangle( (min2_x,min2_y),max2_x-min2_x,max2_y-min2_y,edgecolor='gray',facecolor='gray',alpha=0.3 ) )
        ax1.add_patch( Rectangle( (min3_x,min3_y),max3_x-min3_x,max3_y-min3_y,edgecolor='gray',facecolor='gray',alpha=0.3 ) )
        ax1.add_patch( Rectangle( (min4_x,min4_y),max4_x-min4_x,max4_y-min4_y,edgecolor='gray',facecolor='gray',alpha=0.3 ) )
    
    
        # Cluster 2 / Dark violet
        ax2 = subp.addSubplot()
        line2 = ax2.scatter(c1_x, c1_y, s=1, c=c1_p, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=1.0, vmin=-np.pi, vmax=np.pi, norm=norm)  
        ax2.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=0.2, vmin=-np.pi, vmax=np.pi, norm=norm)                
        ax2.plot(clinfo.comx[index_tracked_clusters[1]],clinfo.comy[index_tracked_clusters[1]],'*',color='darkviolet')                
        ax2.set_title('Dark Violet')
        ax2.set_xlabel('x/b',fontsize=8)        
        ax2.set_ylabel('y/b',fontsize=8)    
        ax2.set_xlim((min1_x_beg,max1_x_beg))
        ax2.set_ylim((min1_y_beg,max1_y_beg))
        plt.setp(ax2.get_xticklabels(),visible=False)
        plt.setp(ax2.get_yticklabels(),visible=False)
        #ax2.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
        #ax2.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
        ax2.tick_params(axis='both', which='major', top='off', bottom='off', left='off', right='off')  
    
        textstr = '$R_g=%.1f$\n$S=%.2f$\n$F_{tr}=%.1f$\n$F_{rot}=%.3f$\n$\Omega=%.4f$'%(clinfo.rgysq[index_tracked_clusters[1]], clinfo.pl[index_tracked_clusters[1]], clinfo.st[index_tracked_clusters[1]], clinfo.sw[index_tracked_clusters[1]], clinfo.ens[index_tracked_clusters[1]])
        ax2.text(1.05, 0.95, textstr, transform=ax2.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
            
      
        # Cluster 1 / Maroon
        ax3 = subp.addSubplot()  
        line3 = ax3.scatter(c0_x, c0_y, s=1, c=c0_p, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=1.0, vmin=-np.pi, vmax=np.pi, norm=norm)  
        ax3.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=0.2, vmin=-np.pi, vmax=np.pi, norm=norm)                               
        ax3.plot(clinfo.comx[index_tracked_clusters[0]],clinfo.comy[index_tracked_clusters[0]],'*',color='maroon')
        ax3.set_title('Maroon')
        ax3.set_xlabel('x/b',fontsize=8)    
        ax3.set_ylabel('y/b',fontsize=8)        
        ax3.set_xlim((min0_x_beg,max0_x_beg))
        ax3.set_ylim((min0_y_beg,max0_y_beg))
        plt.setp(ax3.get_xticklabels(),visible=False)
        plt.setp(ax3.get_yticklabels(),visible=False)
        #ax3.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
        #ax3.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
        ax3.tick_params(axis='both', which='major', top='off', bottom='off', left='off', right='off')
        
        textstr = '$R_g=%.1f$\n$S=%.2f$\n$F_{tr}=%.1f$\n$F_{rot}=%.3f$\n$\Omega=%.4f$'%(clinfo.rgysq[index_tracked_clusters[0]], clinfo.pl[index_tracked_clusters[0]], clinfo.st[index_tracked_clusters[0]], clinfo.sw[index_tracked_clusters[0]], clinfo.ens[index_tracked_clusters[0]])
        ax3.text(1.05, 0.95, textstr, transform=ax3.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
            
            
        # Cluster 4 / Violet
        ax4 = subp.addSubplot()
        line4 = ax4.scatter(c3_x, c3_y, s=1, c=c3_p, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=1.0, vmin=-np.pi, vmax=np.pi, norm=norm)  
        ax4.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=0.2, vmin=-np.pi, vmax=np.pi, norm=norm)                               
        ax4.plot(clinfo.comx[index_tracked_clusters[3]],clinfo.comy[index_tracked_clusters[3]],'*',color='violet')                
        ax4.set_title('Violet')
        ax4.set_xlabel('x/b',fontsize=8)    
        ax4.set_ylabel('y/b',fontsize=8)        
        ax4.set_xlim((min3_x_beg,max3_x_beg))
        ax4.set_ylim((min3_y_beg,max3_y_beg))
        plt.setp(ax4.get_xticklabels(),visible=False)
        plt.setp(ax4.get_yticklabels(),visible=False)
        #ax4.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
        #ax4.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
        ax4.tick_params(axis='both', which='major', top='off', bottom='off', left='off', right='off')  
        
        textstr = '$R_g=%.1f$\n$S=%.2f$\n$F_{tr}=%.1f$\n$F_{rot}=%.3f$\n$\Omega=%.4f$'%(clinfo.rgysq[index_tracked_clusters[3]], clinfo.pl[index_tracked_clusters[3]], clinfo.st[index_tracked_clusters[3]], clinfo.sw[index_tracked_clusters[3]], clinfo.ens[index_tracked_clusters[3]])
        ax4.text(1.05, 0.95, textstr, transform=ax4.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
            
            
        # Cluster 3 / Red
        ax5 = subp.addSubplot()
        line5 = ax5.scatter(c2_x, c2_y, s=1, c=c2_p, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=1.0, vmin=-np.pi, vmax=np.pi, norm=norm) 
        ax5.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=0.2, vmin=-np.pi, vmax=np.pi, norm=norm)                               
        ax5.plot(clinfo.comx[index_tracked_clusters[2]],clinfo.comy[index_tracked_clusters[2]],'*',color='red')                
        ax5.set_title('Red')
        ax5.set_xlabel('x/b',fontsize=8)    
        ax5.set_ylabel('y/b',fontsize=8)        
        ax5.set_xlim((min2_x_beg,max2_x_beg))
        ax5.set_ylim((min2_y_beg,max2_y_beg))
        plt.setp(ax5.get_xticklabels(),visible=False)
        plt.setp(ax5.get_yticklabels(),visible=False)
        #ax5.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
        #ax5.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
        ax5.tick_params(axis='both', which='major', top='off', bottom='off', left='off', right='off')  
        
        textstr = '$R_g=%.1f$\n$S=%.2f$\n$F_{tr}=%.1f$\n$F_{rot}=%.3f$\n$\Omega=%.4f$'%(clinfo.rgysq[index_tracked_clusters[2]], clinfo.pl[index_tracked_clusters[2]], clinfo.st[index_tracked_clusters[2]], clinfo.sw[index_tracked_clusters[2]], clinfo.ens[index_tracked_clusters[2]])
        ax5.text(1.05, 0.95, textstr, transform=ax5.transAxes, fontsize=14,
            verticalalignment='top', bbox=props)
            
        
        # Text
        plt.figtext(subp.beg, subp.ybeg+ax_len+ax_sep/2, '$\\xi_p/L = $' + "{0:.1f}".format(persistence))
        plt.figtext(subp.beg+1.5*ax_sep, subp.ybeg+ax_len+ax_sep/2, '$Pe = $' + "{0:.1f}".format(Pe))
        plt.figtext(subp.beg+3*ax_sep, subp.ybeg+ax_len+ax_sep/2, '$\\tilde F = $' + "{0:.1f}".format(flexure))
        plt.figtext(subp.beg+4.5*ax_sep, subp.ybeg+ax_len+ax_sep/2, '$L = $' + "{0:.1f}".format(body_length))
        plt.figtext(subp.beg+6*ax_sep, subp.ybeg+ax_len+ax_sep/2, '$b = $' + "{0:.1f}".format(B))
        plt.figtext(subp.beg+7.5*ax_sep, subp.ybeg+ax_len+ax_sep/2, '$N = $' + "{0:.1f}".format(M))
        
        plt.figtext(subp.beg+5*ax_len/3, subp.ybeg+ax_len+1.5*ax_sep/2, '$t = $' + str(frame) + ' [timesteps]')
    
    
        plt.savefig(args.savefile + '/climg/'+'frame-'+'{0:05d}'.format(frame)+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
        plt.clf()
        


if __name__ == '__main__':
    main()        
    
    
    

