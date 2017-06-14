#!/usr/bin/python


# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
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

# Particle data
class Particles:
    
    def __init__(self, path):
        file = np.transpose(np.loadtxt(path, dtype=float))
        self.xi = file[0]/B                 # Image particle positions in x
        self.yi = file[1]/B                 # Image particle positions in y 
        self.phi = file[2]                  # Bead orientation 
        self.cidx = file[3]                 # Cluster index
        #self.cidx = np.zeros(len(self.xi))
        
        
# Cluster sizes
class ClusterSize:
    
    def __init__(self, path):
        f = np.loadtxt(path, dtype=float)
        self.cs = f                                  # Cluster sizes
        self.totcs = self.cs                         # Total number of clusters
        
    def getMass(self):
        for i in np.arange(self.max):
            self.cmass[i] = i*self.cs[i]
         
         
# Cluster properties
class ClusterProps:
    
    def __init__(self, path):
        f = np.transpose(np.loadtxt(path, dtype=float))
        self.idx = f[0]
        self.comx = f[1]/B
        self.comy = f[2]/B
        self.str = f[3]
        self.sw = f[4]
        self.ens = f[5]
        
# Cluster information
class ClusterInfo:
    
    def __init__(self, path):
        f = open(path, 'r')
        self.fils = []
        self.idx = []
        self.size = []
        self.comx = []
        self.comy = []
        self.rgysq = []
        self.pl = []
        self.st = []
        self.sw = []
        self.ens = []
        for line in f:
            A = line.split()
            self.idx.append(int(A[0]))
            siz = int(A[1])
            self.size.append(siz)
            fil_list = np.zeros((siz))
            self.comx.append(float(A[2])/B)
            self.comy.append(float(A[3])/B)
            self.rgysq.append(np.sqrt(float(A[4])))
            self.pl.append(float(A[5]))
            self.st.append(float(A[6]))
            self.sw.append(float(A[7]))
            self.ens.append(float(A[8]))
            for el in range(siz):
                fil_list[el] = int(A[9+el])
            self.fils.append(fil_list)
   
        

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
      
      
# Argument parsing (command line options)
parser = argparse.ArgumentParser()
parser.add_argument("initfile", help="File containing initial simulation data")
parser.add_argument("savefile", help="File in which figures should be saved")
args = parser.parse_args()

# Load saved preliminary simulation data into relevant variables
loadSimData(args.initfile+"/init_info.txt")
datafile = args.initfile+"/CLUSTER"

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
tracked_clusters = np.zeros(1, dtype=int)
index_tracked_clusters = np.zeros(1, dtype=int)
fig = plt.figure()

################################################################################

# Time frame loop
frame = 5100000
#frame = 5150000
#frame = 5050000

print frame 

#
# Load and set the data
#  

# Load particle data  -- 1xL --
path = datafile + '/beads_' + str(frame) + '.txt'
p = Particles(path)                    


# Load cluster information
path = datafile + '/cluster_evolution_' + str(frame) + '.txt'
clinfo = ClusterInfo(path)

tracked_clusters = [9]
index_tracked_clusters = []
cnt = 0
for i, ind in enumerate(clinfo.idx):
    if ind == tracked_clusters[0]:
        index_tracked_clusters.append(i)
        cnt += 1
        print ' BURADA '
        if cnt == 2:
            break
else:
    index_tracked_clusters = [9, 1]
    
        
#found = findArrayIndex(clinfo.idx, tracked_clusters, index_tracked_clusters)  
#findIndex(clinfo.size, clinfo.idx, tracked_clusters, index_tracked_clusters)   
        
 
# Highlight the indexed clusters
clidx = []
for i in range(L):
    clidx.append('gray')
    
flags = ['False', 'False', 'False', 'False', 'False']
c0_x, c0_y, c0_p, min0_x, max0_x, min0_y, max0_y, flags[0] = colorToBeads(clidx, clinfo.fils[index_tracked_clusters[0]], 'maroon', p.xi, p.yi, p.phi)
c1_x, c1_y, c1_p, min1_x, max1_x, min1_y, max1_y, flags[1] = colorToBeads(clidx, clinfo.fils[index_tracked_clusters[1]], 'violet', p.xi, p.yi, p.phi)


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
    
  
# Cluster 1 / Maroon
ax3 = subp.addSubplot()  
line3 = ax3.scatter(c0_x, c0_y, s=1, c=c0_p, cmap=plt.cm.get_cmap('hsv',quant_steps), 
            edgecolors='None', alpha=1.0, vmin=-np.pi, vmax=np.pi, norm=norm)  
ax3.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
            edgecolors='None', alpha=0.1, vmin=-np.pi, vmax=np.pi, norm=norm)                               
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
line4 = ax4.scatter(c1_x, c1_y, s=1, c=c1_p, cmap=plt.cm.get_cmap('hsv',quant_steps), 
            edgecolors='None', alpha=1.0, vmin=-np.pi, vmax=np.pi, norm=norm)  
ax4.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
            edgecolors='None', alpha=0.1, vmin=-np.pi, vmax=np.pi, norm=norm)                               
ax4.plot(clinfo.comx[index_tracked_clusters[1]],clinfo.comy[index_tracked_clusters[1]],'*',color='violet')                
ax4.set_title('Violet')
ax4.set_xlabel('x/b',fontsize=8)    
ax4.set_ylabel('y/b',fontsize=8)        
ax4.set_xlim((min1_x_beg,max1_x_beg))
ax4.set_ylim((min1_y_beg,max1_y_beg))
plt.setp(ax4.get_xticklabels(),visible=False)
plt.setp(ax4.get_yticklabels(),visible=False)
#ax4.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
#ax4.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
ax4.tick_params(axis='both', which='major', top='off', bottom='off', left='off', right='off')  

textstr = '$R_g=%.1f$\n$S=%.2f$\n$F_{tr}=%.1f$\n$F_{rot}=%.3f$\n$\Omega=%.4f$'%(clinfo.rgysq[index_tracked_clusters[1]], clinfo.pl[index_tracked_clusters[1]], clinfo.st[index_tracked_clusters[1]], clinfo.sw[index_tracked_clusters[1]], clinfo.ens[index_tracked_clusters[1]])
ax4.text(1.05, 0.95, textstr, transform=ax4.transAxes, fontsize=14,
    verticalalignment='top', bbox=props)
 

plt.savefig(args.savefile + '/climg/'+'frame-'+'{0:05d}'.format(frame)+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
plt.clf()

################################################################################


    

