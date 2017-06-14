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
if os.path.exists(args.savefile+'/hopeful') == False:
    os.mkdir(args.savefile+'/hopeful')
    
# Plot properties
downlim = -5
uplim = max(Lx,Ly)+5
ax_len = 0.8                         # Length of one subplot square box
ax_b = 0.1                            # Beginning/offset of the subplot in the box
ax_sep = 0.1                         # Separation length between two subplots
total_subplots_in_y = 1               # Total number of subplots
tick_interval = int(uplim/5)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

first_time = True
update_limits = True
limit_fac = 150
tracked_clusters = np.zeros(5, dtype=int)
index_tracked_clusters = np.zeros(5, dtype=int)
fig = plt.figure()

# Parameters for the frames

#t0 = 8950000
#t1 = 1210000
#xdown = 800
#xup = 1800
#ydown = 1300
#yup = 2100

#t0 = 53400000
#t1 = 58000000
#xdown = 500
#xup = 1300
#ydown = 0 
#yup = 800

#t0 = 59000000
#t1 = 63000000
#xdown = 2000
#xup = 2740
#ydown = 1000 
#yup = 1740

t0 = 15000000
t1 = 18000000
xdown = 300
xup = 900
ydown = 300 
yup = 900

################################################################################

# Time frame loop
for frame in np.arange(t0,t1,int(nsamp)):
    
    
    print frame 


    #
    # Load and set the data
    #  
    
    # Load particle data  -- 1xL --
    path = datafile + '/beads_' + str(frame) + '.txt'
    p = Particles(path)                    

    
#    # Load cluster information
#    path = datafile + '/cluster_evolution_' + str(frame) + '.txt'
#    clinfo = ClusterInfo(path)

    

################################################################################

    
    #
    # Plot the data
    #   

    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_y)  
   
    if frame == t0 or frame == t1-50000:    
    
        ax1 = subp.addSubplot()  
        quant_steps = 2056
        norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)
        
        line1 = ax1.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, norm=norm)
        ax1.set_xlabel('x/b',fontsize=8)                    
        ax1.set_ylabel('y/b',fontsize=8)
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
     
        ax1.add_patch( Rectangle( (xdown,ydown),xup-xdown,yup-ydown,edgecolor='gray',facecolor='gray',alpha=0.3 ) )
         
    else:
        
        ax1 = subp.addSubplot()  
        quant_steps = 2056
        norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)
        
        line1 = ax1.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, norm=norm)
        ax1.set_xlabel('x/b',fontsize=8)                    
        ax1.set_ylabel('y/b',fontsize=8)
        ax1.set_xlim((xdown,xup))
        ax1.set_ylim((ydown,yup))
        ax1.xaxis.set_ticks( np.linspace(xdown,xup,num=5,endpoint=True) )
        ax1.yaxis.set_ticks( np.linspace(ydown,yup,num=5,endpoint=True) )
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


    # Text
    plt.figtext(subp.beg, subp.ybeg+ax_sep+ax_len, '$\\xi_p/L = $' + "{0:.1f}".format(persistence))
    plt.figtext(subp.beg+1.5*ax_sep, subp.ybeg+ax_sep+ax_len, '$Pe = $' + "{0:.1f}".format(Pe))
    plt.figtext(subp.beg+3*ax_sep, subp.ybeg+ax_sep+ax_len, '$\\tilde F = $' + "{0:.1f}".format(flexure))
    plt.figtext(subp.beg+4.5*ax_sep, subp.ybeg+ax_sep+ax_len, '$L = $' + "{0:.1f}".format(body_length))
    plt.figtext(subp.beg+6*ax_sep, subp.ybeg+ax_sep+ax_len, '$b = $' + "{0:.1f}".format(B))
    plt.figtext(subp.beg+7.5*ax_sep, subp.ybeg+ax_sep+ax_len, '$N = $' + "{0:.1f}".format(M))

    plt.figtext(subp.beg+ax_len/2, subp.ybeg+1.5*ax_sep+ax_len, '$t = $' + str(frame) + ' [timesteps]')

    plt.savefig(args.savefile + '/hopeful/'+'frame-'+'{0:05d}'.format(frame)+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
    plt.clf()
    
################################################################################
    
    
    

