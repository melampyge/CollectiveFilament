#!/usr/bin/python


# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
import os
from matplotlib.patches import Rectangle


#
# Function definitions
#
 

# Load initial simulation data      
def loadSimData(datafile):
    
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


def gen_mol_info(natoms, nfil):
    """ generate information about molecular ID"""
    mol = np.zeros((natoms), dtype = np.int32)
    nmol = natoms/nfil
    k = 0
    for i in range(nmol):
        for j in range(nfil):
            mol[k] = i
            k = k + 1
    return mol,nmol
            

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
        
# Cluster sizes
class ClusterSize:
    
    def __init__(self, path):
        file = np.loadtxt(path, dtype=float)
        self.cs = file                     # Cluster sizes
        #self.max = np.size(self.cs)
        #self.cmass = np.zeros((self.max,1))
        
    def getMass(self):
        for i in np.arange(self.max):
            self.cmass[i] = i*self.cs[i]

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
        

# Argument parsing (command line options)
parser = argparse.ArgumentParser()
parser.add_argument("initfile", help="File containing initial simulation data")
parser.add_argument("savefile", help="File in which figures should be saved")
args = parser.parse_args()

# Load saved preliminary simulation data into relevant variables
loadSimData(args.initfile+"/init_info.txt")
clusterfile = args.initfile+"/CLUSTER"
histofile = args.initfile+"/HISTOGRAMS/tables"
if os.path.exists(args.savefile+'/histo') == False:
    os.mkdir(args.savefile+'/histo')

# Plot properties
downlim = -5
uplim = max(Lx,Ly)+5
ax_len = 0.38                         # Length of one subplot square box
ax_b = 0.1                            # Beginning/offset of the subplot in the box
ax_sep = 0.15                         # Separation length between two subplots
total_subplots_in_y = 2               # Total number of subplots
tick_interval = int(uplim/5)
downlim_zoom = 800
uplim_zoom = 1600
tick_interval_zoom = int((uplim_zoom - downlim_zoom)/4)

quant_steps = 2056
norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)

# Read the preliminary data of each file   
rho_path = histofile + '/rho.data'
rho_file = open(rho_path, 'r')

# Read total number of steps
line = rho_file.readline()   
line = line.split()
nsteps = int(line[-1])

# Read number of bins in x and y directions
line = rho_file.readline()   
line = line.split()
nx = int(line[-1])

line = rho_file.readline()   
line = line.split()
ny = int(line[-1])

xedges = np.zeros((nx, ny))
yedges = np.zeros((nx, ny))

# Read the first frame
rho_file.readline()


line = rho_file.readline()   
line = line.split()
for i in range(nx):
    xedges[i] = float(line[i])

rho_file.readline()

line = rho_file.readline()   
line = line.split()
for i in range(ny):
    yedges[i] = float(line[i])

xlin = np.linspace(0., Lx, nx)
ylin = np.linspace(0., Ly, ny)
xgrid, ygrid = np.meshgrid(xlin, ylin)


# Time frame loop
print "nsteps = ", nsteps
break_flag = False
for frame in np.arange(nsteps):
    
    
    print frame 
    
    # Read the lines corresponding to the current frame 
    rho = np.zeros((nx, ny))
    
    line = rho_file.readline()
    line = line.split()
    timestep = float(line[-1])
    if timestep < 5000000:
        break_flag = True
    else:
        break_flag = False
    
    
    for j in np.arange(nx):
        
        rho_line = rho_file.readline()
        rho_line = rho_line.split()
              
        for k in np.arange(ny):
            rho[j, k] = float(rho_line[k])

    if break_flag:
        continue
    
    #
    # Load and set the data
    #  
    
    # Load particle data  -- 1xL --
    path = clusterfile + '/beads_' + str(int(timestep)) + '.txt'
    p = Particles(path)                    
    mol, nmol = gen_mol_info(int(L), int(N))
    

    
    #
    # Plot the data
    #   
   
    # Set up the plot and the subplots structure
    fig = plt.figure()
    
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_y)  

    
    # (0,0) to (500,0)    
    ax0 = subp.addSubplot()
    ax0.scatter(p.xi, p.yi, s=1, c=mol, cmap=plt.cm.get_cmap('jet'), 
                edgecolors='None', alpha=1)
    ax0.set_title('Yellow region')
    ax0.set_xlim((0,500))
    ax0.set_ylim((0,500))
    ax0.xaxis.set_ticks( np.linspace(0,500,num=6,endpoint=True) )
    ax0.yaxis.set_ticks( np.linspace(0,500,num=6,endpoint=True) )
    ax0.tick_params(axis='both', which='major', labelsize=8)
    ax0.set_xlabel('x/b',fontsize=8)
    ax0.set_ylabel('y/b',fontsize=8)
    
    
    # Bead orientations
    ax1 = subp.addSubplot()
    line1 = ax1.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, norm=norm)
    ax1.set_xlabel('x/b',fontsize=8)
    ax1.set_ylabel('y/b',fontsize=8)
    ax1.set_title('Bead Orientations')
    ax1.set_xlim((downlim,uplim))
    ax1.set_ylim((downlim,uplim))
    #plt.setp(ax1.get_xticklabels(),visible=False)
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
    
#    ax1.add_patch( Rectangle( (0,0),800,800,fill=None,edgecolor='black',alpha=1 ) )
#    ax1.add_patch( Rectangle( (800,0),800,800,fill=None,edgecolor='gray',alpha=1 ) )
#    ax1.add_patch( Rectangle( (800,800),400,400,fill=None,edgecolor='blue',alpha=1 ) )
#    ax1.add_patch( Rectangle( (1600,0),200,200,fill=None,edgecolor='red',alpha=1 ) )
    
    ax1.add_patch( Rectangle( (0,0),500,500,edgecolor='black',facecolor='yellow',alpha=0.3 ) )
    ax1.add_patch( Rectangle( (500,0),300,300,edgecolor='gray',facecolor='gray',alpha=0.3 ) )
    ax1.add_patch( Rectangle( (500,300),300,300,edgecolor='blue',facecolor='blue',alpha=0.3 ) )
    ax1.add_patch( Rectangle( (800,0),200,200,edgecolor='red',facecolor='red',alpha=0.3 ) )



    # Zoom 2 -- (500, 800) to (0, 300)
    ax2 = subp.addSubplot()
    ax2.scatter(p.xi, p.yi, s=1, c=mol, cmap=plt.cm.get_cmap('jet'), 
                edgecolors='None', alpha=0.7)
    ax2.set_title('Gray region')
    ax2.set_xlim((500,800))
    ax2.set_ylim((0,300))
    ax2.xaxis.set_ticks( np.linspace(500,800,num=6,endpoint=True) )
    ax2.yaxis.set_ticks( np.linspace(0,300,num=6,endpoint=True) )
    ax2.tick_params(axis='both', which='major', labelsize=8)
    ax2.set_xlabel('x/b',fontsize=8)
  
  
    # Zoom 3 -- (500, 800) to (300, 600)
    ax3 = subp.addSubplot()  
    ax3.scatter(p.xi, p.yi, s=1, c=mol, cmap=plt.cm.get_cmap('jet'), 
                edgecolors='None', alpha=0.7)
    ax3.set_title('Blue region')
    ax3.set_xlim((500,800))
    ax3.set_ylim((300,600))
    ax3.xaxis.set_ticks( np.linspace(500,800,num=6,endpoint=True) )
    ax3.yaxis.set_ticks( np.linspace(300,600,num=6,endpoint=True) )
    ax3.tick_params(axis='both', which='major', labelsize=8)
    
    
    # More zoomed in
    ax4 = subp.addSubplot()
    ax4.scatter(p.xi, p.yi, s=1, c=mol, cmap=plt.cm.get_cmap('jet'), 
            edgecolors='None', alpha=1.2)
    ax4.set_title('Red region')
    ax4.set_xlim((800,1000))
    ax4.set_ylim((0,200))
    ax4.xaxis.set_ticks( np.linspace(800,1000,num=6,endpoint=True) )
    ax4.yaxis.set_ticks( np.linspace(0,200,num=6,endpoint=True) )
    ax4.tick_params(axis='both', which='major', labelsize=8) 
    ax4.set_xlabel('x/b',fontsize=8)
    
    
    # Density
    ax5 = subp.addSubplot()
    line2 = ax5.pcolor(xgrid, ygrid, rho.transpose(), cmap='jet')
    ax5.set_title('Binned density field')
    ax5.set_xlim((downlim,uplim))
    ax5.set_ylim((downlim,uplim))
    ax5.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
    ax5.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
    #plt.setp(ax5.get_xticklabels(),visible=False)
    ax5.tick_params(axis='both', which='major', labelsize=8)
    

#    ax0.text(800, 400, "Test", size=10, va="center", ha="center", rotation=180,
#             bbox=dict(boxstyle="angled,pad=0.5", alpha=0.2))
#    ax2.text(1600, 400, "Test", size=10, va="center", ha="center", rotation=180,
#        bbox=dict(boxstyle="angled,pad=0.5", alpha=0.2))
#    ax2.text(1200, 800, "Test", size=10, va="center", ha="center", rotation=270,
#        bbox=dict(boxstyle="angled,pad=0.5", alpha=0.2))
    
    
    # Text
    plt.figtext(subp.beg, subp.ybeg+ax_len+ax_sep, '$\\xi_p/L = $' + "{0:.1f}".format(persistence))
    plt.figtext(subp.beg+1.5*ax_sep, subp.ybeg+ax_len+ax_sep, '$Pe = $' + "{0:.1f}".format(Pe))
    plt.figtext(subp.beg+3*ax_sep, subp.ybeg+ax_len+ax_sep, '$\\tilde F = $' + "{0:.1f}".format(flexure))
    plt.figtext(subp.beg+4.5*ax_sep, subp.ybeg+ax_len+ax_sep, '$L = $' + "{0:.1f}".format(body_length))
    plt.figtext(subp.beg+6*ax_sep, subp.ybeg+ax_len+ax_sep, '$b = $' + "{0:.1f}".format(B))
    plt.figtext(subp.beg+7.5*ax_sep, subp.ybeg+ax_len+ax_sep, '$N = $' + "{0:.1f}".format(M))

    
    plt.figtext(subp.beg+ax_len+ax_len/2, subp.ybeg+ax_len+1.5*ax_sep, '$t = $' + str(timestep) + ' [timesteps]')


    plt.savefig(args.savefile + '/histo/'+'frame-'+'{0:05d}'.format(int(timestep))+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
    plt.clf()
    
    
    
    

