#!/usr/bin/python


# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
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

    
################################################################################
                
#
# Class definitions
#

################################################################################

class Particles:
    """ Load particle data"""
    
    def __init__(self, path):
        file = np.transpose(np.loadtxt(path, dtype=float))
        self.xi = file[0]/B                 # Image particle positions in x
        self.yi = file[1]/B                 # Image particle positions in y 
        self.phi = file[2]                  # Bead orientation 
        self.cidx = file[3]                 # Cluster index
        #self.cidx = np.zeros(len(self.xi))
                
################################################################################
                
class Subplots:
    """ Arrange subplot grid structure (square box is assumed) """
    
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
clusterfile = args.initfile+"/CLUSTER"
histofile = args.initfile+"/HISTOGRAMS/tables"
if os.path.exists(args.savefile+'/giant') == False:
    os.mkdir(args.savefile+'/giant')

# Plot properties
quant_steps = 2056
norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)
t0 = 10000000
t1 = 17000000
xdown = 0
xup = 1000
ydown = 0 
yup = 1000
totdown = 0
totup = Lx
fig = plt.figure()

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

break_flag = False

# Time frame loop
for frame in np.arange(nsteps):
    
    print frame 
    
    # Read the lines corresponding to the current frame 
    rho = np.zeros((nx, ny))
    
    line = rho_file.readline()
    line = line.split()
    timestep = float(line[-1])
 
    if timestep < t0:
        break_flag = True
    else:
        break_flag = False
        
    if timestep > t1:
        break
    
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
                      

################################################################################

    
    #
    # Plot the data
    #   


    ## Add the density field subplot
    ax_dens = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    line_dens = ax_dens.pcolor(xgrid, ygrid, rho.transpose(), cmap='jet', vmin=0, vmax=0.7)
    ax_dens.set_title('Density field', fontsize=30)
    ax_dens.set_xlim((totdown,totup))
    ax_dens.set_ylim((totdown,totup))
    ax_dens.set_xlabel('x/b', fontsize=20)
    ax_dens.set_ylabel('y/b', fontsize=20)    
    ax_dens.xaxis.set_ticks(np.linspace(totdown,totup,num=5,endpoint=True))
    ax_dens.yaxis.set_ticks(np.linspace(totdown,totup,num=5,endpoint=True))
    ax_dens.tick_params(axis='both', which='major', labelsize=15)

    cax_dens = plt.axes([0.1+0.8+0.05, 0.1, 0.05, 0.8])
    plt.colorbar(line_dens, cax=cax_dens, ticks=[0, 0.7])
    #cax_dens.set_yticks([0, 0.8])
    cax_dens.set_yticklabels(['0', '0.7'])  
    cax_dens.tick_params(labelsize=20)     
 
    
    ## Add the overview subplot
    ax_overview = fig.add_axes([0.1, 1.2, 0.8, 0.8])
    line_overview = ax_overview.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), \
                edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, norm=norm)
    ax_overview.set_title('Overview', fontsize=30) 
    ax_overview.set_xlim((totdown,totup))
    ax_overview.set_ylim((totdown,totup))
    ax_overview.set_xlabel('x/b', fontsize=20)
    ax_overview.set_ylabel('y/b', fontsize=20)    
    ax_overview.xaxis.set_ticks(np.linspace(totdown,totup,num=5,endpoint=True))
    ax_overview.yaxis.set_ticks(np.linspace(totdown,totup,num=5,endpoint=True))
    ax_overview.tick_params(axis='both', which='major', labelsize=15)

    ax_overview.add_patch( Rectangle( (xdown,ydown),xup-xdown,yup-ydown,edgecolor='yellow',facecolor='yellow',alpha=0.3 ) )
    
    cax_overview = plt.axes([0.1+0.8+0.01, 0.1+0.8+0.3+0.8/3, 0.8/4.6, 0.8/4.6], projection='polar')
    xval = np.arange(-np.pi, np.pi, 0.01)
    yval = np.ones_like(xval)
    cax_overview.scatter(xval, yval, c=xval, s=300, cmap=plt.cm.get_cmap('hsv',quant_steps), norm=norm, linewidths=0)
    cax_overview.set_yticks([])
    cax_overview.set_xticks([])
    cax_overview.set_title('$\\phi$',fontsize=30)
    cax_overview.set_rlim([-1,1])
    cax_overview.set_axis_off()   
    
    
    ## Add the zoom subplot
    ax_zoom = fig.add_axes([1.2, 0.1, 1.9, 1.9])
    line_zoom = ax_zoom.scatter(p.xi, p.yi, s=1, c=mol, cmap=plt.cm.get_cmap('jet'), \
                edgecolors='None', alpha=1.)
    ax_zoom.set_title('Zoom (without color code)', fontsize=30) 
    ax_zoom.set_xlim((xdown,xup))
    ax_zoom.set_ylim((ydown,yup))
#    ax_zoom.set_xlabel('x/b', fontsize=20)
#    ax_zoom.set_ylabel('y/b', fontsize=20)    
    ax_zoom.xaxis.set_ticks(np.linspace(xdown,xup,num=5,endpoint=True))
    ax_zoom.yaxis.set_ticks(np.linspace(ydown,yup,num=5,endpoint=True))
    ax_zoom.tick_params(axis='both', which='major', labelsize=15)

    # Text
    sep_fac = 0.3
    plt.figtext(0.1, 0.1+1.9+0.2, '$\\xi_p/L = $' + "{0:.1f}".format(persistence), fontsize=20)
    plt.figtext(0.1+sep_fac, 0.1+1.9+0.2, '$Pe = $' + "{0:.1f}".format(Pe), fontsize=20)
    plt.figtext(0.1+2*sep_fac, 0.1+1.9+0.2, '$\\tilde F = $' + "{0:.1f}".format(flexure), fontsize=20)
    plt.figtext(0.1+3*sep_fac, 0.1+1.9+0.2, '$\\phi = $' + '0.08', fontsize=20)
    plt.figtext(0.1+4*sep_fac, 0.1+1.9+0.2, '$L = $' + "{0:.1f}".format(body_length), fontsize=20)
    plt.figtext(0.1+5*sep_fac, 0.1+1.9+0.2, '$N = $' + "{0:.1f}".format(M), fontsize=20)

    plt.figtext(0.1+6*sep_fac, 0.1+1.9+0.2, '$t = $' + str(timestep) + ' [timesteps]', fontsize=20)
   
    plt.savefig(args.savefile + '/giant/'+'frame-'+'{0:05d}'.format(frame)+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
    plt.clf()
    
################################################################################
    
    
    

