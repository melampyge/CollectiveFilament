#!/usr/bin/python


# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
import os

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
    body_length = B*N
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    flexure = Pe/persistence
    T = totalStep - ti
    nt = T/nsamp
            

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
ax_sep = 0.15                          # Separation length between two subplots
total_subplots_in_y = 2               # Total number of subplots
tick_interval = int(uplim/5)
downlim_zoom = 800
uplim_zoom = 1600
tick_interval_zoom = int((uplim_zoom - downlim_zoom)/5)

# Read the preliminary data of each file   
rho_path = histofile + '/rho.data'
rho_file = open(rho_path, 'r')

vx_path = histofile + '/vx.data'
vx_file = open(vx_path, 'r')
vy_path = histofile + '/vy.data'
vy_file = open(vy_path, 'r')
wor_path = histofile + '/worticity.data'
wor_file = open(wor_path, 'r')
for i in np.arange(7):
    vx_file.readline()
    vy_file.readline()
    wor_file.readline()

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

#xgrid, ygrid = np.mgrid[slice(0., Lx, nx), slice(0., Ly, ny)]
#xgrid = np.transpose(xedges)
#ygrid = np.transpose(yedges)
xlin = np.linspace(0., Lx, nx)
ylin = np.linspace(0., Ly, ny)
xgrid, ygrid = np.meshgrid(xlin, ylin)


# Time frame loop
print "nsteps = ", nsteps
for frame in np.arange(nsteps):
    
    
    print frame 
    
    # Read the lines corresponding to the current frame 
    rho = np.zeros((nx, ny))
    vx = np.zeros((nx, ny))
    vy = np.zeros((nx, ny))
    wor = np.zeros((nx, ny))
    
    line = rho_file.readline()
    line = line.split()
    timestep = float(line[-1])
    
    vx_file.readline()
    vy_file.readline()
    wor_file.readline()
    
    
    for j in np.arange(nx):
        
        rho_line = rho_file.readline()
        rho_line = rho_line.split()
        
        vx_line = vx_file.readline()
        vx_line = vx_line.split()
        
        vy_line = vy_file.readline()
        vy_line = vy_line.split()
        
        wor_line = wor_file.readline()
        wor_line = wor_line.split()
        
        for k in np.arange(ny):
            rho[j, k] = float(rho_line[k])
            vx[j, k] = float(vx_line[k])
            vy[j, k] = float(vy_line[k])
            wor[j, k] = float(wor_line[k])
    
    
    
    #
    # Load and set the data
    #  
    
    # Load particle data  -- 1xL --
    path = clusterfile + '/beads_' + str(int(timestep)) + '.txt'
    p = Particles(path)                    

    

    
    #
    # Plot the data
    #   
   
    # Set up the plot and the subplots structure
    fig = plt.figure()
    
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_y)  
    

    # vx    
    ax0 = subp.addSubplot()
    line0 = ax0.pcolor(xgrid, ygrid, vx.transpose(), cmap='jet')
    ax0.set_xlabel('x/b',fontsize=8)
    ax0.set_ylabel('y/b',fontsize=8)
    ax0.set_title('Binned velocities in x')
    ax0.set_xlim((downlim,uplim))
    ax0.set_ylim((downlim,uplim))
    ax0.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
    ax0.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
    ax0.tick_params(axis='both', which='major',labelsize=8)
    
#    cax0 = plt.axes([subp.xbeg+ax_len+0.01, subp.ybeg, 0.02, ax_len])
#    plt.colorbar(line0, cax=cax0)
#    cax0.set_title('$v_{x}$',fontsize=10)
#    cax0.tick_params(width=0.4,labelsize=6)
    
    # Bead orientations
    ax1 = subp.addSubplot()
    quant_steps = 2056
    norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)
    
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


    # vy
    ax2 = subp.addSubplot()
    ax2.pcolor(xgrid, ygrid, vy.transpose(), cmap='jet')
    ax2.set_xlabel('x/b',fontsize=8)
    ax2.set_title('Binned velocities in y')
    ax2.set_xlim((downlim,uplim))
    ax2.set_ylim((downlim,uplim))
    ax2.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
    ax2.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
    #plt.setp(ax2.get_yticklabels(),visible=False)
    ax2.tick_params(axis='both', which='major', labelsize=8)
  
  
    # Zoomed in bead orientations
    ax3 = subp.addSubplot()  
    ax3.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, norm=norm)
    ax3.set_title('Zoomed in')
    ax3.set_xlim((downlim_zoom,uplim_zoom))
    ax3.set_ylim((downlim_zoom,uplim_zoom))
    ax3.xaxis.set_ticks( np.arange(downlim_zoom,uplim_zoom,tick_interval_zoom) )
    ax3.yaxis.set_ticks( np.arange(downlim_zoom,uplim_zoom,tick_interval_zoom) )
    ax3.tick_params(axis='both', which='major', labelsize=8)
    
    
    # Rotation field of velocity (worticity)
    ax4 = subp.addSubplot()
    ax4.pcolor(xgrid, ygrid, wor.transpose(), cmap='jet')
    ax4.set_xlabel('x/b',fontsize=8)
    ax4.set_title('Rotation field of velocity')
    ax4.set_xlim((downlim,uplim))
    ax4.set_ylim((downlim,uplim))
    ax4.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
    ax4.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
    #plt.setp(ax4.get_yticklabels(),visible=False)
    ax4.tick_params(axis='both', which='major', labelsize=8)
    
    
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
    
    
    # Text
    plt.figtext(subp.beg, subp.ybeg+ax_len+ax_sep, '$\\xi_p/L = $' + "{0:.1f}".format(persistence))
    plt.figtext(subp.beg+1.5*ax_sep, subp.ybeg+ax_len+ax_sep, '$Pe = $' + "{0:.1f}".format(Pe))
    plt.figtext(subp.beg+3*ax_sep, subp.ybeg+ax_len+ax_sep, '$\\tilde F = $' + "{0:.1f}".format(flexure))
    plt.figtext(subp.beg+4.5*ax_sep, subp.ybeg+ax_len+ax_sep, '$L = $' + "{0:.1f}".format(body_length))
    plt.figtext(subp.beg+6*ax_sep, subp.ybeg+ax_len+ax_sep, '$b = $' + "{0:.1f}".format(B))
    plt.figtext(subp.beg+7.5*ax_sep, subp.ybeg+ax_len+ax_sep, '$N = $' + "{0:.1f}".format(M))

    
    plt.figtext(subp.beg+ax_len, subp.ybeg+ax_len+1.5*ax_sep, '$t = $' + str(timestep) + ' [timesteps]')


    plt.savefig(args.savefile + '/histo/'+'frame-'+'{0:05d}'.format(int(timestep))+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
    plt.clf()
    
    
    
    

