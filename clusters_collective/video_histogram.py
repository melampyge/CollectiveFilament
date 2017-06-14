#!/usr/bin/python


# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys
import bottleneck as bot
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
    T = totalStep - ti
    nt = T/nsamp
    box_area = Lx*Ly
    body_length = B*N
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    flexure = Pe/persistence
    print nsamp, ti, totalStep, T, totalStep-ti
        

#
# Class definitions
#

# Particle data
class Particles:
    
    def __init__(self, path):
        file = np.transpose(np.loadtxt(path, dtype=float))
        self.xi = file[0]/B                 # Image particle positions in x
        self.yi = file[1]/B                 # Image particle positions in y 
        self.phi = file[2] + np.pi          # Bead orientation 
        self.cidx = file[3]                 # Cluster index
        
# Cluster sizes
class ClusterSize:
    
    def __init__(self, path):
        file = np.loadtxt(path, dtype=float)
        self.cs = file                     # Cluster sizes
        self.max = np.size(self.cs)
        self.cmass = np.zeros((self.max,1))
        
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
datafile = args.initfile+"/HISTOGRAMS/tables"
if os.path.exists(args.savefile+'/histo') == False:
    os.mkdir(args.savefile+'/histo')

# Plot properties
downlim = -5
uplim = max(Lx,Ly)+5
ax_len = 0.4                      # Length of one subplot square box
ax_b = 0.1                        # Beginning/offset of the subplot in the box
ax_sep = ax_len/4                 # Separation length between two subplots
total_subplots_in_y = 2           # Total number of subplots
tick_interval = int(uplim/5)
    

# Read the data    
rho_path = datafile + '/rho.data'
rho_file = open(rho_path, 'r')

vx_path = datafile + '/vx.data'
vx_file = open(vx_path, 'r')
vy_path = datafile + '/vy.data'
vy_file = open(vy_path, 'r')
wor_path = datafile + '/worticity.data'
wor_file = open(wor_path, 'r')
for i in np.arange(7):
    vx_file.readline()
    vy_file.readline()
    wor_file.readline()

line = rho_file.readline()   
line = line.split()
nsteps = int(line[-1])

line = rho_file.readline()   
line = line.split()
nx = int(line[-1])

line = rho_file.readline()   
line = line.split()
ny = int(line[-1])

xedges = np.zeros((nx+1, ny+1))
yedges = np.zeros((nx+1, ny+1))

rho_file.readline()

line = rho_file.readline()   
line = line.split()
for i in range(nx+1):
    xedges[i] = float(line[i])

rho_file.readline()

line = rho_file.readline()   
line = line.split()
for i in range(ny+1):
    yedges[i] = float(line[i])

    
# Time frame loop
print 'nsteps = ', nsteps
for i in np.arange(nsteps):
    
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
            
    
    # Plot the data
   
    # Set up the plot and the subplots structure
    fig = plt.figure()
    
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_y)  
    
    
    ax0 = subp.addSubplot()
    ax0.imshow(rho.transpose(), origin = 'lower', interpolation = 'None')
    
    ax1 = subp.addSubplot()
    ax1.imshow(vx.transpose(), origin = 'lower', interpolation = 'None')

    ax2 = subp.addSubplot()
    ax2.imshow(vy.transpose(), origin = 'lower', interpolation = 'None')

    ax3 = subp.addSubplot()
    ax3.imshow(wor.transpose(), origin = 'lower', interpolation = 'None')

    plt.savefig(args.savefile+'/histo/'+'frame-'+'{0:05f}'.format(timestep)+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
    plt.clf()
        
rho_file.close()    
vx_file.close()
vy_file.close()
wor_file.close()

  
#    #
#    # Load and set the data
#    #  
#    
#    path = datafile + 'rho.data'
#    
#    
#    # Load particle data  -- 1xL --
#    path = datafile + '/beads_' + str(frame) + '.txt'
#    p = Particles(path)
#
#    
#    # Load cluster sizes
#    path = datafile + '/cluster_sizes_' + str(frame) + '.txt'
#    clsizes = ClusterSize(path)
#    
#    # Count the cluster indices
#    cnt = np.zeros(len(p.cidx))
#    for pi in p.cidx:
#        cnt[pi] += 1
#    
#    # Get the largest 5 clusters    
#    maxs = -bot.partsort(-cnt, 5)[:5]
#    
#    # Index the largest 5 clusters to visualize them
#    for i, pc in enumerate(p.cidx):
#        if cnt[pc] < maxs[-1]:
#            p.cidx[i] = 200
#        elif cnt[pc] == maxs[0]:
#            p.cidx[i] = 10
#        elif cnt[pc] == maxs[1]:
#            p.cidx[i] = 30
#        elif cnt[pc] == maxs[2]:
#            p.cidx[i] = 60
#        elif cnt[pc] == maxs[3]:
#            p.cidx[i] = 90
#        elif cnt[pc] == maxs[4]:
#            p.cidx[i] = 120
#    
#    
#    # Get the cluster mass
#    clsizes.getMass()
#    
#    #
#    # Plot the data
#    #   
#   
#    # Set up the plot and the subplots structure
#    fig = plt.figure()
#    
#    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_y)  
#    
#    
#    # Histogram of bead orientations 
#    ax0 = subp.addSubplot()
#    weights = np.ones_like(p.phi)/len(p.phi)
#    ax0.hist(p.phi, bins=20, weights=weights, histtype='bar', color='red')
#    ax0.set_xlabel('Orient.',fontsize=8)
#    ax0.set_ylabel('Dist.',fontsize=8)
#    ax0.set_title('Bead orient. dist.')
#    ax0.set_xlim((0, 2.1*np.pi))
#    ax0.set_ylim((0, 1.1))
#    ax0.xaxis.set_ticks( np.arange(0,2.1*np.pi,np.pi) )
#    labels = [item.get_text() for item in ax0.get_xticklabels()]
#    labels[0] = '0'
#    labels[1] = '$\\pi$'
#    labels[2] = '$2\\pi$'
#    ax0.set_xticklabels(labels)
#    ax0.yaxis.set_ticks( np.arange(0,1.1,0.3) )
#    ax0.tick_params(axis='both', which='major', labelsize=8)    
#    
#    
#    
#    # Bead orientations
#    ax1 = subp.addSubplot()
#    norm = mpl.colors.Normalize(0.0, 2*np.pi)
#    line1 = ax1.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.hsv, 
#                edgecolors='None', alpha=0.7, vmin=0.0, vmax=2*np.pi)
#    ax1.set_xlabel('x/b',fontsize=8)
#    ax1.set_ylabel('y/b',fontsize=8)
#    ax1.set_title('Bead Orientations')
#    ax1.set_xlim((downlim,uplim))
#    ax1.set_ylim((downlim,uplim))
#    #plt.setp(ax1.get_xticklabels(),visible=False)
#    ax1.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
#    ax1.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
#    ax1.tick_params(axis='both', which='major', labelsize=8)
#    
#    #cax1 = plt.axes([subp.xbeg+ax_len+0.02, subp.ybeg, 0.8, 0.8], projection='polar')
#    cax1 = plt.axes([subp.xbeg+ax_len+0.01, subp.ybeg+ax_len/3, ax_len/4.6, ax_len/4.6], projection='polar')
#    cax1._direction = 2*np.pi
#    #plt.colorbar(line1,cax=cax1,ticks=[-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
#    quant_steps = 2056
#    cb1 = mpl.colorbar.ColorbarBase(cax1, cmap=plt.cm.get_cmap('hsv',quant_steps), \
#        norm=norm, orientation='horizontal')    
#    plt.colorbar(line1, cax=cax1)
#    cb1.outline.set_visible(False)
#    cax1.set_axis_off()
#    cax1.set_title('$\\phi$',fontsize=10)
#    cax1.set_rlim([-1,1])
#    #cax1.set_yticklabels(['$-\\pi$', '$-\\pi/2$', '$0$', '$\\pi/2$', '$\\pi$'])
#
#
#
#    # Histogram of cluster sizes 
#    ax2 = subp.addSubplot()
#    weights = np.ones_like(clsizes.cmass)/len(clsizes.cmass)
#    ax2.hist(clsizes.cmass, log=True, bins=20, weights=weights, histtype='bar', color='red')
#    ax2.set_xlabel('Cl. mass',fontsize=8)
#    ax2.set_ylabel('log(Dist.)',fontsize=8)
#    ax2.set_title('Cluster mass dist.')
#    ax2.set_xlim((-1,M+500))
#    ax2.set_ylim((-3,3))
#    ax2.xaxis.set_ticks( np.arange(0,M+500,1000) )
#    ax2.yaxis.set_ticks( np.arange(10**(-3),10**0) )
#    ax2.tick_params(axis='both', which='major', labelsize=8)   
#
#    
#  
#    # Cluster id
#    ax3 = subp.addSubplot()  
#    ax3.scatter(p.xi, p.yi, s=1, c=p.cidx, cmap=plt.cm.Set1, edgecolors='None', alpha=0.7)
#    ax3.set_xlabel('x/b',fontsize=8)
#    #ax3.set_ylabel('y/b',fontsize=8)
#    ax3.set_title('Cluster id')
#    ax3.set_xlim((downlim,uplim))
#    ax3.set_ylim((downlim,uplim))
#    #plt.setp(ax3.get_xticklabels(),visible=False)
#    plt.setp(ax3.get_yticklabels(),visible=False)
#    ax3.xaxis.set_ticks( np.arange(0,uplim,tick_interval) )
#    #ax3.yaxis.set_ticks( np.arange(0,uplim,tick_interval) )
#    ax3.tick_params(axis='both', which='major', labelsize=8)
#    
#    # Text
#    plt.figtext(subp.beg, subp.ybeg+ax_len+ax_sep, '$\\xi_p/L = $' + "{0:.1f}".format(persistence))
#    plt.figtext(subp.beg+1.5*ax_sep, subp.ybeg+ax_len+ax_sep, '$Pe = $' + "{0:.1f}".format(Pe))
#    plt.figtext(subp.beg+3*ax_sep, subp.ybeg+ax_len+ax_sep, '$\\tilde F = $' + "{0:.1f}".format(flexure))
#    plt.figtext(subp.beg+4.5*ax_sep, subp.ybeg+ax_len+ax_sep, '$L = $' + "{0:.1f}".format(body_length))
#    plt.figtext(subp.beg+6*ax_sep, subp.ybeg+ax_len+ax_sep, '$b = $' + "{0:.1f}".format(B))
#    plt.figtext(subp.beg+7.5*ax_sep, subp.ybeg+ax_len+ax_sep, '$N = $' + "{0:.1f}".format(M))
#
#    
#    plt.figtext(subp.beg+ax_len, subp.ybeg+ax_len+1.5*ax_sep, '$t = $' + str(frame) + ' [timesteps]')
#
#
#    plt.savefig(args.savefile + '/img/'+'frame-'+'{0:05d}'.format(frame)+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
#    plt.clf()
#    
#    
#    
#    

