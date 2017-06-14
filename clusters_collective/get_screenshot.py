
# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
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
    
    return
    
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
        self.tot = t        # Total number of subplots in the x direction
        
    # Add a subplot in the grid structure
    def addSubplot(self):
        
        # Increase the number of subplots in the figure
        self.totcnt += 1
        
        # Indices of the subplot in the figure
        self.nx = self.totcnt%(self.tot)
        self.ny = self.totcnt/(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])
        
################################################################################
      
def main():
      
    ## Argument parsing (command line options)
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--density", help="Density")
    parser.add_argument("-k", "--kappa", help="Bending rigidity")
    parser.add_argument("-f", "--fp", help="Propulsion strength")
    parser.add_argument("-t", "--time", help="Timestep")    
    parser.add_argument("-s", "--save", help="Save options", action="store_true")      
    parser.add_argument("-i", "--info", help="Get info", action="store_true")          
    args = parser.parse_args()
 
    ## Load saved preliminary simulation data into relevant variables
#    basepath = "/local/duman/SIMULATIONS/many_polymers_5"
#    databasepath = basepath+"/density_"+args.density+"/kappa_"+args.kappa+"/fp_"+args.fp
#    #loadSimData(databasepath+"/init_info.txt")
#    loadSimData(basepath+"/density_"+"0.4"+"/kappa_"+"5.0"+"/fp_"+"1.2"+"/init_info.txt")  
    
    basepath = "/local/duman/SIMULATIONS/long_filaments" 
    databasepath = basepath+"/density_"+args.density+"/kappa_"+args.kappa+"/fp_"+args.fp
    loadSimData(databasepath+"/init_info.txt")  
    
    datapath = databasepath+"/CLUSTER"
    savebasepath = "/usr/users/iff_th2/duman/RolfData/many_polymers_5"
    savepath = savebasepath+"/density_"+args.density+"/kappa_"+args.kappa+"/fp_"+args.fp
    savepath += "/screenshot"
    os.system("mkdir -p " + savepath)
        
    ## Display information if need be    
    if args.info:
        print "First timestep is : ", int(ti)
        print "Last timestep is : ", int(T)
        print "Middle timestep is : ", int((T-ti)/2)
        
        return
        
        
    ## Plot properties
    downlim = -5
    uplim = max(Lx,Ly)+5
    ax_len = 0.4                          # Length of one subplot square box
    ax_b = 0.1                            # Beginning/offset of the subplot in the box
    ax_sep = 0.2                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots in the x direction
    tick_num = 4                          # Number of ticks
    fig = plt.figure()
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x)  
    timestep = int(args.time)
    
    ################################################################################
    
    ## Plot the data
    print 'Plotting'
    print 'Frame = ', args.time
    print 'File = ', datapath
    
    ## Load particle data  -- 1xL --
    path = datapath + '/cluster_' + args.time + '.hdf5'           
    p = Particles(path)
    
    ## Bead orientations
    ax1 = subp.addSubplot()  
    quant_steps = 2056
    norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)
    
    line1 = ax1.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, norm=norm, rasterized=True)
    ax1.set_xlabel('x/r',fontsize=15)                
    ax1.set_ylabel('y/r',fontsize=15)
    #ax1.set_title('Bead Orientations')
    ax1.set_xlim((downlim,uplim))
    ax1.set_ylim((downlim,uplim))
    ax1.xaxis.set_ticks( np.linspace(0,uplim,endpoint=True,num=tick_num) )
    ax1.yaxis.set_ticks( np.linspace(0,uplim,endpoint=True,num=tick_num) )
#    ax1.set_xticklabels(['0', '472', '943', '1414'])    
#    ax1.set_yticklabels(['0', '472', '943', '1414'])
    ax1.tick_params(axis='both', which='major', labelsize=10)
    #plt.setp(ax1.get_yticklabels(),visible=False)
    
#    cax1 = plt.axes([0.5, 0.5, 0.3, 0.3], projection='polar')    
#    #cax1 = plt.axes([subp.xbeg+ax_len+0.01, subp.ybeg+ax_len/3, ax_len/4.6, ax_len/4.6], projection='polar')
#    xval = np.arange(-np.pi, np.pi, 0.01)
#    yval = np.ones_like(xval)
#    cax1.scatter(xval, yval, c=xval, s=300, cmap=plt.cm.get_cmap('hsv',quant_steps), norm=norm, linewidths=0)
#    cax1.set_yticks([])
#    cax1.set_xticks([])
#    cax1.set_title('$\\phi$',fontsize=40)
#    cax1.set_rlim([-1,1])
#    cax1.set_axis_off()
     
    if args.save: 
        plt.savefig(savepath+'/frame-'+'{0:05d}'.format(timestep)+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
        plt.savefig(savepath+'/frame-'+'{0:05d}'.format(timestep)+'.eps',dpi=200,bbox_inches='tight',pad_inches=0.08)
    else:
        plt.savefig('/usr/users/iff_th2/duman/Desktop/frame-test-'+'{0:05d}'.format(timestep)+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
    
    plt.clf()  
    
    return


################################################################################

if __name__ == '__main__':
    main()        
    
    
    

