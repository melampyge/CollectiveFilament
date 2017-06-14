
##############################################################################

## 
## Generate a video of filaments to depict properties of clusters
##

##############################################################################

## load needed libraries and necessary files

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import h5py
import os
from matplotlib.patches import Rectangle
import matplotlib.cm as cm
import matplotlib.colors as mplcolors
import vapory

##############################################################################

def loadSimData(datafile):
    """ load initial simulation data"""
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
        dtSamp, T, box_area, nt, body_length, pe, xil, flexure, tau_D, tau_A, \
            nbeads, nmol, nfil

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
    pe = Fmc*body_length**2/kT
    xil = Kbend/(kT*body_length)
    flexure = pe/xil
    T = totalStep - ti
    nt = T/nsamp
    nsamp = int(nsamp)
    tau_D = body_length**2*(N+1)/4/kT
    if Fmc != 0:
        tau_A = (N+1)/Fmc
    else:
        tau_A = 0.000001
    nmol = int(M)
    nfil = int(N)
    nbeads = int(L)

    print '\n\n*** SIMULATION PARAMETERS FOR ', datafile
    print 'dt = ', dt
    print 'L = ', body_length
    print 'Pe = ', pe
    print 'xi_p/L = ', xil
    print 'T = ', T
    print 'nfil = ', nfil
    print 'nmol = ', nmol
    print 'nbeads = ', nbeads
    print 'lx = ly = ', Lx
    print "Diffusive time scale is ", tau_D
    print "Advective time scale is ", tau_A
    print '***\n\n'
    
    return

##############################################################################

def gen_img_settings_quality(xc, yc):
    """ generate the general povray settings"""
        
    back_out_factor = 100  
    sphere_radius = 0.7
    #sphere_rgbcolor = [0.25,0.65,0.65]
    
    # RESOLUTION
    img_widthpx = 4096
    img_heightpx = 4096

    # includes and defaults
    povray_includes = ["colors.inc", "textures.inc", "shapes.inc"]
    povray_defaults = [vapory.Finish( 'ambient', 0.1, 
	     			  'diffuse', 0.65, 
		    		  'specular', 0.5, 
			    	  'shininess', 0.53, 
				  'opacity', 1.0)]


    # light sources
    sun1 = vapory.LightSource([xc, yc, -1.01*back_out_factor], 'color', 'White')
    sun2 = vapory.LightSource([xc, yc, -1.01*back_out_factor], 'color', [0.7, 0.7, 0.7])

    # background
    background = vapory.Background('color', [1,1,1])

    # camera
    #povray_cam = vapory.Camera('angle', 75, 'location',  [-15 , 15.0+0.5,15.0-0.25],'look_at', [0.25 , 15.0+0.5, 15.0-0.25])
    povray_cam = vapory.Camera('location', [xc, yc, -1.01*back_out_factor], \
        'look_at', [xc, yc, 0], 'angle', 90)

    # text
    # If desired include this in the povray_objects - array declared in the loop
    #text1 = vapory.Text( 'ttf', '"timrom.ttf"' ,'"Division:"', 0.01, 0.0, 'scale', [0.5,0.5,0.5],'rotate', [0,90,0], 'translate' , [0.0 , 15.0+2.75-1 , 15.0+1.5], vapory.Pigment('Black') ) 

    # render quality
    quality = 10
    
    return sphere_radius, img_widthpx, img_heightpx, povray_includes, povray_defaults, sun1, sun2, background, povray_cam, quality

##############################################################################

def gen_colors(phi, nbeads):
    """ generate an array with colors for the spheres"""

    # use colormap (blue, red, green)
    #  blue 0.0, 0.0, 1.0
    #  red: 1.0, 0.0, 0.0
    #  green: 0.0, 1.0, 1.0

    quant_steps = 2056
    my_cmap = cm.get_cmap('hsv',quant_steps)
    minval = -np.pi
    maxval = np.pi
    norm = mpl.colors.Normalize(minval, maxval)

    sphere_rgbcolor = []
    for i in range(nbeads):
        si = my_cmap(norm(phi[i]))
        sphere_rgbcolor.append([si[0], si[1], si[2]])

    return sphere_rgbcolor

##############################################################################

def gen_video(database, savebase, ti, tf):
    """ generate the video"""

    ### create povray settings
    print 'creating povray settings'
    sphere_radius, img_widthpx, img_heightpx, povray_includes, \
        povray_defaults, sun1, sun2, background, povray_cam, quality \
            = gen_img_settings_quality(Lx)
    
    zi = np.zeros((nbeads))
        
    for frame in np.arange(ti, tf, nsamp):
        
        print frame, " of ", tf, " with ", nsamp, " intervals" 
        time = frame*dt
        
        ## data
        
        path = database + "cluster_" + str(frame) + ".hdf5"
        #savepath = savebase + "frame-" + "{0:05d}".format(int(frame)) + ".png"
        savename = "frame-" + "{0:05d}".format(int(frame)) + ".png"
        
        p = Particles(path)

        ### define the color for the spheres
        
        print 'defining colors'
        sphere_rgbcolor = gen_colors(p.phi, nbeads)
      
        ### create povray items
        
        print 'generating povray item'
        particles = vapory.Object( \
            vapory.Union( \
                *[ vapory.Sphere([p.xi[j]/B, p.yi[j]/B, zi[j]], \
                    sphere_radius, vapory.Texture( \
                        vapory.Pigment('color', sphere_rgbcolor[j]), \
                            vapory.Finish('phong',1)) ) for j in range(0, nbeads) ] ) )

        ### generate povray objects

        print 'generating povray objects'
        povray_objects = [sun1, sun2, background, particles]
        
        ### create the scene
        
        scene = vapory.Scene( camera = povray_cam,
                       objects = povray_objects, 
                       included = povray_includes, 
                       defaults = povray_defaults )
                       
        ### render image
                           
        print 'rendering scene'
        scene.render(outfile=savename, width=img_widthpx, height=img_heightpx, \
            antialiasing=0.001, quality=quality, remove_temp=True)

    return   
    
##############################################################################
    
def generate_povray_zoom(x, y, phi, xc, yc, savename):
    """ zoom in on a specific part via povray"""
    
    natoms = len(x)
    z = np.zeros((natoms))
    
    ### define the color for the spheres
        
    print 'defining colors'
    sphere_rgbcolor = gen_colors(phi, natoms)
    
    ### create povray settings

    print 'creating povray settings'
    sphere_radius, img_widthpx, img_heightpx, povray_includes, \
        povray_defaults, sun1, sun2, background, povray_cam, quality \
            = gen_img_settings_quality(xc, yc)
  
    ### create povray items
    
    print 'generating povray item'
    particles = vapory.Object( \
        vapory.Union( \
            *[ vapory.Sphere([x[j], y[j], z[j]], \
                sphere_radius, vapory.Texture( \
                    vapory.Pigment('color', sphere_rgbcolor[j]), \
                        vapory.Finish('phong',1)) ) for j in range(0, natoms) ] ) )

    ### generate povray objects

    print 'generating povray objects'
    povray_objects = [sun1, sun2, background, particles]
    
    ### create the scene
    
    scene = vapory.Scene( camera = povray_cam,
                   objects = povray_objects, 
                   included = povray_includes, 
                   defaults = povray_defaults )
                   
    ### render image
                       
    print 'rendering scene'
    scene.render(outfile=savename, width=img_widthpx, height=img_heightpx, \
        antialiasing=0.001, quality=quality, remove_temp=True)
            
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
        self.xi = np.asarray(beads['x'], dtype=float)/B                 # Image particle positions in x
        self.yi = np.asarray(beads['y'], dtype=float)/B                 # Image particle positions in y 
        self.phi = np.asarray(beads['phi'], dtype=float)                # Bead orientation 
        self.cidx = np.asarray(beads['idx'], dtype=int)                 # Cluster index
        f.close()
        
        return
        
    ###
        
    def assign_index(self):
        """ assign filament index to beads"""
        
        nbeads = len(self.xi)
        self.idx = np.zeros((nbeads))
        nfil = 51
        for j in range(nbeads):
            self.idx[j] = int(j/nfil)
            
        return
    
    ###
        
    def mask_beads(self, xdown, ydown, xup, yup):
        """ find the beads in the given range"""
        
        xm = np.ma.masked_outside(self.xi, xdown, xup)
        ym = np.ma.masked_outside(self.yi, ydown, yup)
        
        #x = xm[~xm.mask].data
        #y = ym[~ym.mask].data
        xatomids = np.where(xm.mask == False)
        yatomids = np.where(ym.mask == False)
        atomids = np.intersect1d(xatomids[0], yatomids[0])
        
        return self.xi[atomids], self.yi[atomids], atomids
        
    ### 

    def track_bead(self, bead_idx, xtraj, ytraj):
        """ track the motion of the specified bead"""
        
        xtraj.append(self.xi[bead_idx])
        ytraj.append(self.yi[bead_idx])
        
        return

################################################################################         
        
class ClusterInfo:
    """ Data structure for cluster information"""
    
    def __init__(self, path):
        f = h5py.File(path, 'r')
        props = f['props']    
        self.idx = np.array(props['idx'], dtype=int)
        self.size = np.array(props['size'], dtype=int)
        self.comx = np.array(props['comx'], dtype=float)
        self.comy = np.array(props['comy'], dtype=float)
        self.rgysq = np.array(props['rgysq'], dtype=float)
        self.pl = np.array(props['planarity'], dtype=float)
        self.st = np.array(props['straightness'], dtype=float)
        self.sw = np.array(props['swirliness'], dtype=float)
        self.ens = np.array(props['enstrophy'], dtype=float) 
        fil_grp = props['filament_list']
        self.fils = []
        for el in np.arange(len(self.size)):
            fil_list = np.array(fil_grp[str(el)], dtype=int)
            self.fils.append(fil_list)
       
        f.close()
        
        return
        
    ###
        
    def find_index(self, xdown, ydown, xup, yup):
        """ find the index of the cluster within the specified range"""

        xm = np.ma.masked_outside(self.comx, xdown, xup)
        ym = np.ma.masked_outside(self.comy, ydown, yup)

        xids = np.where(xm.mask == False)
        yids = np.where(ym.mask == False)
        cids = np.intersect1d(xids[0], yids[0])

        print cids
        print self.idx[cids]        
        
        return 
        
    ### 
    
    def track_com(self, cidx, comx, comy):
        """ track the center of mass of a given cluster with index cidx"""
        
        array_idx = np.where(self.idx == cidx)
        j = array_idx[0]
        print "Array index = ", j
        print "Cluster index = ", self.idx[j]
        print "Cluster size = ", self.size[j]
        print "comx = ", self.comx[j]
        print "comy = ", self.comy[j], '\n\n'
        
        comx.append(self.comx[j])
        comy.append(self.comy[j])

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

def plot_frames(ti, tf, database, savebase, xd, yd, xu, yu):
    """ plot some of the time frames"""

    ## set plot properties

    ax_len = 1.0                          # Length of one subplot square box
    ax_b = 0.0                            # Beginning/offset of the subplot in the box
    ax_sep = 0.0                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots    
    fig = plt.figure()
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
    
    ## set more plot properties

    downlim = 0
    uplim = max(Lx,Ly)
    quant_steps = 2056
    norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)  
    num_ticks = 5
    markersize = 9.0
    
    ## plot the frames
    
    comx = []
    comy = []
    beadx = []
    beady = []
    
    for frame in np.arange(ti, tf+nsamp, nsamp):
        
        print frame, " of ", tf, " with ", nsamp, " intervals" 
        time = frame*dt
        
        ## data
        
        path = database + "cluster_" + str(frame) + ".hdf5"
        savepath1 = savebase + "frame-" + "{0:05d}".format(int(frame)) + ".png"
        savepath2 = savebase + "frame-" + "{0:05d}".format(int(frame)) + ".eps"
        savename = "povray_zoom_" + "{0:05d}".format(int(frame)) + ".png"
        
        p = Particles(path)
        cl = ClusterInfo(path) 
        #cl.find_index(xd, yd, xu, yu)
        ##### clidx = 2624
        cl.track_com(2624, comx, comy)
        #p.mask_beads(301, 555, 305, 605)
        p.track_bead(6911, beadx, beady)
        
        if frame == tf or frame == ti or frame == 23000000 or frame == 22100000:
            
            print np.shape(comx)
        
            comx_arr = np.array(comx) + Lx
            beadx_arr = np.array(beadx) + Lx
            
            ## plot
            
            subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
            ax0 = subp.addSubplot()
                  
            line0 = ax0.scatter(p.xi, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                        edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, norm=norm, rasterized=True)
            line0 = ax0.scatter(p.xi+Lx, p.yi, s=1, c=p.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                        edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, norm=norm, rasterized=True)
            
            ax0.axis('scaled')
            
            if frame == ti:
                ax0.plot(comx_arr, comy, 'o', color='red', alpha=1.0, markersize=markersize)
#                xprobe = p.xi[6886:6936]
#                yprobe = p.yi[6886:6936]
                ax0.plot(p.xi[6911]+Lx, p.yi[6911], 'o', color='yellow', alpha=1.0, markersize=markersize)
                #generate_povray_zoom(p.xi, p.yi, p.phi, 302, 570, savename)
                ax0.plot(comx, comy, 'o', color='red', alpha=1.0, markersize=markersize)                
                ax0.plot(p.xi[6911], p.yi[6911], 'o', color='yellow', alpha=1.0, markersize=markersize)

                
            elif frame == 22100000:
                ax0.plot(comx_arr[0], comy[0], 'o', color='red', alpha=1.0, markersize=markersize)                
                ax0.plot(beadx_arr[0], beady[0], 'o', color='yellow', alpha=1.0, markersize=markersize)                
                ax0.plot(comx_arr, comy, color='red', alpha=1.0, linewidth=2.0)
                ax0.plot(beadx_arr, beady, color='yellow', alpha=1.0, linewidth=2.0)
                
                ax0.plot(comx[0], comy[0], 'o', color='red', alpha=1.0, markersize=markersize)                
                ax0.plot(beadx[0], beady[0], 'o', color='yellow', alpha=1.0, markersize=markersize)                
                ax0.plot(comx, comy, color='red', alpha=1.0, linewidth=2.0)
                ax0.plot(beadx, beady, color='yellow', alpha=1.0, linewidth=2.0)
                
            elif frame == 23000000:
                ax0.plot(comx_arr[0], comy[0], 'o', color='red', alpha=1.0, markersize=markersize)                
                ax0.plot(beadx_arr[0], beady[0], 'o', color='yellow', alpha=1.0, markersize=markersize)                
                ax0.plot(comx_arr, comy, color='red', alpha=1.0, linewidth=2.0)
                ax0.plot(beadx_arr, beady, color='yellow', alpha=1.0, linewidth=2.0)
 
                ax0.plot(comx[0], comy[0], 'o', color='red', alpha=1.0, markersize=markersize)                
                ax0.plot(beadx[0], beady[0], 'o', color='yellow', alpha=1.0, markersize=markersize)                
                ax0.plot(comx, comy, color='red', alpha=1.0, linewidth=2.0)
                ax0.plot(beadx, beady, color='yellow', alpha=1.0, linewidth=2.0)
                
            else:
                ax0.plot(comx_arr[0], comy[0], 'o', color='red', alpha=1.0, markersize=markersize)                
                ax0.plot(beadx_arr[0], beady[0], 'o', color='yellow', alpha=1.0, markersize=markersize)                 
                ax0.plot(comx_arr, comy, color='red', alpha=1.0, linewidth=2.0)
                ax0.plot(beadx_arr, beady, color='yellow', alpha=1.0, linewidth=2.0)

                ax0.plot(comx[0], comy[0], 'o', color='red', alpha=1.0, markersize=markersize)                
                ax0.plot(beadx[0], beady[0], 'o', color='yellow', alpha=1.0, markersize=markersize)                 
                ax0.plot(comx, comy, color='red', alpha=1.0, linewidth=2.0)
                ax0.plot(beadx, beady, color='yellow', alpha=1.0, linewidth=2.0)
                
            ## title
            
            ax0.set_title("$t/\\tau_{D}$ = " + "{0:.2f}".format(time/tau_D) + \
                ", $t/\\tau_{A}$ = " + "{0:.2f}".format(time/tau_A), fontsize=30)
            
            ## labels 
                
            ax0.set_xlabel("$x/r_{0}$", fontsize=30)
            ax0.set_ylabel("$y/r_{0}$", fontsize=30)
    
            ## limits
    
            ax0.set_xlim((600, 600+Lx))
            ax0.set_ylim((downlim, uplim))
            
            ## ticks
            
            #ax0.xaxis.set_ticks(np.linspace(600, 600+Lx, num_ticks, endpoint=True))
            #ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
            ax0.xaxis.set_ticks([600, 1100, 1600, 2100])
            ax0.yaxis.set_ticks([0, 500, 1000, 1500])
            ax0.set_xticklabels([0, 500, 1000, 1500])
            ax0.tick_params(axis='both', which='major', labelsize=30)
    
            #ax0.add_patch( Rectangle( (301,555),4,50,edgecolor='gray',facecolor='gray',alpha=0.7 ) )
            
            ## colorbar
            
#            cax0 = plt.axes([subp.xbeg+ax_len+0.01, subp.ybeg+ax_len/3, ax_len/4.6, ax_len/4.6], projection='polar')
#            xval = np.arange(-np.pi, np.pi, 0.01)
#            yval = np.ones_like(xval)
#            cax0.scatter(xval, yval, c=xval, s=300, cmap=plt.cm.get_cmap('hsv',quant_steps), norm=norm, linewidths=0)
#            cax0.set_xticks([])
#            cax0.set_yticks([])
#            cax0.set_title('$\\phi$',fontsize=20)
#            cax0.set_rlim([-1,1])
#            cax0.set_axis_off()
                        
            ## save

            plt.savefig(savepath1, dpi=200, bbox_inches='tight', pad_inches=0.08)
            plt.savefig(savepath2, dpi=200, bbox_inches='tight', pad_inches=0.08)        
            fig.clf()
            
            if frame == tf:
                return
            


    return

##############################################################################        

def main():
    
    ## time frames
    
    #ti = 17300000
    #ti = 16950000
    #tf = 19200000
    ti = 21000000
    tf = 24500000
    xdown = 434.0
    ydown = 0.0
    xup = 1700.0
    yup = 434.0
    
    ## load saved preliminary simulation data into relevant variables
    
    datafolder = "/local/duman/SIMULATIONS/many_polymers_5/density_0.2/kappa_400.0/fp_0.08/"
    initfilepath = datafolder + "init_info.txt"
    assert os.path.exists(initfilepath), "init_info.txt does NOT exist for the following file " + initfilepath
    loadSimData(initfilepath)
        
    ## set data folder paths 
        
    clusterpath = datafolder + "CLUSTER/"
    assert os.path.exists(clusterpath), "Cluster analysis has NOT been conducted for : " + clusterpath

    ## set save folder paths 

    savebase = "/usr/users/iff_th2/duman/RolfData/MOVIES/many_polymers_5/"
    os.system("mkdir -p " + savebase)
    savepath = savebase + "density_0.2_kappa_400.0_fp_0.08/"
    os.system("mkdir -p " + savepath)
    savepath = savepath + "ANALYSIS/"
    os.system("mkdir -p " + savepath)
    
    
    ## plot certain time frames to generate a movie later on
    
    print "Generating images for the following file : " + datafolder
    print "From frame : " + str(ti) + " to frame : " + str(tf) + " with nsamp : " + str(nsamp)  
    plot_frames(ti, tf, clusterpath, savepath, xdown, ydown, xup, yup)
    
    return
        
##############################################################################
        
if __name__ == "__main__":
    main()
    
##############################################################################    
    

