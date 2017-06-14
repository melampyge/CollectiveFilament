
#########################################################################

import vapory
import argparse
import numpy as np
import math
import h5py
import os
import matplotlib
import matplotlib.cm as cm
import matplotlib.colors as mplcolors

#########################################################################

def loadSimData(datafile):
    """ load initial simulation data"""
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
        dtSamp, T, box_area, nt, body_length, Pe, persistence, flexure, tau_D, tau_A, \
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
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    flexure = Pe/persistence
    T = totalStep - ti
    nt = T/nsamp
    nsamp = int(nsamp)
    nbeads = int(L)
    nfil = int(N)
    nmol = int(M)
    tau_D = body_length**2*(N+1)/4/kT
    tau_A = (N+1)/Fmc
    print "Diffusive time scale is ", tau_D
    print "Advective time scale is ", tau_A
    
    return

#########################################################################
    
def gen_img_settings_quality(l):
    """ generate the general povray settings"""
    
    lhalf = 0.5*l
    
    # sphere radius
    sphere_radius = 0.7
    #sphere_rgbcolor = [0.25,0.65,0.65]
    
    # RESOLUTION
    img_widthpx = 1024
    img_heightpx = 1024

    # includes and defaults
    povray_includes = ["colors.inc", "textures.inc", "shapes.inc"]
    povray_defaults = [vapory.Finish( 'ambient', 0.1,
	     			  'diffuse', 0.65,
		    		  'specular', 0.5,
			    	  'shininess', 0.53,
				  'opacity', 1.0)]


    # light sources
    sun1 = vapory.LightSource([lhalf, lhalf, -1.01*lhalf], 'color', 'White')
    sun2 = vapory.LightSource([lhalf, lhalf, -1.01*lhalf], 'color', [0.7, 0.7, 0.7])

    # background
    background = vapory.Background('color', [1,1,1])

    # camera
    #povray_cam = vapory.Camera('angle', 75, 'location',  [-15 , 15.0+0.5,15.0-0.25],'look_at', [0.25 , 15.0+0.5, 15.0-0.25])
    povray_cam = vapory.Camera('location', [lhalf, lhalf, -1.01*lhalf], 'look_at', [lhalf,lhalf,0], 'angle', 90)

    # text
    # If desired include this in the povray_objects - array declared in the loop
    #text1 = vapory.Text( 'ttf', '"timrom.ttf"' ,'"Division:"', 0.01, 0.0, 'scale', [0.5,0.5,0.5],'rotate', [0,90,0], 'translate' , [0.0 , 15.0+2.75-1 , 15.0+1.5], vapory.Pigment('Black') ) 

    # render quality
    quality = 10
    
    return sphere_radius, img_widthpx, img_heightpx, povray_includes, povray_defaults, sun1, sun2, background, povray_cam, quality

#########################################################################

def gen_colors(nbeads):
    """ generate an array with colors for the spheres"""

    # use colormap (blue, red, green)
    #  blue 0.0, 0.0, 1.0
    #  red: 1.0, 0.0, 0.0
    #  green: 0.0, 1.0, 1.0

    my_cmap = cm.get_cmap('gnuplot')
    minval = 0
    maxval = nbeads - 1
    norm = matplotlib.colors.Normalize(minval, maxval)
    
    sphere_rgbcolor = []
    for i in range(nbeads):
        si = my_cmap(norm(i))
        sphere_rgbcolor.append([si[0], si[1], si[2]])

    return sphere_rgbcolor

#########################################################################

def gen_video(database, savebase, ti, tf):
    """ generate the video"""

    ### define the color for the spheres

    print 'defining colors'
    sphere_rgbcolor = gen_colors(nbeads)

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
        savepath = savebase + "pov-frame-" + "{0:05d}".format(int(frame)) + ".png"
        savename = "pov-frame-" + "{0:05d}".format(int(frame)) + ".png"
        
        p = Particles(path)
        p.assign_index()
      
        ### create povray items
        
        print 'generating povray item'
        particles = vapory.Object( \
            vapory.Union( \
                *[ vapory.Sphere([p.xi[j],p.yi[j],zi[j]], \
                    sphere_radius, vapory.Texture( \
                        vapory.Pigment('color', sphere_rgbcolor[j]), \
                            vapory.Finish('phong',1)) ) for j in range(0, nbeads ) ] ) )

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

class Particles:
    """ data structure for storing particle information"""
    
    def __init__(self, path):
        f = h5py.File(path, 'r')
        beads = f['beads']
        self.xi = np.asarray(beads['x'], dtype=float)/B                   # Image particle positions in x
        self.yi = np.asarray(beads['y'], dtype=float)/B                   # Image particle positions in y 
        self.phi = np.asarray(beads['phi'], dtype=float)                # Bead orientation 
        self.cidx = np.asarray(beads['idx'], dtype=int)                 # Cluster index
        f.close()
        
        return
        
    ###
        
    def assign_index(self):
        nbeads = len(self.xi)
        self.idx = np.zeros((nbeads))
        self.color = []
        nfil = 201
        for j in range(nbeads):
            molidx = int(j/nfil)
            self.idx[j] = molidx
            self.color.append('grey')
            if molidx == 23:
                self.color.append('red')
            
        return
        
#########################################################################

def main():
    
    ## do argument parsing
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-ti","--init_time", nargs="?", const=10000000, type=int, help="First time frame of the video, in timesteps, you can also leave it empty")
    parser.add_argument("-tf","--fin_time", nargs="?", const=100000000, type=int, help="Last time frame of the video, in timesteps, you can also leave it empty")
    args = parser.parse_args()
    
    ## time frames
    
    ti = args.init_time
    tf = args.fin_time
    
    ## load saved preliminary simulation data into relevant variables
    
    folderbase = "/local/duman/SIMULATIONS/"     # the data is in iff416 (assumption!)
    datafolder = folderbase + "long_filaments/density_0.2/kappa_5.0/fp_1.0/"
    initfilepath = datafolder + "init_info.txt"
    assert os.path.exists(initfilepath), "init_info.txt does NOT exist for the following file " + initfilepath
    loadSimData(initfilepath)
        
    ## set data folder paths 
        
    clusterpath = datafolder + "CLUSTER/"
    assert os.path.exists(clusterpath), "Cluster analysis has NOT been conducted for : " + clusterpath

    ## set save folder paths 

    savebase = "/usr/users/iff_th2/duman/RolfData/MOVIES/long_filaments/"
    os.system("mkdir -p " + savebase)
    savepath = savebase + "/density_0.2_kappa_5.0_fp_1.0"
    os.system("mkdir -p " + savepath)
    savepath = savepath + '/POV/'    

    gen_video(clusterpath, savepath, ti, tf)

#########################################################################

if __name__ == '__main__':
    main()

#########################################################################
