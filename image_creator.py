
# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os.path
import pandas as pd
from string import atof
import matplotlib.cm as cm
#import vapory

#########################################################################

#
# Function definitions
#

#########################################################################

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
    
#########################################################################    
    
def formatter(density, kappa, fp, datafile):
    """ Format the data into files and folders"""
    
    ## Format the data from the folder names   
    folders = []
    for d in density:
        for k in kappa:
            for f in fp:
                folders.append( datafile + '/density_' + str(d) + \
                '/kappa_' + str(k) + '/fp_' + str(f) )

    ## Format the data            
    files = []
    tmp = {}
    for folder in folders:
        
        tmp = {}
        dname = folder.split('/')[-3]
        ddict = dict(zip(dname.split('_')[::2],map(atof,dname.split('_')[1::2])))
        
        kname = folder.split('/')[-2]
        kdict = dict(zip(kname.split('_')[::2],map(atof,kname.split('_')[1::2])))
        
        fname = folder.split('/')[-1]
        fdict = dict(zip(fname.split('_')[::2],map(atof,fname.split('_')[1::2])))
    
        tmp.update(ddict)
        tmp.update(kdict)
        tmp.update(fdict)
    
        
        tmp.update({'folder':folder})
        files.append(tmp)
        
    Data = pd.DataFrame(files,columns=files[0].keys())
    
    return Data

#########################################################################
    
def get_data(path):
    """ Load the particle trajectory information"""
    
    if os.path.exists(path) == False:
        return 0
    contents = []
    for content in os.listdir(path):
        if content[:5] == 'beads':
            contents.append(int(content.split('_')[-1].split('.')[0]))
    
    ## Get the identity of the last image
    if np.size(contents) == 0:
        return 0   
    last = max(contents)
    first = min(contents)
    print first, '\t', last, '\n'
    
    last_imag_path = path + "/beads_" + str(last) + ".txt"
    #last_imag_path = path + "/beads_99950000.txt"
    
    p_last = Particles(last_imag_path)  
    
    return p_last
    
#########################################################################

def plot_with_mpl(Data, savef):
    """ Plot using matplotlib"""

    for folder in Data['folder']:
        
        ## Load preliminary information about the simulation
        print folder
        loadSimData(folder+"/init_info.txt")
        
        ## Load the data files for the last image
        path = folder + "/CLUSTER"
        p_last = get_data(path)
        if p_last == 0:
            continue
        
        ## Load plot properties
        downlim = -5
        uplim = max(Lx,Ly)+5
        tick_interval = int(uplim/5)            
        quant_steps = 2056
        norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)
        
        ## Plot
        plt.scatter(p_last.xi, p_last.yi, s=1, c=p_last.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, norm=norm)
        #ax.set_xlabel('x/b',fontsize=8)
        #ax.set_ylabel('y/b',fontsize=8)
        plt.xlim((downlim,uplim))
        plt.ylim((downlim,uplim))
        plt.xticks( np.arange(0,uplim,tick_interval) )
        plt.yticks( np.arange(0,uplim,tick_interval) )
        plt.tick_params(axis='both', which='major', labelsize=8)
        
        folder_detail = folder.split('/')
        save_path = savef + "/" + folder_detail[-3] + "/" + folder_detail[-2] + "/" + folder_detail[-1] + "/img"
        if os.path.exists(save_path):
            plt.savefig(save_path+'/fin.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
        else:
            os.mkdir(save_path)
            plt.savefig(save_path+'/fin.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
        plt.clf()    

#########################################################################

def gen_img_settings_quality(l):
    """ generate the general povray settings"""
    lhalf = 0.5*l
    # sphere radius
    sphere_radius=0.7
    #sphere_rgbcolor = [0.25,0.65,0.65]
    # RESOLUTION
    img_widthpx=1024
    img_heightpx=1024

    # includes and defaults
    povray_includes = ["colors.inc", "textures.inc", "shapes.inc"]
    povray_defaults = [vapory.Finish( 'ambient', 0.1,
	     			  'diffuse', 0.65,
		    		  'specular', 0.5,
			    	  'shininess', 0.53,
                         'opacity', 1.0)]


    # lightsources
    sun1 = vapory.LightSource([lhalf, lhalf, -1.01*lhalf], 'color', 'White')
    sun2 = vapory.LightSource([lhalf, lhalf, -1.01*lhalf], 'color', [0.7, 0.7, 0.7])

    # background
    background = vapory.Background('color', [1,1,1])

    # camera
    #povray_cam = vapory.Camera('angle', 75,'location',  [-15 , 15.0+0.5,15.0-0.25],'look_at', [0.25 , 15.0+0.5, 15.0-0.25])
    povray_cam = vapory.Camera('location', [lhalf, lhalf, -1.01*lhalf], 'look_at', [lhalf,lhalf,0], 'angle', 90)

    # text
    # If desired include this in the povray_objects - array declared in the loop
    #text1 = vapory.Text( 'ttf', '"timrom.ttf"' ,'"Division:"', 0.01, 0.0, 'scale', [0.5,0.5,0.5],'rotate', [0,90,0], 'translate' , [0.0 , 15.0+2.75-1 , 15.0+1.5], vapory.Pigment('Black') ) 

    # render quality
    quality = 10
    
    return sphere_radius, img_widthpx, img_heightpx, povray_includes, povray_defaults, sun1, sun2, background, povray_cam, quality

#########################################################################

def gen_img_settings_efficiency(l):
    """ generate the general povray settings"""
    lhalf = 0.5*l
    # sphere radius
    sphere_radius=0.7
    #sphere_rgbcolor = [0.25,0.65,0.65]
    # RESOLUTION
    img_widthpx=512
    img_heightpx=512

    # includes and defaults
    povray_includes = ["colors.inc", "textures.inc", "shapes.inc"]
    povray_defaults = [vapory.Finish( 'ambient', 0.0,
	     			  'diffuse', 0.65,
		    		  'specular', 0.0,
			    	  'shininess', 0.0,
                         'reflection', 0.0,
                         'refraction', 0.0,
				  'opacity', 1.0)]


    # lightsources
    sun1 = vapory.LightSource([lhalf, lhalf, -1.01*lhalf], 'color', 'White')
    sun2 = vapory.LightSource([lhalf, lhalf, -1.01*lhalf], 'color', [0.7, 0.7, 0.7])

    # background
    background = vapory.Background('color', [1,1,1])

    # camera
    #povray_cam = vapory.Camera('angle', 75,'location',  [-15 , 15.0+0.5,15.0-0.25],'look_at', [0.25 , 15.0+0.5, 15.0-0.25])
    povray_cam = vapory.Camera('location', [lhalf, lhalf, -1.01*lhalf], 'look_at', [lhalf,lhalf,0], 'angle', 90)

    # text
    # If desired include this in the povray_objects - array declared in the loop
    #text1 = vapory.Text( 'ttf', '"timrom.ttf"' ,'"Division:"', 0.01, 0.0, 'scale', [0.5,0.5,0.5],'rotate', [0,90,0], 'translate' , [0.0 , 15.0+2.75-1 , 15.0+1.5], vapory.Pigment('Black') ) 

    quality = 2
    
    return sphere_radius, img_widthpx, img_heightpx, povray_includes, povray_defaults, sun1, sun2, background, povray_cam, quality

#########################################################################

def gen_colors(natoms):
    """ generate an array with colors for the spheres"""
    # use colormap (blue, red, green)
    #  blue 0.0, 0.0, 1.0
    #  red: 1.0, 0.0, 0.0
    #  green: 0.0, 1.0, 1.0
    my_cmap = cm.get_cmap('gnuplot')
    minval = 0
    maxval = natoms - 1
    norm = mpl.colors.Normalize(minval, maxval)
    
    sphere_rgbcolor = []
    for i in range(natoms):
        si = my_cmap(norm(i))
        sphere_rgbcolor.append([si[0], si[1], si[2]])

    return sphere_rgbcolor

#########################################################################
    
def plot_with_povray(Data, savef, frac):
    """ Plot with povray"""

    for folder in Data['folder']:
        
        ## Load preliminary information about the simulation
        print folder
        loadSimData(folder+"/init_info.txt")
        
        ## Load the data files for the last image
        path = folder + "/CLUSTER"
        p_last = get_data(path)
        if p_last == 0:
            continue
        
        ## define the color for the spheres
        print '     defining colors'
        sphere_rgbcolor = gen_colors(L)
        z = np.zeros((L))
        
        ## create povray settings
        print '      creating povray settings'
        sphere_radius, img_widthpx, img_heightpx, povray_includes, povray_defaults, sun1, sun2, background, povray_cam, quality = gen_img_settings_quality(frac*Lx)
            
        ## create povray items
        print 'generating povray item'
        particles = vapory.Object( vapory.Union( *[ vapory.Sphere([p_last.xi[j],p_last.yi[j],z[j]], sphere_radius, vapory.Texture( vapory.Pigment('color', sphere_rgbcolor[j]), vapory.Finish('phong',1)) ) for j in range(0, len(p_last.xi) ) ] ) )
    
        ## generate povray objects
        print '      generating povray objects'
        povray_objects = [sun1, sun2, background, particles]
        
        ### create the scene
        scene = vapory.Scene( camera = povray_cam,
                       objects = povray_objects, 
                       included = povray_includes, 
                       defaults = povray_defaults )
        
        ## save path details
        folder_detail = folder.split('/')
        save_path = savef + "/" + folder_detail[-3] + "/" + folder_detail[-2] + "/" + folder_detail[-1] + "/img"
        if os.path.exists(save_path) == False:
            os.mkdir(save_path)
               
        ## render image
        print 'rendering scene'
        scene.render(outfile=save_path+'/fin_pov.png', width=img_widthpx, height=img_heightpx, antialiasing=0.001, quality=quality,remove_temp=True,show_window=False)


#########################################################################

#
# Class definitions
#

#########################################################################

class Particles:
    """ Data structure for particles"""
    
    def __init__(self, path):
        file = np.transpose(np.loadtxt(path, dtype=float))
        self.xi = file[0]/B                 # Image particle positions in x
        self.yi = file[1]/B                 # Image particle positions in y 
        self.phi = file[2]                  # Bead orientation 
        self.cidx = file[3]                 # Cluster index        

#########################################################################    

def main(): 
    
    # Argument parsing (command line options)
    parser = argparse.ArgumentParser()
    parser.add_argument("datafile", help="Folder containing data")
    parser.add_argument("savefile", help="Folder in which data should be saved")
    parser.add_argument("frac", help="Fraction of the image to show")
    args = parser.parse_args()
    
    ## Index the data
    #density = [0.08, 0.2, 0.4, 0.6, 0.8]
    density = [0.08, 0.2, 0.4]
    kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 400.0]
    fp = [0.0, 0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
    ## Format the data 
    Data = formatter(density, kappa, fp, args.datafile)
    
    ## Plot the data in Matplotlib
    plot_with_mpl(Data, args.savefile)
    
    ## Plot with Povray
    #plot_with_povray(Data, args.savefile, args.frac)




if __name__ == '__main__':
    main()



