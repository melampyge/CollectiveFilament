#!/usr/users/iff_th2/isele/Applications/Anaconda/bin/python2.7

import vapory
#from subprocess import check_output
#from itertools import islice

import numpy as np
import math
import sys
import codecs
import os
import read_char

try:
    ifname = sys.argv[1]
    hname = sys.argv[2]
    ofname = sys.argv[3]
    nevery = int(sys.argv[4])
    fast = int(sys.argv[5])
except:
    print 'Usage: ' + sys.argv[0] + '     char fiel        header file          outfilename           nevery        fast (1=fast, 0=high quality)         ending for first frame (optional)'
    exit()

try:
    nzero = int(sys.argv[6])
except:
    nzero = 0

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
    # use colormap (red, white, blue)
    #  red: 1.0, 0.0, 0.0
    #  white: 1.0, 1.0, 1.0
    #  blue: 0.0, 0.0, 1.0
    sphere_rgbcolor = []
    for i in range(natoms/2):
        si = [1.0, float(2*i)/natoms, float(2*i)/natoms]
        sphere_rgbcolor.append(si)
    for i in range(natoms/2):
        si = [1.0 - float(2*i)/natoms, 1.0 - float(2*i)/natoms, 1.0]
	sphere_rgbcolor.append(si)
    sphere_rgbcolor.append(sphere_rgbcolor[-1])
    sphere_rgbcolor.append(sphere_rgbcolor[-1])
    return sphere_rgbcolor

#########################################################################

def gen_video():
    """ generate the video"""
    ### open files for reading
    hfile = open(hname)
    ifile = codecs.open(ifname, 'r', 'UTF-8')
    ### get information from the header file
    natoms, nsteps = read_char.read_first(hfile)
    ### compute the number of images
    nimg = nsteps / nevery
    ### define the color for the spheres
    sphere_rgbcolor = gen_colors(natoms)
    ### loop over all images
    for i in range(nimg):
        ### print stats
        print 'current image / all images', i, '/', nimg
        print '    reading atoms'
        ### read current step and skip nevery-1 other steps
        x,y,lx,ly,tstep,natoms = read_char.read_snapshot(hfile, ifile)
        x *= lx
        y *= ly
        read_char.skip_snapshots(hfile, ifile, nevery-1)

        
        ### on first step: generate general povray settings and z array
        if i == 0:
            print 'defining general settings'
            z = np.zeros((natoms))
            # create povray settings
            if fast == 0:
                print 'quality'
                sphere_radius, img_widthpx, img_heightpx, povray_includes, povray_defaults, sun1, sun2, background, povray_cam, quality = gen_img_settings_quality(lx)
            if fast == 1:
                print 'efficiency'
                sphere_radius, img_widthpx, img_heightpx, povray_includes, povray_defaults, sun1, sun2, background, povray_cam, quality = gen_img_settings_efficiency(lx)
        
        ### create povray items
        print 'generating povray item'
        particles = vapory.Object( vapory.Union( *[ vapory.Sphere([x[j],y[j],z[j]], sphere_radius, vapory.Texture( vapory.Pigment('color', sphere_rgbcolor[j]), vapory.Finish('phong',1)) ) for j in range(0, len(x) ) ] ) )

        ### generate povray objects
        povray_objects = [sun1, sun2, background, particles]
        ### create the scene
        scene = vapory.Scene( camera = povray_cam,
                   objects = povray_objects, 
                   included = povray_includes, 
                   defaults = povray_defaults )
        ### render image
        cur_name = ofname + "%04d"%(i+nzero) + ".png"
        print 'rendering scene'
        scene.render(outfile=cur_name, width=img_widthpx, height=img_heightpx, antialiasing=0.001, quality=quality,remove_temp=True,show_window=False)
    # close input file
    ifile.close()
    # create subfolder and move images
    os.system('mkdir IMG')
    os.system('mv img*.png IMG')

    return


#########################################################################

def main():
    """ main function, called when the script is started"""
    # generate the images
    gen_video()

#########################################################################

if __name__ == '__main__':
    main()

