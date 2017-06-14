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
import matplotlib
import matplotlib.cm as cm

try:
    ifname = sys.argv[1]
    hname = sys.argv[2]
    ofname = sys.argv[3]
    frac = float(sys.argv[4])
except:
    print 'Usage: ' + sys.argv[0] + '     char fiel        header file          outfilename        fraction to show'
    exit()

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
    norm = matplotlib.colors.Normalize(minval, maxval)
    
    sphere_rgbcolor = []
    for i in range(natoms):
        si = my_cmap(norm(i))
        sphere_rgbcolor.append([si[0], si[1], si[2]])

    return sphere_rgbcolor

#########################################################################

def gen_video():
    """ generate the video"""
    ### open files for reading
    hfile = open(hname)
    ifile = codecs.open(ifname, 'r', 'UTF-8')
    ### get information from the header file
    natoms, nsteps = read_char.read_first(hfile)
    ### skip all but the last snapshots
    print '     skipping snapshots'
    read_char.skip_snapshots(hfile, ifile, nsteps - 1)
    ### define the color for the spheres
    print '     defining colors'
    sphere_rgbcolor = gen_colors(natoms)
    ### read the last step
    print '      reading final snapshot'
    x,y,lx,ly,tstep,natoms = read_char.read_snapshot(hfile, ifile)
    x *= lx
    y *= ly

    z = np.zeros((natoms))
    # create povray settings
    print '      creating povray settings'
    sphere_radius, img_widthpx, img_heightpx, povray_includes, povray_defaults, sun1, sun2, background, povray_cam, quality = gen_img_settings_quality(frac*lx)
        
    ### create povray items
    print 'generating povray item'
    particles = vapory.Object( vapory.Union( *[ vapory.Sphere([x[j],y[j],z[j]], sphere_radius, vapory.Texture( vapory.Pigment('color', sphere_rgbcolor[j]), vapory.Finish('phong',1)) ) for j in range(0, len(x) ) ] ) )

    ### generate povray objects
    print '      generating povray objects'
    povray_objects = [sun1, sun2, background, particles]
    ### create the scene
    scene = vapory.Scene( camera = povray_cam,
                   objects = povray_objects, 
                   included = povray_includes, 
                   defaults = povray_defaults )
    ### render image
    cur_name = ofname
    print 'rendering scene'
    scene.render(outfile=cur_name, width=img_widthpx, height=img_heightpx, antialiasing=0.001, quality=quality,remove_temp=True,show_window=False)
    # close input file
    ifile.close()
    hfile.close()

    return


#########################################################################

def main():
    """ main function, called when the script is started"""
    # generate the images
    gen_video()

#########################################################################

if __name__ == '__main__':
    main()

