#!/usr/users/iff_th2/isele/Applications/Anaconda/bin/python2.7


import vapory

import numpy as np
import math
import sys
import os
import BackwardsReader

try:
    ifname = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '     dump file'
    exit()

#########################################################################
### General parameters

def define_povray_settings(lx):
    lhalf = lx/2
    # SPHERES TO USE FOR PARTICLES
    sphere_radius=1.0
    #sphere_rgbcolor = [0.25,0.65,0.65]

    # RESOLUTION
    img_widthpx=2048
    img_heightpx=2048


    #INCLUDES & DEFAULTS
    povray_includes = ["colors.inc", "textures.inc", "shapes.inc"]
    povray_defaults = [vapory.Finish( 'ambient', 0.1,
	    			  'diffuse', 0.65,
		    		  'specular', 0.5,
			    	  'shininess', 0.53,
				      'opacity', 1.0)]


    # LIGHTSOURCE
    sun1 = vapory.LightSource([lhalf, lhalf, -1.01*lhalf], 'color', 'White')
    sun2 = vapory.LightSource([lhalf, lhalf, -1.01*lhalf], 'color', [0.7, 0.7, 0.7])


    # BACKGROUND
    background = vapory.Background('color', [1,1,1])

    # CAMERA
    #povray_cam = vapory.Camera('angle', 75,'location',  [-15 , 15.0+0.5,15.0-0.25],'look_at', [0.25 , 15.0+0.5, 15.0-0.25])
    povray_cam = vapory.Camera('location', [lhalf, lhalf, -1.01*lhalf], 'look_at', [lhalf,lhalf,0], 'angle', 90)



    # TEXT
    # If desired include this in the povray_objects - array declared in the loop
    text1 = vapory.Text( 'ttf', '"timrom.ttf"' ,'"Division:"', 0.01, 0.0, 'scale', [0.5,0.5,0.5],'rotate', [0,90,0], 'translate' , [0.0 , 15.0+2.75-1 , 15.0+1.5], vapory.Pigment('Black') ) 
    return sphere_radius, img_widthpx, img_heightpx, povray_includes, povray_defaults, sun1, sun2, background, povray_cam, text1
    

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

def gen_image():
    """ generate the video"""
    print 'Reading Data'
    ### generate temporary file for readine
    ifile = open(ifname, 'r')
    ### get information from the header file
    # box dimensions in x and y plane & number of atoms
    for i in range(3):
        ifile.readline()
    line = ifile.readline()
    line = line.split()
    natoms = int(line[0])
    ifile.readline()
    line = ifile.readline()
    line = line.split()
    xlo = float(line[0])
    xhi = float(line[1])
    lx = xhi - xlo
    line = ifile.readline()
    line = line.split()
    ylo = float(line[0])
    yhi = float(line[1])
    ly = yhi - ylo
    ifile.close()
    ### define scene objects
    sphere_radius, img_widthpx, img_heightpx, povray_includes, povray_defaults, sun1, sun2, background, povray_cam, text1 = define_povray_settings(lx)
    ### define the color for the spheres
    sphere_rgbcolor = gen_colors(natoms)
    ### allocate arrays to store the data
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    z = np.zeros((natoms))
    # read infile backwards
    ifile = BackwardsReader.BackwardsReader(open(ifname))
    for i in range(natoms):
        line = ifile.readline()
        line = line.split()
        aID = int(line[0]) - 1
        xs = float(line[2])
        ys = float(line[3])
        xs = xs - math.floor(xs)
        ys = ys - math.floor(ys)
        x[aID] = xs*lx
        y[aID] = ys*ly

    print 'Generating Povray'
    # create povray items
    particles = vapory.Object( vapory.Union( *[ vapory.Sphere([x[j],y[j],z[j]], sphere_radius, vapory.Texture( vapory.Pigment('color', sphere_rgbcolor[j]), vapory.Finish('phong',1)) ) for j in range(0, len(x) ) ] ) )

    # generate povray objects
    povray_objects = [sun1, sun2, background, particles]
    # create the scene
    scene = vapory.Scene( camera = povray_cam,
                   objects = povray_objects, 
                   included = povray_includes, 
                   defaults = povray_defaults )

    print 'Rendering'
    # render image
    cur_name = "final.png"
    scene.render(outfile=cur_name, width=img_widthpx, height=img_heightpx, antialiasing=0.001, quality=10,remove_temp=True,show_window=False)

    return


#########################################################################

def main():
    """ main function, called when the script is started"""
    # generate an image
    gen_image()

#########################################################################

if __name__ == '__main__':
    main()

