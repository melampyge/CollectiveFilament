#!/usr/local/bin/python2.7

import matplotlib as mpl
mpl.use('Agg')
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy import optimize


#####################################################################
### define / read in some global variables
gamma = 2.0  # viscosity
kT = 1.0     # thermal energy

ifile = open('scale_params.data')
line = ifile.readline()
line = line.split()
L = float(line[-1])    # polymer length
line = ifile.readline()
line = line.split()
dt = float(line[-1])     # simulation timestep
ttrans = gamma*L**3/6./kT

#####################################################################

def decompose_forces():
    #try:
    #    ifile = open('run2.xyz')
    #except:
    #    ifile = open('vis2.xyz')
    ifile = open('vis2.xyz')
    fext = []
    ftrans = []
    frot = []
    t = []
    # loop over the entire file
    while True:
        # check whether file end has been reached
        line = ifile.readline()
        # read in the x and y coordinates
        line = line.split()
        try:
            natoms = int(line[0])
        except:
            break
        x = np.zeros((natoms))
        y = np.zeros((natoms))
        line = ifile.readline() 
        line = line.split()
        t.append(float(line[-1]))
        for i in range(natoms):
            line = ifile.readline()
            line = line.split()
            x[i] = float(line[1])
            y[i] = float(line[2])
        # compute the center of mass
        cxi = np.average(x)
        cyi = np.average(y)
        # compute the total external force
        fe = 0
        for i in range(natoms - 1):
            x0 = x[i]
            y0 = y[i]
            x1 = x[i+1]
            y1 = y[i+1]
            dx = x0 - x1
            dy = y0 - y1
            fe += np.sqrt(dx**2 + dy**2)
        # compute the translational force
        fix = 0
        fiy = 0
        for i in range(natoms - 1):
            x0 = x[i]
            y0 = y[i]
            x1 = x[i+1]
            y1 = y[i+1]
            dx = x0 - x1
            dy = y0 - y1
            fix += dx
            fiy += dy
        ft = np.sqrt(fix**2 + fiy**2)
        # compute the rotational force
        fr = 0
        dfx_com = fix/ft  # direction of the com force
        dfy_com = fiy/ft
        for i in range(natoms - 1):
            x0 = x[i]
            y0 = y[i]
            x1 = x[i+1]
            y1 = y[i+1]
            dx = x0 - x1
            dy = y0 - y1
            # project bond force on translational force
            pr = dx*dfx_com + dy*dfy_com
            dx = dx*(1-pr)
            dy = dy*(1-pr)
            xc = 0.5*(x0 + x1)
            yc = 0.5*(y0 + y1)
            fr += ((xc-cxi)*dy - (yc-cyi)*dx)/np.sqrt((xc-cxi)**2 + (yc-cyi)**2)
        fext.append(fe)
        ftrans.append(ft)
        frot.append(fr)
    ifile.close()
    t = np.array(t)
    fext = np.array(fext)
    ftrans = np.array(ftrans)
    frot = np.array(frot)
    return t, fext, ftrans, frot


#####################################################################

def main():
    """ main function"""
    # get ftrans, frot
    t,fext, ftrans, frot = decompose_forces()
    # write the results to a file
    ofile = open('force_decomposition.data', 'w')
    ofile.write('Results of the force decomposition\n\n')
    ofile.write('t\tf_total\tf_trans\tf_rot\n')
    for i in range(len(t)):
        ofile.write(str(t[i]) + '\t' + str(fext[i]) + '\t' + str(ftrans[i]) + '\t' + str(frot[i]) + '\n')
    ofile.close()


    # create histograms

    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left+width+0.02

    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular Figure
    plt.figure(1, figsize=(8,8))

    axScatter = plt.axes(rect_scatter)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # the scatter plot:
    axScatter.hist2d(ftrans/fext, np.fabs(frot)/fext, bins = [50,50], range = [[0,1],[0,1]], cmin = 1, cmap = 'binary', norm=colors.LogNorm())

    # now determine nice limits by hand:
    binwidth = 0.02
    xymax = 1
    lim = ( int(xymax/binwidth) + 1) * binwidth

    axScatter.set_xlim( (0, lim) )
    axScatter.set_ylim( (0, lim) )

    bins = np.arange(0, lim + binwidth, binwidth)
    axHistx.hist(ftrans/fext, bins=bins)
    axHisty.hist(np.fabs(frot)/fext, bins=bins, orientation='horizontal')

    axHistx.set_xlim( axScatter.get_xlim() )
    axHisty.set_ylim( axScatter.get_ylim() )

    plt.savefig('hist_forces.png')
    plt.close()
    return

#####################################################################

if __name__ == '__main__':
    main()
