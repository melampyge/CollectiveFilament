#!/usr/local/bin/python2.7

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import read_char
import codecs
import os

try:
    ifname = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '      parameter file'
    exit()
    

#################################################################

def read_settings():
    """ read in the settings from the parameter file"""
    # open file for reading
    ifile = open(ifname, 'r')
    # skip comment line
    ifile.readline()
    # char-file name
    line = ifile.readline()
    line = line.split()
    charfile = line[-1]
    # header-file name
    line = ifile.readline()
    line = line.split()
    headerfile = line[-1]
    # output folder
    line = ifile.readline()
    line = line.split()
    ofname = line[-1]
    # length of the filaments
    line = ifile.readline()
    line = line.split()
    nfil = int(line[-1])
    # number of snapshots to skip
    line = ifile.readline()
    line = line.split()
    nskip = int(line[-1])
    # number of histogram bins
    line = ifile.readline()
    line = line.split()
    nbins = int(line[-1])
    # close file
    ifile.close()
    # return input values
    return charfile, headerfile, ofname, nfil, nskip, nbins

##################################################################

def nearest_neighbor(x1,x3,lx):
    """ compute vector to nearest neighbor"""
    dx1 = x3 - x1
    dx2 = x3 - x1 + lx
    dx3 = x3 - x1 - lx
    if dx1**2 < dx2**2 and dx1**2 < dx3**2:
        return dx1
    if dx2**2 < dx3**2:
        return dx2
    return dx3

##################################################################

def get_velocity(x1,y1,x3,y3,lx,ly,tstep1,tstep3,natoms):
    """ compute the velocity of each bead, consider pbc"""
    # define target arrays
    vx = np.zeros(natoms)
    vy = np.zeros(natoms)
    # loop over all beads
    for i in range(natoms):
        xi = x1[i]
        yi = y1[i]
        xj = x3[i]
        yj = y3[i]
        dx = nearest_neighbor(xi,xj,lx)
        dy = nearest_neighbor(yi,yj,ly)
        vx[i] = dx/(tstep3-tstep1)
        vy[i] = dy/(tstep3-tstep1)
    return vx, vy

##################################################################

def get_forces(x,y,lx,ly,natoms,nfil):
    """ compute the external force that acts onto each bead, consider pbc"""
    # define target arrays
    fx = np.zeros((natoms))
    fy = np.zeros((natoms))
    # loop over all beads
    for i in range(natoms):
        # is current atom the last atom?
        if i % nfil == nfil - 1:
            x1 = x[i-1]
            x2 = x[i]
            y1 = y[i-1]
            y2 = y[i]     
        # is current atom the first atom?
        elif i % nfil == 0:
            x1 = x[i]
            x2 = x[i+1]
            y1 = y[i]
            y2 = y[i+1]
        # else
        else:
            x1 = x[i-1]
            x2 = x[i+1]
            y1 = y[i-1]
            y2 = y[i+1]
        dx = nearest_neighbor(x2,x1,lx)
        dy = nearest_neighbor(y2,y1,ly)
        fx[i] = 0.5*dx
        fy[i] = 0.5*dy
    # return results
    return fx,fy

##################################################################

def bin_data(x,y,vx,vy,fx,fy,natoms,lx,ly,bx,by):
    """ created binned densities and velocities"""
    # allocate output arrays
    rho = np.zeros((bx,by))
    vx_b = np.zeros((bx,by))
    vy_b = np.zeros((bx,by))
    px_b = np.zeros((bx,by))
    py_b = np.zeros((bx,by))
    fx_b = np.zeros((bx,by))
    fy_b = np.zeros((bx,by))
    # fill all atoms into bins
    for i in range(natoms):
        # get coordinates
        xi = x[i]
        yi = y[i]
        # get current bin
        segx = int(xi/lx*bx)
        segy = int(yi/ly*by)
        # add data to bin
        rho[segx,segy] += 1
        px_b[segx,segy] += vx[i]
        py_b[segx,segy] += vy[i]
        fx_b[segx,segy] += fx[i]
        fy_b[segx,segy] += fy[i]
    # transform moments to velocities
    for i in range(bx):
        for j in range(by):
            if rho[i,j] > 1:
                vx_b[i,j] = px_b[i,j]/rho[i,j]
                vy_b[i,j] = py_b[i,j]/rho[i,j]
    # transform number counts to densities
    wx = lx/bx
    wy = ly/by
    rho /= wx*wy
    # generate an array that contains the edges
    xedges = np.linspace(0., 1., bx+1)*lx
    yedges = np.linspace(0., 1., by+1)*ly
    return xedges, yedges, rho, vx_b, vy_b, px_b, py_b, fx_b, fy_b

##################################################################

def compute_rotation(vx,vy,bx,by,lx,ly):
    """ compute the rotation of a vector field"""
    # allocate output array
    rotation = np.zeros((bx,by))
    # bin width
    wx = lx/bx
    wy = ly/by
    # compute rotation in each bin
    for i in range(bx):
        for j in range(by):
            # compute velocity gradients using finite differences
            duy_dx = (vy[(i+1)%bx,j] - vy[i-1,j])/(2*wx)
            dux_dy = (vx[i,(j+1)%by] - vx[i,j-1])/(2*wy)
            rotation[i,j] = duy_dx - dux_dy
    return rotation


##################################################################

def analyse_data(charfile, headerfile, nfil, nskip, nbins):
    """ analyse the simulation data"""
    ### open input file
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile)
    natoms, nsteps = read_char.read_first(hfile)
    ### skip initial snapshots and correct nsteps value
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - 2 - nskip
    ### allocate histograms and time output array
    bx = nbins
    by = nbins
    time = np.zeros((nsteps), dtype = int)
    rho_hist = np.zeros((nsteps,bx,by))
    vx_hist = np.zeros((nsteps,bx,by))
    vy_hist = np.zeros((nsteps,bx,by))
    worticity = np.zeros((nsteps,bx,by))
    px_hist = np.zeros((nsteps,bx,by))
    py_hist = np.zeros((nsteps,bx,by))
    momticity = np.zeros((nsteps,bx,by))
    fx_hist = np.zeros((nsteps,bx,by))
    fy_hist = np.zeros((nsteps,bx,by))
    forticity = np.zeros((nsteps,bx,by))
    ### read in the first two steps
    xs, ys, lx, ly, tstep2, natoms = read_char.read_snapshot(hfile, ifile)
    x2 = xs*lx
    y2 = ys*ly
    xs, ys, lx, ly, tstep3, natoms = read_char.read_snapshot(hfile, ifile)
    x3 = xs*lx
    y3 = ys*ly
    ### loop over all steps
    for i in range(nsteps):
        # print some information
        print 'Process:', i+1, '/', nsteps
        # move previous data to other arrays
        x1 = np.copy(x2)
        y1 = np.copy(y2)
        tstep1 = tstep2
        x2 = np.copy(x3)
        y2 = np.copy(y3)
        tstep2 = tstep3
        # read in the new coordinates
        xs,ys,lx,ly,tstep3,natoms = read_char.read_snapshot(hfile, ifile)
        x3 = xs*lx
        y3 = ys*ly
        # approximate the velocity of the particles (finite difference scheme)
        vx,vy = get_velocity(x1,y1,x3,y3,lx,ly,tstep1,tstep3,natoms)
        # compute the force of the particles
        fx, fy = get_forces(x2,y2,lx,ly,natoms,nfil)
        # bin densities, velocities, and forces
        xedges, yedges, rho_hist[i], vx_hist[i], vy_hist[i], px_hist[i], py_hist[i], fx_hist[i], fy_hist[i] = bin_data(x2,y2,vx,vy,fx,fy,natoms,lx,ly,bx,by)
        # compute worticity
        worticity[i] = compute_rotation(vx_hist[i],vy_hist[i],bx,by,lx,ly)
        # compute momticity
        momticity[i] = compute_rotation(px_hist[i],py_hist[i],bx,by,lx,ly)
        # compute forticity
        forticity[i] = compute_rotation(fx_hist[i],fy_hist[i],bx,by,lx,ly)

        # add current step to time array
        time[i] = tstep2

    ### close input files
    ifile.close()
    hfile.close()
    return nsteps, time, xedges, yedges, rho_hist, vx_hist, vy_hist, px_hist, py_hist, fx_hist, fy_hist, worticity, momticity, forticity

##################################################################

def write_table(nsteps, time, xedges, yedges, array, fname):
    """ write table to file"""
    ofile = open(fname, 'w')
    nx = len(xedges) - 1
    ny = len(yedges) - 1
    # write header information
    ofile.write('nsteps = ' + str(nsteps) + '\n')
    ofile.write('nx = ' + str(nx) + '\n')
    ofile.write('ny = ' + str(ny) + '\n')
    # write edges of the bins
    ofile.write('xedges:\n')
    for i in range(nx+1):
        ofile.write(str(xedges[i]) + '\t')
    ofile.write('\n')
    ofile.write('yedges:\n')
    for i in range(ny+1):
        ofile.write(str(yedges[i]) + '\t')
    ofile.write('\n')
    # loop over all steps and write data
    for i in range(nsteps):
        ofile.write('timestep = ' + str(time[i]) + '\n')
        for j in range(nx):
            for k in range(ny):
                ofile.write(str(array[i,j,k]) + '\t')
            ofile.write('\n')
    ofile.close()
    return

##################################################################

def main():
    """ main function, called when script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip, nbins = read_settings()
    ### analyse the data
    nsteps, time, xedges, yedges, rho, vx, vy, px, py, fx, fy, worticity, momticity, forticity= analyse_data(charfile, headerfile, nfil, nskip, nbins)
    ### create folders to store data
    os.system('mkdir ' + ofname)
    os.system('mkdir ' + ofname + '/tables')
    os.system('mkdir ' + ofname + '/figs')
    write_table(nsteps, time, xedges, yedges, rho, ofname + '/tables/rho.data')
    write_table(nsteps, time, xedges, yedges, vx, ofname + '/tables/vx.data')
    write_table(nsteps, time, xedges, yedges, vy, ofname + '/tables/vy.data')
    write_table(nsteps, time, xedges, yedges, px, ofname + '/tables/px.data')
    write_table(nsteps, time, xedges, yedges, py, ofname + '/tables/py.data')
    write_table(nsteps, time, xedges, yedges, fx, ofname + '/tables/fx.data')
    write_table(nsteps, time, xedges, yedges, fy, ofname + '/tables/fy.data')
    write_table(nsteps, time, xedges, yedges, worticity, ofname + '/tables/worticity.data')
    write_table(nsteps, time, xedges, yedges, momticity, ofname + '/tables/momticity.data')
    write_table(nsteps, time, xedges, yedges, forticity, ofname + '/tables/forticity.data')
    
#    # create plots of the data
#    for i in range(nsteps):
#        fig = plt.figure(figsize = (10,5))
#        ax1 = plt.subplot(251)
#        ax2 = plt.subplot(252)
#        ax3 = plt.subplot(253)
#        ax4 = plt.subplot(254)
#        ax5 = plt.subplot(255)
#        ax6 = plt.subplot(256)
#        ax7 = plt.subplot(257)
#        ax8 = plt.subplot(258)
#        ax9 = plt.subplot(259)
#        ax10 = plt.subplot(2,5,10)    
#        # add the data
#        ax1.imshow(rho[i].transpose(), origin = 'lower', interpolation = 'None')
#        ax1.set_title('Density')
#        
#        ax2.imshow(vx[i].transpose(), origin = 'lower', interpolation = 'None')
#        ax2.set_title('vx')
#
#        ax3.imshow(vy[i].transpose(), origin = 'lower', interpolation = 'None')
#        ax3.set_title('vy')
#
#        ax4.imshow(worticity[i].transpose(), origin = 'lower', interpolation = 'None')
#        ax4.set_title('worticity')
#
#        ax5.imshow(fx[i].transpose(), origin = 'lower', interpolation = 'None')
#        ax5.set_title('fx')
#
#        ax6.imshow(fy[i].transpose(), origin = 'lower', interpolation = 'None')
#        ax6.set_title('fy')
#
#        ax7.imshow(forticity[i].transpose(), origin = 'lower', interpolation = 'None')
#        ax7.set_title('forticity')
#
#        ax8.imshow(px[i].transpose(), origin = 'lower', interpolation = 'None')
#        ax8.set_title('px')
#
#        ax9.imshow(py[i].transpose(), origin = 'lower', interpolation = 'None')
#        ax9.set_title('py')
#
#        ax10.imshow(momticity[i].transpose(), origin = 'lower', interpolation = 'None')
#        ax10.set_title('momticity')
#    
#        # save and close the plots
#        plt.savefig(ofname + '/figs/step_' + str(time[i]) + '.png')
#        plt.close()
    return

##################################################################

if __name__ == '__main__':
    main()
