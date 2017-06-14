#!/usr/local/bin/python2.7

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import codecs
import read_char
import os
import performance_toolsWrapper

try:
    infilename = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '       parameter file'
    exit()

##################################################################

def read_settings():
    """ read in the settings from the parameter file"""
    # open file for reading
    ifile = open(infilename, 'r')
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
    # number of points in the MSD function
    line = ifile.readline()
    line = line.split()
    nmsd = int(line[-1])
    # maximum value of dt compared to entire simulation time
    line = ifile.readline()
    line = line.split()
    limit = float(line[-1])
    # read com every that many steps
    line = ifile.readline()
    line = line.split()
    nevery = int(line[-1])
    # close file
    ifile.close()
    # return input values
    return charfile, headerfile, ofname, nfil, nskip, nmsd, limit, nevery

##################################################################

def neigh_min(dx,lx):
    ### compute minimum distance
    dx1 = dx + lx
    dx2 = dx - lx
    if dx**2 < dx1**2 and dx**2 < dx2**2:
        return dx
    if dx1**2 < dx2**2:
        return dx1
    return dx2

##################################################################

def compute_com_single(x,y,lx,ly):
    """ compute the absolute of the spiral number of a filament,
        correct pbc on the fly"""
    nbeads = len(x)
    # correct pbcs
    for i in range(1,nbeads):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        x[i] = x[i-1] + neigh_min(dx,lx)
        y[i] = y[i-1] + neigh_min(dy,ly)
    # compute comx and comy
    comx = np.average(x)
    comy = np.average(y)
    # put comx and comy into simulation box
    comx /= lx
    comx = comx - math.floor(comx)
    comx *= lx
    comy /= ly
    comy = comy - math.floor(comy)
    comy *= ly
    # return center of mass
    return comx, comy

##################################################################

def compute_com(charfile, headerfile, nfil, nskip):
    """ compute the center of mass of all filaments"""
    ### open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
    ### skip initial snapshots
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    ### allocate arrays, including output: spiral histogram + spiral counter
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    # center of mass coordinates
    comx = np.zeros((nmol, nsteps), dtype = np.float64)
    comy = np.zeros((nmol, nsteps), dtype = np.float64)
    # displacement between two timesteps
    dx = np.zeros((nsteps), dtype = np.float64)
    dy = np.zeros((nsteps), dtype = np.float64)
    # time
    time = np.zeros((nsteps), dtype = np.int32)
    ### loop over all timesteps
    for step in range(nsteps):
        ### print progress to screen
        print 'Current Step / All Steps', step, '/', nsteps
        ### read in coordinates
        xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
        x = xs*lx
        y = ys*ly
        ### store time
        time[step] = tstep
        ### compute com for each filament, correct pbc on the fly
        for i in range(nmol):
            comx[i,step], comy[i,step] = compute_com_single(x[i*nfil:(i+1)*nfil-1], y[i*nfil:(i+1)*nfil-1], lx, ly)
    ### close the input files
    hfile.close()
    ifile.close()
    ### correct pbc errors
    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.correct_pbc(comx,dx,lx,nsteps,nmol)
    performance_tools.correct_pbc(comy,dy,ly,nsteps,nmol)
    return time, comx, comy

##################################################################

def compute_msd(time, comx, comy, nmsd, limit, nevery):
    """ compute and average the MSD"""
    ### nsteps and nmol
    nmol = len(comx)
    nsteps = len(comx[0])
    ### allocate output array to store the results
    t_msd = np.zeros((nmsd), dtype = np.int32)
    msd = np.zeros((nmsd), dtype = np.float64)
    msd_tmp = np.zeros((nmsd), dtype = np.float64)
    counter = np.zeros((nmsd), dtype = np.int32)
    logval = np.zeros((nmsd), dtype = np.int32)
    ### compute the maximum number of steps
    limit = int(limit*nsteps)
    ### loop over all molecules, submit MSD computation, read and average
    ###   the results
    performance_tools = performance_toolsWrapper.Performance_tools()
    for i in range(nmol):
        print 'Computing MSD', i, '/', nmol
        # submit MSD computation
        performance_tools.compute_MSD(t_msd, time, comx[i], comy[i], msd_tmp, counter, logval, nsteps, nmsd, limit, nevery)     
        # average the results using the on-line algorithm
        delta = msd_tmp - msd
        msd += delta/(i+1)
    return t_msd, msd

##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip, nmsd, limit, nevery = read_settings()
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    time, comx, comy = compute_com(charfile, headerfile, nfil, nskip)
    ### start msd computations with cpp code for all files and average
    t_msd, msd = compute_msd(time, comx, comy, nmsd, limit, nevery)
    ### write results to file and generage a plot
    #generate subfolder
    os.system('mkdir -p ' + ofname)
    ofile = open(ofname + '/msd.data', 'w')
    ofile.write('Averaged MSD over all molecules\n\n')
    ofile.write('Timestep\tMSD\n')
    for i in range(len(t_msd)):
        ofile.write(str(t_msd[i]) + '\t' + str(msd[i]) + '\n')
    ofile.close()
#    fig = plt.figure()
#    ax = plt.subplot(111)
#    ax.loglog(t_msd,msd)
#    ax.set_xlabel('t [steps]')
#    ax.set_ylabel(r'msd [$\sigma^2$]')
#    ax.set_title('MSD')
#    plt.savefig(ofname + '/msd.png')
#    plt.close()
    
##################################################################

if __name__ == '__main__':
    main()
    
