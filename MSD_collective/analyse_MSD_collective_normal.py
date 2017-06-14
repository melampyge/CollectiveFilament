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
import h5py

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
    
    return charfile, headerfile, ofname, nfil, nskip, nmsd, limit, nevery

##################################################################

def neigh_min(dx,lx):
    """ compute minimum image distance"""
    
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
    
    return time, comx, comy, nsteps, natoms, nmol

##################################################################

def write_com_data(time, comx, comy, path, nsteps, natoms, nmol):
    """ write pbc corrected com data as a function of time"""
    
    f = h5py.File(path, 'w')
    f.create_dataset('time', data=time, dtype='int32', compression='gzip')
    f.create_dataset('comx', (nmol, nsteps), data=comx, dtype='float64', compression='gzip')   
    f.create_dataset('comy', (nmol, nsteps), data=comy, dtype='float64', compression='gzip')    
    f.close()
    
    return

##################################################################

def read_com_data(path):
    """ read pbc corrected com data as a function of time"""
    
    f = h5py.File(path, 'r')
    
    time = f['time']
    comx = f['comx']
    comy = f['comy']
    
    time = np.asarray(time, dtype=np.int32)
    comx = np.asarray(comx, dtype=np.float64)
    comy = np.asarray(comy, dtype=np.float64)
    
    nmol, nsteps = np.shape(comx)
            
    return time, comx, comy, nsteps, nmol

##################################################################

def write_msd_data(delay, msd, path):
    """ write msd data"""
    
    f = open(path, 'w')
    for j in range(len(delay)):
        f.write(str(delay[j]) + '\t\t' + str(msd[j]) + '\n')
        
    f.close()
    
    return
    
##################################################################

def compute_msd(time, comx, comy, nsteps, nmol):
    """ compute and average the MSD"""
        
    ndelay = int(nsteps/2)
    delay = np.zeros((ndelay), dtype = np.int32)
    msd = np.zeros((ndelay), dtype = np.float64)
    
    for d in range(1,ndelay):
        delay[d] = d*5
        msd[d] = np.mean(np.mean((comx[:,d:]-comx[:,:-d])**2 + (comy[:,d:]-comy[:,:-d])**2,axis=0))
        
    return delay, msd

##################################################################

def main():
    
    ### read parameters from input file
    
    charfile, headerfile, ofname, nfil, nskip, nmsd, limit, nevery = read_settings()
    
#    ### loop over the entire trajectory and compute the center of mass, correct for pbc
#    
#    time, comx, comy, nsteps, natoms, nmol  = compute_com(charfile, headerfile, nfil, nskip)
    
    ### write the pbc corrected coms as a fnc of time
    
#    ofilepath = ofname + '/com.data'
#    write_com_data(time, comx, comy, ofilepath, nsteps, natoms, nmol)
    
    ### compute msd
    
    ifilepath = ofname + '/com.data'
    time, comx, comy, nsteps, nmol = read_com_data(ifilepath)
    delay, msd = compute_msd(time, comx, comy, nsteps, nmol)
    
    ### write msd data
    
    ofilepath = ofname + '/msd_alternative.data'
    write_msd_data(delay, msd, ofilepath)
    
   
##################################################################

if __name__ == '__main__':
    main()
    
