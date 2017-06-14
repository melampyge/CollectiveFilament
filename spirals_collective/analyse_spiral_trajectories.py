
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import math
import os
import h5py
import argparse

##############################################################################

def get_sim_info(path):
    """ load initial simulation data"""
    
    f = h5py.File(path, 'r')
    
    info = f['info']
    box = info['box']
    
    lx = box['lx']
    lx = lx[()]
    ly = box['ly']
    ly = ly[()]
    
    nsteps = info['nsteps']
    nsteps = nsteps[()]
    natoms = info['natoms']
    natoms = natoms[()]
    nmol = info['nmol']
    nmol = nmol[()]
    
    #f.close()
    
    return lx, ly, nsteps, natoms, nmol

##############################################################################

def load_data(path):
    """ load simulation data"""
    
    f = h5py.File(path, 'r')
    
    pos = f['pos']
    x = pos['x']
    y = pos['y']
    
    time = f['time']
  
    #f.close()
    
    return time, x, y
    
##############################################################################

def neigh_min(dx, lx):
    """ compute minimum distance"""
    
    dx1 = dx + lx
    dx2 = dx - lx
    if dx**2 < dx1**2 and dx**2 < dx2**2:
        return dx
    if dx1**2 < dx2**2:
        return dx1
    return dx2

##############################################################################

def compute_com(x, y, lx, ly):
    """ compute the center of mass of filaments,
        correct pbc on the fly"""
    
    ### NOT CORRECTED FOR UNWRAPPED COORDINATES!!!
    
    nbeads = len(x)

    comx = x[0]
    comy = y[0]
    
    # correct pbcs
    for i in range(1,nbeads):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        x[i] = x[i-1] + neigh_min(dx,lx)
        y[i] = y[i-1] + neigh_min(dy,ly)
        
    # compute center of mass
    for i in range(1,nbeads):
        comx += x[i]
        comy += y[i]
        
    comx /= nbeads
    comy /= nbeads
    
    return comx, comy
    
##############################################################################
    
def analyse_com_trajectory(xt, yt, lx, ly, nsteps, nmol, nfil):
    
    comx = ((nsteps, nmol))
    comy = ((nsteps, nmol))
    
    for tstep in range(nsteps):
        
        print "step ", tstep, " / nsteps ", nsteps
        
        x = xt[tstep]
        y = yt[tstep]
        
        for i in range(nmol):
            comx[tstep][i], comy[tstep][i] = compute_com(x[i*nfil:(i+1)*nfil-1], \
                y[i*nfil:(i+1)*nfil-1], lx, ly) 

    return comx, comy 
            
##############################################################################

def main():
    
    ### generate contextual info about the data
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", help="Path to data file as in density_0.2/kappa_5.0/fp_1.0/")
    args = parser.parse_args()  
    fname = 'out1.hdf5'
    base = '/local/duman/SIMULATIONS/long_filaments/'
    path = base + args.folder + fname
    lx, ly, nsteps, natoms, nmol = get_sim_info(path)
    nfil = 201      ## number of beads per filaments
  
    ### load trajectories
    
    time, x, y = load_data(path)
    
    ### print out general information about the simulation
    
    print "lx = ", lx, " / ly = ", ly
    print "nsteps = ", nsteps
    print "natoms = ", natoms, " / nmol = ", nmol
    print "nfil = ", nfil 
    
    comx, comy = analyse_com_trajectory(x, y, lx, ly, nsteps, nmol, nfil)
    
    ### write results to files
    
    ofname = 'SPIRAL/'
    ofolder = base + args.folder + ofname
    os.system('mkdir -p ' + ofolder)
    ofilename = ofolder + 'trajectories.hdf5'
    
    ofile = h5py.File(ofilename, 'w')
    ofile.create_dataset('num_spirals', data=num_spirals, dtype='float64')  
    ofile.create_dataset('time', data=time, dtype='float64')    
    ofile.close()
    
    return
    
    
##############################################################################

if __name__ == '__main__':
    main()

##############################################################################    
