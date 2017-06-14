
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

def compute_spiral_number(x, y, lx, ly):
    """ compute the absolute of the spiral number of a filament,
        correct pbc on the fly"""
        
    nbeads = len(x)
    
    # correct pbcs
    for i in range(1,nbeads):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        x[i] = x[i-1] + neigh_min(dx,lx)
        y[i] = y[i-1] + neigh_min(dy,ly)
        
    # compute all bond orientations and center of mass
    phi = np.zeros((nbeads - 1))
    for i in range(1,nbeads):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        dphi = math.atan2(dy,dx)
        phi[i-1] = dphi
                        
    # correct for 2*pi periodicity
    phi2 = np.copy(phi)
    nbonds = len(phi)
    for i in range(1,nbonds):
        dphi = phi[i] - phi[i-1]
        if dphi < -np.pi:
            dphi += 2*np.pi
        elif dphi > np.pi:
            dphi -= 2*np.pi
        phi2[i] = phi2[i-1] + dphi
        
    # compute the spiral number
    s = (phi2[-1] - phi2[0])/2/np.pi
    
    return s
        
##############################################################################

def compute_num_of_spirals(xt, yt, lx, ly, nsteps, nmol, nfil):
    """ compute number of spirals as a fnc. of time"""
    
    s7 = np.zeros((nsteps), dtype='float64')
    s5 = np.zeros((nsteps), dtype='float64')
    s3 = np.zeros((nsteps), dtype='float64')
    s2 = np.zeros((nsteps), dtype='float64')
    
    for tstep in range(nsteps):
        
        print 'step ', tstep, ' / nsteps ', nsteps
        
        x = xt[tstep]
        y = yt[tstep]
        
        for i in range(nmol):
        
            s = compute_spiral_number(x[i*nfil:(i+1)*nfil-1], \
                y[i*nfil:(i+1)*nfil-1], lx, ly)     
            s = math.fabs(s)
            if s > 6.0:
                s7[tstep] += 1
            elif s >= 1.6 and s < 2.2:
                s2[tstep] += 1
            elif s >= 2.2 and s < 3.5:
                s3[tstep] += 1
            elif s >= 3.5 and s < 6.0:
                s5[tstep] += 1
            
    return s2, s3, s5, s7
    
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
    
    s2, s3, s5, s7 = compute_num_of_spirals(x, y, lx, ly, nsteps, nmol, nfil)
    
    ### write results to files
    
    ofname = 'SPIRAL/'
    ofolder = base + args.folder + ofname
    os.system('mkdir -p ' + ofolder)
    ofilename = ofolder + 'num_spirals.hdf5'
    
    ofile = h5py.File(ofilename, 'w')
    ofile.create_dataset('s2', data=s2, dtype='float64')  
    ofile.create_dataset('s3', data=s3, dtype='float64')  
    ofile.create_dataset('s5', data=s5, dtype='float64')  
    ofile.create_dataset('s7', data=s7, dtype='float64')      
    ofile.create_dataset('time', data=time, dtype='float64')    
    ofile.close()
    
    return
    
    
##############################################################################

if __name__ == '__main__':
    main()

##############################################################################    
