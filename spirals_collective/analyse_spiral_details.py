
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
    
def compute_spiral_number_and_com(x, y, lx, ly):
    """ compute the absolute of the spiral number of a filament and center of mass,
        correct pbc on the fly"""
        
    nbeads = len(x)
    comx = x[0]
    comy = y[0]
    
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
        
        comx += x[i]
        comy += y[i]
        
    comx /= nbeads
    comy /= nbeads
        
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
    
    return s, comx, comy
    
##############################################################################

def compute_tagged_info(sid, fid, mode2, mode3, mode5, xt, yt, lx, ly, nsteps, nmol, nfil):
    """ compute and save information about the tagged filaments"""
    
    s = np.zeros((5, nsteps))
    comx = np.zeros((5, nsteps))
    comy = np.zeros((5, nsteps))
    
    m2 = mode2[0]
    m3 = mode3[0]
    m5 = mode5[0]
    
    print "\nCOMPUTING TAGGED INFO"
    
    for tstep in range(nsteps):
        
        print 'step ', tstep, ' / nsteps ', nsteps
        
        x = xt[tstep]
        y = yt[tstep]
        
        s[0][tstep], comx[0][tstep], comy[0][tstep] = compute_spiral_number_and_com(x[sid*nfil:(sid+1)*nfil-1], \
            y[sid*nfil:(sid+1)*nfil-1], lx, ly)     
        s[1][tstep], comx[1][tstep], comy[1][tstep] = compute_spiral_number_and_com(x[fid*nfil:(fid+1)*nfil-1], \
            y[fid*nfil:(fid+1)*nfil-1], lx, ly)
        s[2][tstep], comx[2][tstep], comy[2][tstep] = compute_spiral_number_and_com(x[m2*nfil:(m2+1)*nfil-1], \
            y[m2*nfil:(m2+1)*nfil-1], lx, ly)     
        s[3][tstep], comx[3][tstep], comy[3][tstep] = compute_spiral_number_and_com(x[m3*nfil:(m3+1)*nfil-1], \
            y[m3*nfil:(m3+1)*nfil-1], lx, ly)   
        s[4][tstep], comx[4][tstep], comy[4][tstep] = compute_spiral_number_and_com(x[m5*nfil:(m5+1)*nfil-1], \
            y[m5*nfil:(m5+1)*nfil-1], lx, ly)     
            
    return s, comx, comy

##############################################################################

def detect_first_spiral(xt, yt, lx, ly, nsteps, nmol, nfil):
    """ detect the identity of the first spiral"""
    
    ### compute spiral number for each filament, correct pbc on the fly
    
    for tstep in range(nsteps):
        
        x = xt[tstep]
        y = yt[tstep]
        
        for i in range(nmol):
            
            s = math.fabs(compute_spiral_number(x[i*nfil:(i+1)*nfil-1], \
                y[i*nfil:(i+1)*nfil-1], lx, ly))  
                
            if s > 6.:
                
                print "\nFIRST SPIRAL DETECTED"
                print "tstep = ", tstep
                print "filament id = ", i
                print "spiral number = ", s
                
                return tstep, i
        
    return

##############################################################################

def detect_spirals(xt, yt, lx, ly, nsteps, nmol, nfil):
    """ detect the identity of the several spiral modes at the end of the simulation"""
    
    ### compute spiral number for each filament, correct pbc on the fly
    
    T = len(xt)
    s = np.zeros((nmol))

    for tstep in range(T):
        
        x = xt[tstep]
        y = yt[tstep]
            
        for i in range(nmol):
            
            s[i] += math.fabs(compute_spiral_number(x[i*nfil:(i+1)*nfil-1], \
                y[i*nfil:(i+1)*nfil-1], lx, ly))
                
    s /= nmol
            
    s_min = min(s)
    s_min_id = np.argmin(s)
    
    s_mode_2 = []
    s_mode_3 = []
    s_mode_5 = []
    for j in range(nmol):
        if s[j] >= 1.6 and s[j] < 2.2:
            s_mode_2.append(j)
        elif s[j] >= 2.2 and s[j] < 3.5:
            s_mode_3.append(j)
        elif s[j] >= 3.5 and s[j] < 6.0:
            s_mode_5.append(j)
    
    print "\nFREE FILAMENT DETECTED"
    print "filament id = ", s_min_id
    print "s = ", s_min
            
    return s_min_id, s_mode_2, s_mode_3, s_mode_5

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
    
    ### get the id of the first spiral
    
    #tstep, sid = detect_first_spiral(x, y, lx, ly, nsteps, nmol, nfil)
    
    ### get the id of the free spiral at the end of the simulation
    
    slice_index = nsteps-1500
    fid, s_mode_2, s_mode_3, s_mode_5 = detect_spirals(x[slice_index:], y[slice_index:], lx, ly, nsteps, nmol, nfil)

    ofname = 'SPIRAL/'
    ofolder = base + args.folder + ofname
    
    mode2file = open(ofolder + 'mode_2.data', 'w')
    for j in range(len(s_mode_2)):
        mode2file.write(str(s_mode_2[j]) + '\n')
    mode2file.close()
    
    mode3file = open(ofolder + 'mode_3.data', 'w')
    for j in range(len(s_mode_3)):
        mode3file.write(str(s_mode_3[j]) + '\n')  
    mode3file.close()        
        
    mode5file = open(ofolder + 'mode_5.data', 'w')
    for j in range(len(s_mode_5)):
        mode5file.write(str(s_mode_5[j]) + '\n') 
    mode5file.close()
    
    ### compute information about the tagged filaments
    sid = 1176
    #fid = 621
    spiral_number, comx, comy = compute_tagged_info(sid, fid, s_mode_2, s_mode_3, s_mode_5, x, y, lx, ly, nsteps, nmol, nfil)
    
    ### write results to files
    
    os.system('mkdir -p ' + ofolder)
    ofilename = ofolder + 'tagged_filament_info.hdf5'
    
    ofile = h5py.File(ofilename, 'w')
    ofile.create_dataset('comx', (5, nsteps), data=comx, dtype='float64')
    ofile.create_dataset('comy', (5, nsteps), data=comy, dtype='float64')
    ofile.create_dataset('spiral_number', (5, nsteps), data=spiral_number, dtype='float64')  
    ofile.create_dataset('time', data=time, dtype='float64')    
    ofile.close()
    
    return
    
    
##############################################################################

if __name__ == '__main__':
    main()

##############################################################################    
