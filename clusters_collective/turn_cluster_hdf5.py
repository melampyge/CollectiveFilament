#!/usr/bin/python

# Load needed libraries and necessary files
import argparse
import numpy as np
import os
import h5py
import gzip

################################################################################

#
# Function definitions
#
 
################################################################################
 
def loadSimData(datafile):
    """ Load initial simulation data"""
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
    dtSamp, T, box_area, nt, body_length, Pe, persistence, flexure 

    datafile = open(datafile,"r")
    for line in datafile:
        A = line.split()
        if A[0] == "dt":                    # Time interval between MD steps
            dt = float(A[-1])
        elif A[0] == "ti":                  # Beginning time for data acquisition
            ti = float(A[-1])
        elif A[0] == "Lx":                  # Box size in x
            Lx = float(A[-1])            
        elif A[0] == "Ly":                  # Box size in y
            Ly = float(A[-1])
        elif A[0] == "totalStep":           # Total MD steps
            totalStep = float(A[-1])
        elif A[0] == "nsamp":               # Data sampling frequency
            nsamp = float(A[-1])
        elif A[0] == "nfil":                # Number of particles per polymer
            N = float(A[-1])
        elif A[0] == "L":                   # Number of particles
            L = float(A[-1])
        elif A[0] == "B":                   # Bond length between particles of a body
            B = float(A[-1])
        elif A[0] == "kT":                  # Boltzmann constant*Temperature
            kT = float(A[-1])
        elif A[0] == "Fmc":                 # Self propulsion force constant
            Fmc = float(A[-1])     
        elif A[0] == "Kbend":               # Bending constant
            Kbend = float(A[-1])
    
    Lx /= B
    Ly /= B
    M = L/N
    dtSamp = dt*nsamp
    T = totalStep - ti
    nt = T/nsamp
    box_area = Lx*Ly
    body_length = B*(N-1)
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    flexure = Pe/persistence

################################################################################            

def loadSizes(dataf, address, frame):
    """ Load cluster sizes"""
    
    path = dataf + '/cluster_' + address + '_' + str(frame) + '.txt' 
    f_gz = np.loadtxt(path, dtype=int)
    cs = np.transpose(f_gz)

    return cs    
    
################################################################################
    
def loadBeads(dataf, address, frame):
    """ Load bead information"""
    
    path = dataf + '/' + address + '_' + str(frame) + '.txt' 
    f_gz = np.transpose(np.loadtxt(path, dtype=float))
    xi = f_gz[0]
    yi = f_gz[1]
    phi = f_gz[2]
    cidx = f_gz[3]

    return xi, yi, phi, cidx    
    
################################################################################               

def loadCluster(dataf, address, frame):
    """ Load cluster information"""
    
    path = dataf + '/cluster_' + address + '_' + str(frame) + '.txt' 
    
    idx = []
    size = []
    comx = []
    comy = []
    rgysq = []
    pl = []
    st = []
    sw = []
    ens = []
    fils = []
    
    #f_gz = gzip.open(path, 'rt')
    f_gz = open(path, 'r')    
    for line in f_gz:
        A = line.split()
        idx.append(int(A[0]))
        siz = int(A[1])
        size.append(siz)
        fil_list = np.zeros((siz))
        comx.append(float(A[2])/B)
        comy.append(float(A[3])/B)
        rgysq.append(float(A[4]))
        pl.append(float(A[5]))
        st.append(float(A[6]))
        sw.append(float(A[7]))
        ens.append(float(A[8]))
        for el in range(siz):
            fil_list[el] = int(A[9+el])    
        fils.append(fil_list)    
        
    f_gz.close()

    return idx, size, comx, comy, rgysq, pl, st, sw, ens, fils   

################################################################################

def saveHDF5(ifile, cs, xi, yi, phi, cidx, idx, size, comx, comy, rgysq, pl, st, sw, ens, fils):
    """ Save in HDF5 format"""
    
    ## Groups
    size_grp = ifile.create_group('size')
    beads_grp = ifile.create_group('beads')
    props_grp = ifile.create_group('props')
    
    ## Datasets
    size_grp.create_dataset('size', data=cs, compression='gzip')
    
    beads_grp.create_dataset('x', data=xi, compression='gzip')
    beads_grp.create_dataset('y', data=yi, compression='gzip')
    beads_grp.create_dataset('phi', data=phi, compression='gzip')
    beads_grp.create_dataset('idx', data=cidx, compression='gzip')
    
    idx_h5 = np.asarray(idx, dtype=int)
    size_h5 = np.asarray(size, dtype=int)
    comx_h5 = np.asarray(comx, dtype=float)
    comy_h5 = np.asarray(comy, dtype=float)
    rgysq_h5 = np.asarray(rgysq, dtype=float)
    planarity_h5 = np.asarray(pl, dtype=float)
    straightness_h5 = np.asarray(st, dtype=float)
    swirliness_h5 = np.asarray(sw, dtype=float)
    enstrophy_h5 = np.asarray(ens, dtype=float)
           
    props_grp.create_dataset('idx', data=idx_h5, compression='gzip')
    props_grp.create_dataset('size', data=size_h5, compression='gzip')
    props_grp.create_dataset('comx', data=comx_h5, compression='gzip')
    props_grp.create_dataset('comy', data=comy_h5, compression='gzip')
    props_grp.create_dataset('rgysq', data=rgysq_h5, compression='gzip')
    props_grp.create_dataset('planarity', data=planarity_h5, compression='gzip')
    props_grp.create_dataset('straightness', data=straightness_h5, compression='gzip')
    props_grp.create_dataset('swirliness', data=swirliness_h5, compression='gzip')
    props_grp.create_dataset('enstrophy', data=enstrophy_h5, compression='gzip')
    
    ## Filament list
    fil_grp = props_grp.create_group('filament_list')
    for sz_idx in np.arange(len(size_h5)):
        fil_list = np.asarray(fils[sz_idx], dtype=int)
        fil_grp.create_dataset(str(sz_idx), data=fil_list, compression='gzip')
     
    return

################################################################################

def loadHDF5(ofile):
    """ Load HDF5 file for testing purposes"""
    
    ## Groups    
    size_grp = ofile['size']
    beads_grp = ofile['beads']
    props_grp = ofile['props']
    
    ## Datasets
    size = np.asarray(size_grp['size'])
    x = np.asarray(beads_grp['x'])
    comx = np.asarray(props_grp['comx'])
    
    ## Filament list
    fil_grp = props_grp['filament_list']
    fil_list_1 = np.asarray(fil_grp['1'])
    
    print size
    print x
    print comx
    print fil_list_1
    
    return 

################################################################################
    
def main(): 

    ## Argument parsing (command line options)
    parser = argparse.ArgumentParser()
    parser.add_argument("simfolder", help="File containing simulation data")
    args = parser.parse_args()
    
    ## Load saved preliminary simulation data into relevant variables
    loadSimData('/local/duman/SIMULATIONS/many_polymers_5/density_0.08/kappa_5.0/fp_0.24/init_info.txt')    
    
    ## Address of the folder containing simulation data
    datafolder = args.simfolder+"/CLUSTER"
    
    ## Time frame loop
    for frame in np.arange(int(ti),int(totalStep),int(nsamp)):
        
        print frame , ' of ', totalStep
        
        ## LOAD
        
        ## Load cluster sizes
        cs = loadSizes(datafolder, 'sizes', frame)
        
        ## Load beads information
        xi, yi, phi, cidx = loadBeads(datafolder, 'beads', frame)

        ## Load cluster information
        idx, size, comx, comy, rgysq, pl, st, sw, ens, fils = loadCluster(datafolder, 'evolution', frame)
        
        ## SAVE
        ifile = h5py.File(datafolder + '/cluster_' + str(frame) + '.hdf5', 'w')
        saveHDF5(ifile, cs, xi, yi, phi, cidx, idx, size, comx, comy, rgysq, pl, st, sw, ens, fils)
        ifile.close()
        

#    ## TEST
#    ofile = h5py.File(datafolder + '/cluster_10000000.hdf5', 'r')
#    loadHDF5(ofile)
#    ofile.close()
        


    
if __name__ == '__main__':
    main()
    
    