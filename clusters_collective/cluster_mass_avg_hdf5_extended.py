
################################################################################

import argparse
import numpy as np
import os
import h5py
 
################################################################################
 
def loadSimData(datafile):
    """ load initial simulation data"""
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
        dtSamp, T, box_area, nt, body_length, Pe, persistence, flexure, tau_D, tau_A, \
            nbeads, nmol, nfil

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
    box_area = Lx*Ly
    body_length = B*(N-1)
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    flexure = Pe/persistence
    T = totalStep - ti
    nt = T/nsamp
    nsamp = int(nsamp)
    nbeads = int(L)
    nfil = int(N)
    nmol = int(M)
    tau_D = body_length**2*(N+1)/4/kT
    if Fmc == 0:
        Fmc = 0.0000001
    tau_A = (N+1)/Fmc
    print "Diffusive time scale is ", tau_D
    print "Advective time scale is ", tau_A
    
    return
        
################################################################################        
        
class ClusterSize:
    """ Data structure for cluster sizes"""
    
    def __init__(self, path):
        f = h5py.File(path, 'r')
        size_grp = f['size']
        self.cs = np.asarray(size_grp['size'], dtype=int)       # Cluster sizes
        self.totcs = len(self.cs)                               # Total number of clusters
        f.close()
     
################################################################################ 


def main(): 

    ### argument parsing (command line options)

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", help="Folder containing simulation data")
    args = parser.parse_args()
    
    ### load saved preliminary simulation data into relevant variables
    
    loadSimData(args.folder+"init_info.txt")
    tic = int(ti)
    tf = int(totalStep)
    nsampc = int(nsamp)
    total_num_frames = (tf - tic)/nsampc
    print "ti = ", tic, " tf = ", tf, " nsamp = ", nsampc
    print "Total number of frames = ", total_num_frames
    
    datafolder = args.folder+"CLUSTER/"
    
    num_of_clusters_at_size_at_time = np.zeros((total_num_frames, nmol+1), dtype=np.float64)
    frame_cnt = -1
    
    ### count the number of clusters with a given size at each time frame
    
    for frame in np.arange(tic, tf, nsampc):
        
        print frame , ' of ', totalStep
        
        ### load cluster sizes at each time frame
        
        path = datafolder + 'cluster_' + str(frame) + '.hdf5'
        
        if os.path.exists(path):
            frame_cnt += 1
            clsizes = ClusterSize(path)
            
            ### histogram the cluster sizes
                    
            for size in clsizes.cs:
                
                size_idx = int(size)
                num_of_clusters_at_size_at_time[frame_cnt][size_idx] += 1.
            
        else:
            print 'Cluster size data does NOT exist for ', path
            
    ### compute cluster size distribution -- take a time average
            
    cluster_size_dist = np.mean(num_of_clusters_at_size_at_time, axis=0)
    thalf = len(num_of_clusters_at_size_at_time)/2
    cluster_size_dist_first_half = np.mean(num_of_clusters_at_size_at_time[:thalf], axis=0)
    cluster_size_dist_sec_half = np.mean(num_of_clusters_at_size_at_time[thalf:], axis=0)
    
    ### save the data
        
    total_num_frames_after = frame_cnt + 1
    print "Total number of frames before / after \t", total_num_frames, total_num_frames_after
    if total_num_frames_after <= 0:
        total_num_frames_after = 1.0
        
    savefile = datafolder + "avg_cl_size_cnt.txt"
    sfile = open(savefile, 'w')

    avg_mass = 0.0
    avg_mass_first_half = 0.0
    avg_mass_sec_half = 0.0
    for j in np.arange(1,nmol+1):
       
        cluster_size_dist[j] = cluster_size_dist[j]*j/nmol
        cluster_size_dist_first_half[j] = cluster_size_dist_first_half[j]*j/nmol
        cluster_size_dist_sec_half[j] = cluster_size_dist_sec_half[j]*j/nmol        
        sfile.write(str(j) + "\t" + str(cluster_size_dist[j]) + "\n")
        avg_mass += cluster_size_dist[j]*j
        avg_mass_first_half += cluster_size_dist_first_half[j]*j
        avg_mass_sec_half += cluster_size_dist_sec_half[j]*j
        
    sfile.close()
    
    avg_mass_halves = [avg_mass_first_half, avg_mass_sec_half]
    avg_mass_std = np.std(avg_mass_halves)/np.sqrt(2)

    savefile = datafolder + "avg_size_extended.txt"
    sfile = open(savefile, 'w')
    sfile.write(str(avg_mass) + '\t\t' + str(avg_mass_std) + '\n')    
    sfile.close()
    
################################################################################
    
if __name__ == '__main__':
    main()

################################################################################