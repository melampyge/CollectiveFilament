
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
    datafolder = args.folder+"CLUSTER/"
    
    size_max = nmol+1
    size_min = 0
    bin_width = 10
    nbins = (size_max-size_min)/bin_width
    
    size_cnt = np.zeros((nbins), dtype=np.float64)
    cnt = 0
    
    for frame in np.arange(int(ti),int(totalStep),int(nsamp)):
        
        
        print frame , ' of ', totalStep
        cnt += 1.
        
        ### load cluster sizes at each time frame
        
        path = datafolder + 'cluster_' + str(frame) + '.hdf5'
        if os.path.exists(path):
            clsizes = ClusterSize(path)
            
            ### histogram the cluster sizes
            
            for size in clsizes.cs:
                
                seg = int(size/bin_width)
                size_cnt[seg] += 1
                if seg > 2:
                    print size, '\t\t', seg, '\t\t', size_cnt[seg]
        
    ### save the data
        
    savefile = datafolder + "avg_cl_size_cnt.txt"
    sfile = open(savefile, 'w')
    for j in np.arange(nbins):
        size_cnt[j] /= cnt
    
    sumt = np.sum(size_cnt)
    for j in np.arange(nbins):
        size_cnt[j] /= (bin_width*sumt)
        lower_edge = j*bin_width
        upper_edge = (j+1)*bin_width
        middle_pt = (upper_edge+lower_edge)/2.
        sfile.write(str(middle_pt) + "\t" + str(size_cnt[j]) + "\n")
        
    sfile.close()

################################################################################
    
if __name__ == '__main__':
    main()

################################################################################