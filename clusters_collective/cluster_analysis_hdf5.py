#!/usr/bin/python

# Load needed libraries and necessary files
import argparse
import numpy as np
import os
import h5py

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
      
def cluster_lifetime(cl_idx, old_cl_idx, lifetime, ids):
    """ Determine and make a histogram from the lifetime of clusters"""
    
    # If one of the clusters in the old cluster indices are lost, mark the number of steps it has taken to get lost
    # If the old cluster index is retained, increase the count by one more
    # ids[ cluster id j ] = number of steps the cluster with id j is living (dictionary)
    # lifetime[ number of steps of life times of clusters ] = counts (array)
    
    ## For all the clusters from the previous time step
    for j in range(len(old_cl_idx)):
            
        ## Search through all the new cluster indices to find a matching new cluster
        for nj in range(len(cl_idx)):
            
            ## A match: Old cluster index is retained, increase the counter (the cluster is living for one more step) and go to the next old cluster
            if old_cl_idx[j] == cl_idx[nj]:
                ids[old_cl_idx[j]] += 1
                break
            
        ## No match: Old cluster index is lost (cluster from previous time died!), increase the lifetime counter (mark the death of the old cluster!) and delete its key value
        else:
            ## Add the lifetime to the corresponding cluster size life counter
            ## ??
            lifetime[ids[old_cl_idx[j]]] += 1
            del ids[old_cl_idx[j]]    
        
        
################################################################################

#
# Class definitions
#

################################################################################         
        
class ClusterInfo:
    """ Data structure for cluster information"""
    
    def __init__(self, path):
        f = h5py.File(path, 'r')
        props = f['props']    
        self.idx = np.asarray(props['idx'], dtype=int)
        self.size = np.asarray(props['size'], dtype=int)
        self.comx = np.asarray(props['comx'], dtype=float)
        self.comy = np.asarray(props['comy'], dtype=float)
        self.rgysq = np.asarray(props['rgysq'], dtype=float)
        self.pl = np.asarray(props['planarity'], dtype=float)
        self.st = np.asarray(props['straightness'], dtype=float)
        self.sw = np.asarray(props['swirliness'], dtype=float)
        self.ens = np.asarray(props['enstrophy'], dtype=float) 
        fil_grp = props['filament_list']
        self.fils = []
        for el in np.arange(len(self.size)):
            fil_list = np.asarray(fil_grp[str(el)], dtype=int)
            self.fils.append(fil_list)
       
        f.close()
        
################################################################################ 

def main():
    
    ## Argument parsing (command line options)
    parser = argparse.ArgumentParser()
    parser.add_argument("initfile", help="Folder for initialization of data")
    parser.add_argument("datafile", help="Folder containing data")
    args = parser.parse_args()
    
    ## Load saved preliminary simulation data into relevant variables
    loadSimData(args.initfile+'/density_0.08/kappa_5.0/fp_0.24/init_info.txt')
    datafile = args.datafile+"/CLUSTER"
    
    ## Time averages
    cnt_frame = -1      # Counter for number of total time frames (-1, because first frame doesn't count)
    st_avg = 0          # Straightness / Asphericity
    ftr_avg = 0         # Translational force
    frot_avg = 0        # Rotational force
    ens_avg = 0         # Enstrophy
    rgysq_avg = 0       # Radius of gyration squared
    size_avg = 0        # Size of cluster
    
    lifetime = np.zeros(int(totalStep-ti), dtype=float)     # Array keeping track of histogram of lifetimes
    ids = {}            # Dictionary keeping track of number of times the id'ed cluster has been living
    old_cl_idx = []     # List keeping track of cluster ids from the previous time frame
    
    first_time = True
    
    size_rgy = np.zeros(M+1,dtype=float)      # Distribution of radius of gyration and cluster size
    size_st = np.zeros(M+1,dtype=float)       # Distribution of straigtness and cluster size
    size_ftr = np.zeros(M+1,dtype=float)      # Distribution of translational force and cluster size
    size_frot = np.zeros(M+1,dtype=float)     # Distribution of rotational force and cluster size
    size_ens = np.zeros(M+1,dtype=float)      # Distribution of enstrophy and cluster size    
    #size_life = np.zeros(M+1, dtype=float)   # Distribution of cluster size with cluster lifetime
    cnt_size = np.zeros(M+1,dtype=float)      # Number of times the cluster size is observed for normalization
    
    ## Time frame loop
    for frame in np.arange(int(ti),int(totalStep),int(nsamp)):
        
        print frame , ' of ', totalStep
        
        cnt_frame += 1
    
        ## Load cluster sizes
        path = datafile + '/cluster_' + str(frame) + '.hdf5'
        
        ## If the file doesn't exist, break the for loop
        if os.path.exists(path):
            
            ## Load the data
            clinfo = ClusterInfo(path)
            
            ## If it's the first time frame, save the cluster indices
            if first_time:
                first_time = False
                
                ## Save the old cluster ids
                old_cl_idx = np.asarray(clinfo.idx)
                
                ## Initiate the living time counters for each observed id in the first frame
                for j in range(len(clinfo.idx)):
                    ids[clinfo.idx[j]] = 1 
                    
                ## Go on to the second time step without going any further
                continue

            ## Cluster indices as array
            cl_idx = np.asarray(clinfo.idx)
            
            ## Find the cluster lifetime
            cluster_lifetime(cl_idx, old_cl_idx, lifetime, ids)    
                        
            ## Update the ids for the newly generated clusters
            old_cl_idx = np.asarray(clinfo.idx)
            for j in range(len(clinfo.idx)):
                key = clinfo.idx[j]
                ## Add new living time counters for the newly created clusters with unique ids
                if key not in ids:
                    ids[key] = 1
            
            ## Calculate the radius of gyration - cluster size distribution
            for j, siz in enumerate(clinfo.size):
                size_rgy[siz] += np.sqrt(clinfo.rgysq[j])
                size_st[siz] += clinfo.pl[j]
                size_ftr[siz] += clinfo.st[j]
                size_frot[siz] += clinfo.sw[j]
                size_ens[siz] += clinfo.ens[j]
                cnt_size[siz] += 1
            
            ## Average of observables over the total number of clusters with size larger than 10 in a time frame
            if len(clinfo.pl) != 0:
                st_avg += np.sum(clinfo.pl)/len(clinfo.pl)
                ftr_avg += np.sum(clinfo.st)/len(clinfo.st)
                frot_avg += np.sum(clinfo.sw)/len(clinfo.sw)
                ens_avg += np.sum(clinfo.ens)/len(clinfo.ens)
                rgysq_avg += np.sum(clinfo.rgysq)/len(clinfo.rgysq)
                size_avg += np.sum(clinfo.size)/len(clinfo.size)
#                for j, life in enumerate(lifetime):
#                    lifetime[j] += life/len(clinfo.st)
    
        ## Break the loop if the file cannot be found    
        else:
            break
    
    
    ## Average of observables over the total number of time frames    
    st_avg /= cnt_frame
    ftr_avg /= cnt_frame
    frot_avg /= cnt_frame
    ens_avg /= cnt_frame
    rgysq_avg /= cnt_frame
    size_avg /= cnt_frame
    
    
    ## Average of radius of gyration cluster size distribution
    for j in range(int(M+1)):
        if cnt_size[j] != 0:
            size_rgy[j] /= cnt_size[j]
            size_st[j] /= cnt_size[j]
            size_ftr[j] /= cnt_size[j]
            size_frot[j] /= cnt_size[j]
            size_ens[j] /= cnt_size[j]
    
    
    ## Calculate the average lifetime of a cluster
    lifetime_avg = 0
    for life, count_of_life in enumerate(lifetime):
        lifetime_avg += life*count_of_life
    lifetime_avg /= np.sum(lifetime)    
    
    ## Write data
    path = datafile + '/avg_straightness.txt'
    fl = open(path, 'w')
    fl.write(str(st_avg))
    fl.close()
    
    path = datafile + '/avg_translational.txt'
    fl = open(path, 'w')
    fl.write(str(ftr_avg))
    fl.close()
    
    path = datafile + '/avg_rotational.txt'
    fl = open(path, 'w')
    fl.write(str(frot_avg))
    fl.close()
    
    path = datafile + '/avg_enstrophy.txt'
    fl = open(path, 'w')
    fl.write(str(ens_avg))
    fl.close()
    
    path = datafile + '/avg_rgysq.txt'
    fl = open(path, 'w')
    fl.write(str(np.sqrt(rgysq_avg)))
    fl.close()
    
    path = datafile + '/avg_size.txt'
    fl = open(path, 'w')
    fl.write(str(size_avg))
    fl.close()
    
    path = datafile + '/avg_lifetime.txt'
    fl = open(path, 'w')
    fl.write(str(lifetime_avg))
    fl.close()
    
    path = datafile + '/lifetime.txt'
    fl = open(path, 'w')
    for el in lifetime:
        fl.write(str(el) + '\n')
    fl.close()
    
    path = datafile + '/size_rgy.txt'
    fl = open(path, 'w')
    for j in range(int(M+1)):
        fl.write(str(size_rgy[j]) + '\n')
    fl.close()

    path = datafile + '/size_st.txt'
    fl = open(path, 'w')
    for j in range(int(M+1)):
        fl.write(str(size_st[j]) + '\n')
    fl.close()
    
    path = datafile + '/size_ftr.txt'
    fl = open(path, 'w')
    for j in range(int(M+1)):
        fl.write(str(size_ftr[j]) + '\n')
    fl.close()

    path = datafile + '/size_frot.txt'
    fl = open(path, 'w')
    for j in range(int(M+1)):
        fl.write(str(size_frot[j]) + '\n')
    fl.close()

    path = datafile + '/size_ens.txt'
    fl = open(path, 'w')
    for j in range(int(M+1)):
        fl.write(str(size_ens[j]) + '\n')
    fl.close()

    
if __name__ == '__main__':
    main()

