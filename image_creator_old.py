
# Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os.path
import glob
import pandas as pd
from string import atof
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib._png import read_png
import colorsys 
import sys

#
# Function definitions
#


# Load initial simulation data      
def loadSimData(datafile):
    
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
    body_length = B*N
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    flexure = Pe/persistence


#
# Class definitions
#

# Particle data
class Particles:
    
    def __init__(self, path):
        file = np.transpose(np.loadtxt(path, dtype=float))
        self.xi = file[0]/B                 # Image particle positions in x
        self.yi = file[1]/B                 # Image particle positions in y 
        self.phi = file[2]                  # Bead orientation 
        self.cidx = file[3]                 # Cluster index        
    

def main(): 
    
    # Argument parsing (command line options)
    parser = argparse.ArgumentParser()
    parser.add_argument("datafile", help="Folder containing data")
    parser.add_argument("savefile", help="Folder in which data should be saved")
    args = parser.parse_args()
    
    
    # Index the data
    #density = [0.08, 0.2, 0.4, 0.6, 0.8]
    density = [0.08, 0.2, 0.4]
    kappa = [2.5, 5.0, 25.0, 62.5, 125.0, 200.0]
    fp = [0.0, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
    # Format the data from the folder names   
    folders = []
    for d in density:
        for k in kappa:
            for f in fp:
#                folders.append( args.datafile + '/job_d_' + str(d) + '_k_' + str(k) + '_f_' + str(f) )
                folders.append( args.datafile + '/density_' + str(d) + \
                '/kappa_' + str(k) + '/fp_' + str(f) )
    
    
    # Format the data            
    files = []
    tmp = {}
    for folder in folders:
        
        tmp = {}
        dname = folder.split('/')[-3]
        ddict = dict(zip(dname.split('_')[::2],map(atof,dname.split('_')[1::2])))
        
        kname = folder.split('/')[-2]
        kdict = dict(zip(kname.split('_')[::2],map(atof,kname.split('_')[1::2])))
        
        fname = folder.split('/')[-1]
        fdict = dict(zip(fname.split('_')[::2],map(atof,fname.split('_')[1::2])))
    
        tmp.update(ddict)
        tmp.update(kdict)
        tmp.update(fdict)
    
        
        tmp.update({'folder':folder})
        files.append(tmp)
            
    Data = pd.DataFrame(files,columns=files[0].keys())

    # For each folder
    for folder in Data['folder']:
        
        # Load preliminary information about the simulation
        print folder
        loadSimData(folder+"/init_info.txt")
        
        # Load the data files for the first and last images
        path = folder + "/CLUSTER"
        if os.path.exists(path) == False:
            continue
        contents = []
        for content in os.listdir(path):
            if content[:5] == 'beads':
                contents.append(int(content.split('_')[-1].split('.')[0]))
        
        # Get the identity of the first and last images
        if np.size(contents) == 0:
            continue
        last = max(contents)
        first = min(contents)
        print first, '\t', last, '\n'
        
        last_imag_path = path + "/beads_" + str(last) + ".txt"
        #first_imag_path = path + "/beads_" + str(first) + ".txt"
        
#        last_imag_path = path + "/beads_99950000.txt"
        
        # Load plot properties
        downlim = -5
        uplim = max(Lx,Ly)+5
        tick_interval = int(uplim/5)

        # Load particle data
        p_last = Particles(last_imag_path)   
        #p_first = Particles(first_imag_path)
         
        #plt.figure()
        quant_steps = 2056
        norm = mpl.colors.Normalize(vmin=-np.pi, vmax=np.pi)
        
        plt.scatter(p_last.xi, p_last.yi, s=1, c=p_last.phi, cmap=plt.cm.get_cmap('hsv',quant_steps), 
                    edgecolors='None', alpha=0.7, vmin=-np.pi, vmax=np.pi, norm=norm)
        #ax.set_xlabel('x/b',fontsize=8)
        #ax.set_ylabel('y/b',fontsize=8)
        plt.xlim((downlim,uplim))
        plt.ylim((downlim,uplim))
        plt.xticks( np.arange(0,uplim,tick_interval) )
        plt.yticks( np.arange(0,uplim,tick_interval) )
        plt.tick_params(axis='both', which='major', labelsize=8)
        
        folder_detail = folder.split('/')
        save_path = args.savefile + "/" + folder_detail[-3] + "/" + folder_detail[-2] + "/" + folder_detail[-1] + "/img"
        if os.path.exists(save_path):
            plt.savefig(save_path+'/fin.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
        else:
            os.mkdir(save_path)
            plt.savefig(save_path+'/fin.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
        plt.clf()


if __name__ == '__main__':
    main()



