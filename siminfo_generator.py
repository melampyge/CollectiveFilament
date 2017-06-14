# Load needed libraries and necessary files
import argparse
import numpy as np
#import matplotlib as mpl
#mpl.use('Agg')
#import matplotlib.pyplot as plt
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

######################################################
 
# Load initial simulation data      
def saveSimData(datafile):
       
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
        dtSamp, T, box_area, nt, body_length, Pe, persistence 

    dens = datafile.split('/')[-4]
    dens = dens.split('_')[-1]
    ka = datafile.split('/')[-3]
    ka = ka.split('_')[-1]    
    fmc = datafile.split('/')[-2]
    fmc = fmc.split('_')[-1]
    print datafile, dens, ka, fmc
    
#    print dens
#    if dens == '0.08':
#        Lx = 1.36e+03
#        Ly = 1.36e+03
#    elif dens == '0.2':
#        Lx = 8.66+02
#        Ly = 8.66+02
#    elif dens == '0.4':
#        Lx = 612.37
#        Ly = 612.37
#    
#    L = 306000
        
#    datafile = open(datafile,"w")
#    datafile.write( "dt = " + str(1e-3) + "\n" )
#    datafile.write( "ti = " + str(1e+7) + "\n" )
#    datafile.write( "Lx = " + str(Lx) + "\n" )
#    datafile.write( "Ly = " + str(Ly) + "\n" )
#    datafile.write( "totalStep = " + str(totalStep) + "\n" )
#    datafile.write( "nsamp = " + str(nsamp) + "\n" )
#    datafile.write( "nfil = " + str(N) + "\n" )
#    datafile.write( "L = " + str(L) + "\n" )
#    datafile.write( "B = " + str(B) + "\n" )
#    datafile.write( "kT = " + str(kT) + "\n" )
#    datafile.write( "Fmc = " + str(Fmc) + "\n" )
#    datafile.write( "Kbend = " + str(Kbend) + "\n" )

    datafile = open(datafile,"w")
    datafile.write( "dt = " + str(1e-3) + "\n" )
    datafile.write( "ti = " + str(10000000) + "\n" )
    datafile.write( "Lx = " + str(866.025) + "\n" )
    datafile.write( "Ly = " + str(866.025) + "\n" )
    datafile.write( "totalStep = " + str(99950000) + "\n" )
    datafile.write( "nsamp = " + str(50000) + "\n" )
    datafile.write( "nfil = " + str(51) + "\n" )
    datafile.write( "L = " + str(306000) + "\n" )
    datafile.write( "B = " + str(0.5) + "\n" )
    datafile.write( "kT = " + str(1) + "\n" )
    datafile.write( "Fmc = " + str(fmc) + "\n" )
    datafile.write( "Kbend = " + str(ka) + "\n" )
    
    return

######################################################

# Load initial simulation data 
def loadSimData(datafile):
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
    dtSamp, T, box_area, nt, body_length, Pe, persistence 
    
    datafile = open(datafile,"r")
    for line in datafile:
        A = line.split()
        if len(A) > 0:
            if A[0] == "variable":
                if A[1] == "bl":
                    B = float(A[-1])
                elif A[1] == "kT":
                    kT = float(A[-1])
                elif A[1] == "kappa":
                    Kbend = float(A[-1])
                elif A[1] == "fp":
                    Fmc = float(A[-1])
            if A[0] == "run":
                totalStep = float(A[-1])
            if A[-1] == "xhi":
                Lx = float(A[1]) - float(A[0])
            if A[-1] == "yhi":
                Ly = float(A[1]) - float(A[0])
            if A[0] == "nfil":
                N = float(A[-1]) 
            if A[-1] == "atoms":
                L = float(A[0])
            if A[0] == "dump":
                nsamp = float(A[-2])
                    
    return
                
######################################################    

def main(): 
    
    # Argument parsing (command line options)
    parser = argparse.ArgumentParser()
    parser.add_argument("datafile", help="/local/duman/SIMULATIONS/long_filaments") 
    parser.add_argument("datafile1", help="/input1.lammps")
    parser.add_argument("datafile2", help="/input.data")
    parser.add_argument("datafile3", help="/analyse_long_cluster.data")
    parser.add_argument("savefile", help="/init_info.txt")
    args = parser.parse_args()
    
    
    # Index the data  
#    density = [0.08, 0.2, 0.4]
#    kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 400.0]
#    fp = [0.0, 0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
#    density = [0.2]
#    kappa = [5.0, 20.0, 100.0, 250.0, 1600.0]
#    fp = [0.0003, 0.015, 0.075, 1.0]    
    density = [0.2]
    kappa = [300.0]    
    fp = [0.0, 0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0] 
    
    # Create the paths to the folders contanining data
    folders = []
    for d in density:
        for k in kappa:
            for f in fp:
                folders.append( args.datafile + '/density_' + str(d) + \
                '/kappa_' + str(k) + '/fp_' + str(f) )
    
    # Format the data from the folder names   
    files = []
    tmp = {}
    for folder in folders:
        
        tmp = {}
        
        # Densities
        dname = folder.split('/')[-3]
        ddict = dict(zip(dname.split('_')[::2],map(atof,dname.split('_')[1::2])))
        
        # Bending rigidities
        kname = folder.split('/')[-2]
        kdict = dict(zip(kname.split('_')[::2],map(atof,kname.split('_')[1::2])))
        
        # Propulsion forces
        fname = folder.split('/')[-1]
        fdict = dict(zip(fname.split('_')[::2],map(atof,fname.split('_')[1::2])))
    
        # Put the variables to an auxiliary dictionary
        tmp.update(ddict)
        tmp.update(kdict)
        tmp.update(fdict)
    
        # Add the folders to the dictionary
        tmp.update({'folder':folder})
        files.append(tmp)
            
    # Data containing folder, density, kappa, fp values in Pandas DataFrame format
    Data = pd.DataFrame(files,columns=files[0].keys())
    
    # For each folder
    for folder in Data['folder']:
        
#        # Load the data
#        if os.path.exists(folder+args.datafile1):
#            loadSimData(folder+args.datafile1)
#        if os.path.exists(folder+args.datafile2):
#            loadSimData(folder+args.datafile2)
#        if os.path.exists(folder+args.datafile3):
#            loadSimData(folder+args.datafile3)
        
        # save the data
        saveSimData(folder+args.savefile)
        


if __name__ == '__main__':
    main()



