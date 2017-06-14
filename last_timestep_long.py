# Load needed libraries and necessary files
import argparse
import numpy as np
import os
import pandas as pd
from string import atof
import fnmatch

#
# Function definitions
#
 
# Load initial simulation data      
def loadSimData(datafile):
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
    dtSamp, T, box_area, nt, body_length, Pe, persistence 

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


def main(): 
    
    ## Argument parsing (command line options)
    parser = argparse.ArgumentParser()
    parser.add_argument("datafile", help="Folder containing data")
    args = parser.parse_args()
    
#    ## Load saved preliminary simulation data into relevant variables
#    loadSimData(args.datafile+'/density_0.2/kappa_5.0/fp_0.24/init_info.txt')

    ## Index the data
    density = [0.2]
    kappa = [5.0, 20.0, 100.0, 250.0, 1600.0]
    fp = [0.0003, 0.015, 0.075]
    
    ## Index the folders
    folders = []
    for d in density:
        for k in kappa:
            for f in fp:
                folders.append( args.datafile + '/density_' + str(d) + \
                '/kappa_' + str(k) + '/fp_' + str(f) )
                
    last_dumps = {}         # key->folder name, value->last dump folder  
    finished_folders = {}   # key->folder name, value->1(finished),0(unfinished)      
    cnt = 0
    
    for folder in folders:  
        
        if os.path.exists(folder):
            
            print folder    
            
            ## Index dump folders and get the largest dump file 
            dump_numbers = []
            for file in os.listdir(folder):
                if fnmatch.fnmatch(file, '*.dump'):
                    dump_numbers.append(int(file.split('.')[0][3:]))
            if len(dump_numbers) != 0:
                max_dump_number = max(dump_numbers)   
            else:
                continue
            last_dump_file = folder + '/out' + str(max_dump_number) + '.dump'
            last_dumps[folder] = last_dump_file
            
            ## Check if the last dump folder containts the last time step
            os.system('grep 99950000 ' + last_dump_file + ' > tmp.txt')
            #os.system('tail -1 ' + last_dump_file + ' | awk \"{print $1}\" > tmp.txt')
            ifile = open('tmp.txt')
            line = ifile.readline()
            if len(line) > 0:
                finished_folders[folder] = 1
            else:
                finished_folders[folder] = 0
            print finished_folders[folder]
            ifile.close()
            os.system('rm tmp.txt')
    
            ## If the simulation is not finished, create a restart job
            if finished_folders[folder] == 0:
                cnt += 1
                os.system('cp ' + folder + '/job' + str(max_dump_number) + '.cmd ' + folder + '/job' + str(max_dump_number+1) + '.cmd')
                os.system('cp ' + folder + '/input' + str(max_dump_number) + '.lammps ' + folder + '/input' + str(max_dump_number+1) + '.lammps')
                os.system('sed -i \"s/backup' + str(2*max_dump_number-1) + '.rs/backup' + str(2*max_dump_number+1) + '.rs/g\" ' + folder + '/input' + str(max_dump_number+1) + '.lammps')
                os.system('sed -i \"s/backup' + str(2*max_dump_number) + '.rs/backup' + str(2*max_dump_number+2) + '.rs/g\" ' + folder + '/input' + str(max_dump_number+1) + '.lammps')
                os.system('sed -i \"s/out' + str(max_dump_number) + '.dump/out' + str(max_dump_number+1) + '.dump/g\" ' + folder + '/input' + str(max_dump_number+1) + '.lammps') 
                os.system('sed -i \"s/out' + str(max_dump_number) + '.rs/out' + str(max_dump_number+1) + '.rs/g\" ' + folder + '/input' + str(max_dump_number+1) + '.lammps')       
                os.system('sed -i \"s/out' + str(max_dump_number) + '.data/out' + str(max_dump_number+1) + '.data/g\" ' + folder + '/input' + str(max_dump_number+1) + '.lammps')             
                if max_dump_number == 1:
                    os.system('sed -i \"s/read_data/read_restart/g\" ' + folder + '/input' + str(max_dump_number+1) + '.lammps')
                    os.system('sed -i \"s/input.data/backup' + str(2*max_dump_number) + '.rs/g\" ' + folder + '/input' + str(max_dump_number+1) + '.lammps')            
                    os.system('sed -i \"s/100000000/50000000/g\" ' + folder + '/input' + str(max_dump_number+1) + '.lammps')           
                os.system('sed -i \"s/log' + str(max_dump_number) + '.lammps/log' + str(max_dump_number+1) + '.lammps/g\" ' + folder + '/job' + str(max_dump_number+1) + '.cmd')             
                os.system('sed -i \"s/input' + str(max_dump_number) + '.lammps/input' + str(max_dump_number+1) + '.lammps/g\" ' + folder + '/job' + str(max_dump_number+1) + '.cmd')   
                os.system('sed -i \"s/time=24:00:00/time=6:00:00/g\" ' + folder + '/job' + str(max_dump_number+1) + '.cmd')   
                
                opath = 'joblist_restart/job' + str(cnt) + '.sh'
                ofile = open(opath, 'w')
                ofile.write('#!/bin/bash\n\n')
                ofile.write('cd ' + folder + '\n')
                ofile.write('chmod 755 job' + str(max_dump_number+1) + '.cmd\n')
                ofile.write('sbatch job' + str(max_dump_number+1) + '.cmd\n')
                ofile.write('cd -\n')
                ofile.close()
    #            os.system('cp joblist_restart/job0.sh joblist_restart/job' + str(cnt) + '.sh')
    #            os.system('sed -i \"s/job/job' + str(max_dump_number+1) + '.cmd/g\" joblist_restart/job' + str(cnt) + '.sh')
    #            os.system('sed -i \"s/destination/' + folder + '/g\" joblist_restart/job' + str(cnt) + '.sh')
                
        else:
            continue
            
            
        

if __name__ == '__main__':
    main()



