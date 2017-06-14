
import numpy as np
import sys
import math
import codecs
import read_char
import os

try:
    infilename = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '       parameter file'
    exit()

##################################################################

def read_settings():
    """ read in the settings from the parameter file"""
    
    # open file for reading
    ifile = open(infilename, 'r')
    
    # skip comment line
    ifile.readline()
    
    # char-file name
    line = ifile.readline()
    line = line.split()
    charfile = line[-1]
    
    # header-file name
    line = ifile.readline()
    line = line.split()
    headerfile = line[-1]
    
    # output file
    line = ifile.readline()
    line = line.split()
    ofname = line[-1]
    
    # length of the filaments
    line = ifile.readline()
    line = line.split()
    nfil = int(line[-1])
    
    # number of snapshots to skip
    line = ifile.readline()
    line = line.split()
    nskip = int(line[-1])
    
    # close file
    ifile.close()
    
    return charfile, headerfile, ofname, nfil, nskip

##################################################################

def neigh_min(dx,lx):
    """ compute minimum distance"""
    
    dx1 = dx + lx
    dx2 = dx - lx
    if dx**2 < dx1**2 and dx**2 < dx2**2:
        return dx
    if dx1**2 < dx2**2:
        return dx1
    return dx2

##################################################################

def compute_rgy_single(x,y,lx,ly):
    """ compute the gyration tensor of a single filament"""
      
    ### allocate arrays
      
    nbeads = len(x)
    x = np.array(x)
    y = np.array(y)
    rgy = np.zeros((3))
    
    ### correct pbcs
    
    for i in range(1,nbeads):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        x[i] = x[i-1] + neigh_min(dx,lx)
        y[i] = y[i-1] + neigh_min(dy,ly)
        
    ### compute comx and comy
        
    comx = np.average(x)
    comy = np.average(y)
    
    ### compute gyration tensor
    
    xd = x-comx
    yd = y-comy
    
    rgy[0] = np.mean(xd**2)
    rgy[1] = np.mean(xd*yd)
    rgy[2] = np.mean(yd**2)
    
    return rgy

##################################################################

def compute_rgy(charfile, headerfile, nfil, nskip):
    """ compute the gyration tensor"""
    
    ### open files for reading
    
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    
    ### get general information from header file
    
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
    
    ### skip initial snapshots
    
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    
    ### allocate arrays
    
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    time = np.zeros((nsteps))
    
    ### rgy
    
    rgy = np.zeros((3,nsteps))
    eig = np.zeros((2,nsteps))
    rgy_temp = np.zeros((3))
    
    ### loop over all timesteps
    
    for step in range(nsteps):
        
        ### print progress to screen
        
        print 'Current Step / All Steps', step, '/', nsteps
        
        ### read in coordinates
        
        xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
        x = xs*lx
        y = ys*ly
        
        time[step] = tstep
        rgy_avg = np.zeros((3))
        eig_squared_avg = np.zeros((2))
        
        ### calculate gyration tensor for each filament, correct pbc on the fly
        
        for i in range(nmol):
            rgy_temp = compute_rgy_single(x[i*nfil:(i+1)*nfil-1], y[i*nfil:(i+1)*nfil-1], lx, ly)
            rgy_avg += rgy_temp
            rgy_matrix = [ [rgy_temp[0], rgy_temp[1]], [rgy_temp[1], rgy_temp[2]] ]
            eig_squared_avg += np.linalg.eigvalsh(rgy_matrix)
            
        rgy[0][step] = rgy_avg[0]/nmol
        rgy[1][step] = rgy_avg[1]/nmol
        rgy[2][step] = rgy_avg[2]/nmol
        eig[0][step] = eig_squared_avg[0]/nmol
        eig[1][step] = eig_squared_avg[1]/nmol
    
    ### close the input files
        
    hfile.close()
    ifile.close()
    
    return time, rgy, eig

##################################################################

def main():
    
    ### read parameters from input file
    
    charfile, headerfile, ofname, nfil, nskip = read_settings()
    
    ### loop over the entire trajectory and compute the gyration tensor, correct for pbc
    
    time, rgy, eig = compute_rgy(charfile, headerfile, nfil, nskip)
    
    ### write data to file
    
    os.system('mkdir -p ' + ofname)
    ofile = open(ofname + '/rgyr.data', 'w')
    for i in range(len(time)):
        ofile.write(str(time[i]) + '\t\t' + str(rgy[0][i]) + '\t\t' + str(rgy[1][i]) \
            + '\t\t' + str(rgy[2][i]) + '\t\t' + str(eig[0][i]) + '\t\t' + str(eig[1][i]) + '\n')
    
    ofile.close()
 

    
##################################################################

if __name__ == '__main__':
    main()
    
##################################################################
