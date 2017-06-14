
import numpy as np
import sys
import math
import codecs
import read_char
import os
import performance_toolsWrapper

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
    
    # output folder
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
    
    # box bound
    line = ifile.readline()
    line = line.split()
    box = float(line[-1])
    
    # close file
    ifile.close()
    
    return charfile, headerfile, ofname, nfil, nskip, box

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

def analyse_static_structure(charfile, headerfile, nfil, nskip, box):
    """ load the position data at once"""
    
    ### open files for reading
    
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    
    ### get general information from header file
    
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
    
    ### skip the initial snapshots
    
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    
    nsteps = 51
    
    ### allocate arrays
    
    x = np.zeros((nsteps*natoms))
    y = np.zeros((nsteps*natoms))
            
    ### read all the coordinates into a coordinate array
    
    for step in range(nsteps):
        
        ### print progress to screen
        
        print 'READING - Current Step / All Steps', step, '/', nsteps
        
        ### read in coordinates
        
        xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
        for i in range(len(xs)):
            x[step*natoms+i] = xs[i]*lx
            y[step*natoms+i] = ys[i]*ly
         
    ### close the input files

    hfile.close()
    ifile.close()

    ### compute the structure factor

    N = int(box/2.0)
    S = np.zeros((N**2))

    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.compute_static_structure(S, x, y, nsteps, natoms, box)     
    
    return S

##################################################################

def main():
    
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip, box = read_settings()
    
    ### loop over the entire trajectory and count spirals
    S = analyse_static_structure(charfile, headerfile, nfil, nskip, box)
    N = int(box/2.0)
    delk = 2*np.pi/box
    
    ### write results to files
    
    # create output folder
    os.system('mkdir ' + ofname)
    
    # write the wavevector and static structure factor data
    ofile = open(ofname + '/static_structure_factor_numba.data', 'w')
    for nx in range(N):
        kx = delk*nx
        
        for ny in range(N):
            ky = delk*ny
            
            k = int(np.sqrt(kx**2 + ky**2))
            ofile.write(str(k) + '\t\t' + str(S[k]) + '\n')
        
    # close the file
    ofile.close()
        
    
##################################################################

if __name__ == '__main__':
    main()
    
##################################################################