
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

def compute_static_structure(charfile, headerfile, nfil, nskip, box):
    """ compute the static structure factor"""
    
    ### open files for reading
    
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    
    ### get general information from header file
    
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
        
    ### skip the initial snapshots
    
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    
    ### allocate arrays

    N = int(box/2.0)
    delta_k = 2*np.pi/box
    
    kx = ky = np.arange(0,N)*delta_k
    k = np.transpose(np.vstack((kx,ky)))

    term_to_average = np.zeros((N))
    S = np.zeros((N))
    
    x = np.zeros((natoms))
    y = np.zeros((natoms))
            
    ### read all the coordinates into a coordinate array

    for step in range(nsteps):
        
        ### print progress to screen
        
        print 'Current Step / All Steps', step, '/', nsteps
        
        ### read in coordinates
        
        xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
        x = xs*lx
        y = ys*ly
                        
        ### take the dot product k*r
        
        r = np.vstack((x,y))        
        dot_product = np.dot(k,r)
        term_to_average += np.sum(np.cos(dot_product),axis=1)**2 + np.sum(np.sin(dot_product),axis=1)**2
        
    ### normalize the average
        
    for j in range(N):
        S[j] = term_to_average[j]/natoms/nsteps

    ### close the input files

    hfile.close()
    ifile.close()
    
    return N, S

##################################################################

def main():
    
    ### read parameters from input file
    
    charfile, headerfile, ofname, nfil, nskip, box = read_settings()
    
    ### loop over the entire trajectory and compute static structure factor
    N, S = compute_static_structure(charfile, headerfile, nfil, nskip, box)
    delta_k = 2*np.pi/box
    
    ### write results to files
    
    # create output folder
    os.system('mkdir -p ' + ofname)
    
    # write the wavevector and static structure factor data
    ofile = open(ofname + '/static_structure_factor.data', 'w')
    for j in range(N):
        k = delta_k*j
        ofile.write(str(k) + '\t\t' + str(S[j]) + '\n')
        
    # close the file
    ofile.close()
        
    
##################################################################

if __name__ == '__main__':
    main()
    
##################################################################