
##################################################################

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
    # number of msr values to compute
    line = ifile.readline()
    line = line.split()
    nmsr = int(line[-1])
    # maximum value of dt compared to entire simulation time
    line = ifile.readline()
    line = line.split()
    limit = float(line[-1])
    # read data every that many steps
    line = ifile.readline()
    line = line.split()
    nevery = int(line[-1])
    # close file
    ifile.close()
    # return input values
    return charfile, headerfile, ofname, nfil, nskip, nmsr, limit, nevery

##################################################################

def neigh_min(dx,lx):
    ### compute minimum distance
    dx1 = dx + lx
    dx2 = dx - lx
    if dx**2 < dx1**2 and dx**2 < dx2**2:
        return dx
    if dx1**2 < dx2**2:
        return dx1
    return dx2

##################################################################

def compute_ee_single(x,y,lx,ly):
    """ compute the absolute of the spiral number of a filament,
        correct pbc on the fly"""
    nbeads = len(x)
    # correct pbcs
    for i in range(1,nbeads):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        x[i] = x[i-1] + neigh_min(dx,lx)
        y[i] = y[i-1] + neigh_min(dy,ly)
    # compute ex and ey
    ex = x[-1] - x[0]
    ey = y[-1] - y[0]
    return ex, ey

##################################################################

def get_ee(charfile, headerfile, nfil, nskip):
    """ extract the end-to-end vector of all filaments"""
    ### open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
    ### skip the initial snapshots
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    ### allocate arrays, including output: spiral histogram + spiral counter
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    # center of mass coordinates
    ex = np.zeros((nmol, nsteps), dtype = np.float64)
    ey = np.zeros((nmol, nsteps), dtype = np.float64)
    # time
    time = np.zeros((nsteps), dtype = np.int32)
    ### loop over all timesteps
    for step in range(nsteps):
        ### print progress to screen
        print 'Current Step / All Steps', step, '/', nsteps
        ### read in coordinates
        xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
        x = xs*lx
        y = ys*ly
        ### store time
        time[step] = tstep
        ### extract differences in endpoint coordinates for each filament, correct pbc on the fly
        for i in range(nmol):
            ex[i,step], ey[i,step] = compute_ee_single(x[i*nfil:(i+1)*nfil-1], y[i*nfil:(i+1)*nfil-1], lx, ly)
    ### close the input files
    hfile.close()
    ifile.close()
    return time, ex, ey

##################################################################

def compute_angle(ex,ey):
    """ compute phi"""
    ### allocate phi array
    nmol = len(ex)
    nsteps = len(ex[0])
    phi = np.zeros((nmol, nsteps), dtype = np.float64)
    dp = np.zeros((nsteps), dtype = np.float64)
    ### compute angle and correct angle on the fly
    performance_tools = performance_toolsWrapper.Performance_tools()
    for i in range(nmol):
        performance_tools.compute_angle(ex[i],ey[i],phi[i],dp,nsteps)
    return phi

##################################################################

def compute_msr(time, phi, nmsr, limit, nevery):
    """ compute and average the MSR"""
    ### compute required values
    nmol = len(phi)
    nsteps = len(phi[0])
    limit = int(limit*nsteps)
    ### allocate output array to store the results
    t_msr = np.zeros((nmsr), dtype = np.int32)
    msr = np.zeros((nmsr), dtype = np.float64)
    msr_tmp = np.zeros((nmsr), dtype = np.float64)
    counter = np.zeros((nmsr), dtype = np.int32)
    logval = np.zeros((nmsr), dtype = np.int32)
    ### loop over all molecules, submit MSD computation, read and average
    ###   the results
    performance_tools = performance_toolsWrapper.Performance_tools()
    for i in range(nmol):
        # submit MSR computation
        performance_tools.compute_MSR(t_msr, time, phi[i], msr_tmp, counter, logval, nsteps, nmsr, limit, nevery)
        # average the results using the on-line algorithm
        delta = msr_tmp - msr
        msr += delta/(i+1)
    return t_msr, msr

##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip, nmsr, limit, nevery = read_settings()
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    time, ex, ey = get_ee(charfile, headerfile, nfil, nskip)
    ### compute local angle phi based on overall rotation
    phi = compute_angle(ex,ey)
    ### start msr computations
    t_msr, msr = compute_msr(time, phi, nmsr, limit, nevery)
    ### remove temporary files
    os.system('rm ee_*.data')
    ### write results to file and generate a plot
    # output directory
    os.system('mkdir ' + ofname)
    # write table
    ofile = open(ofname + '/msr.data', 'w')
    ofile.write('Averaged MSRD over all molecules\n\n')
    ofile.write('Timestep\tMSR\n')
    for i in range(len(t_msr)):
        ofile.write(str(t_msr[i]) + '\t' + str(msr[i]) + '\n')
    ofile.close()
    os.system('mkdir ' + ofname)
#    fig = plt.figure()
#    ax = plt.subplot(111)
#    ax.loglog(t_msr,msr)
#    ax.set_xlabel('t [steps]')
#    ax.set_ylabel(r'msr')
#    ax.set_title('MSR')
#    plt.savefig(ofname + '/msr.png')
#    plt.close()
    
##################################################################

if __name__ == '__main__':
    main()
    
