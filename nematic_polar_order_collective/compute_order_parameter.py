
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import os
import codecs
import read_char

try:
    infilename = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '      paraneter file'
    exit()
    

##################################################################
#    global variables

# number of fractions into which to decompose the box
frac = [1000., 700., 500., 300., 200., 100., 70., 50., 30., 20., 10., 7., 5., 3., 2., 1.]
frac = np.array(frac)
nfrac = len(frac)

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
    # close file
    ifile.close()
    # return input values
    return charfile, headerfile, ofname, nfil, nskip

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

############################################################################

def compute_orientation(x,y,lx,ly,nfil):
    """ compute orientation of all beads from bond vectores"""
    # number of molecules
    natoms = len(x)
    nmol = natoms/nfil
    # allocate aray for results
    phi = np.zeros((natoms), dtype = np.float64)
    tx = np.zeros((natoms), dtype = np.float64)
    ty = np.zeros((natoms), dtype = np.float64)
    tx2 = np.zeros((natoms), dtype = np.float64)
    ty2 =  np.zeros((natoms), dtype = np.float64)
    # loop over all polymers
    k = 0
    for i in range(nmol):
        for j in range(nfil):
            if j == 0:
                x1 = x[k]
                y1 = y[k]
                x2 = x[k+1]
                y2 = y[k+1]
            elif j == nfil-1:
                x1 = x[k-1]
                y1 = y[k-1]
                x2 = x[k]
                y2 = y[k]
            else:
                x1 = x[k-1]
                y1 = y[k-1]
                x2 = x[k+1]
                y2 = y[k+1]
            # compute nearest neighbor
            dx = neigh_min(x2-x1,lx)
            dy = neigh_min(y2-y1,ly)
            # compute angle using atan2
            pi = math.atan2(dy,dx)
            pi2 = pi
            if pi2 < 0:
                pi2 = np.pi - pi2
            phi[k] = pi
            tx[k] = dx / np.sqrt(dx**2 + dy**2)
            ty[k] = dy / np.sqrt(dx**2 + dy**2)
            tx2[k] = math.cos(2*pi2)
            ty2[k] = math.sin(2*pi2)
            # increment k
            k = k + 1
    return phi, tx, ty, tx2, ty2

##################################################################

def compute_order(x,y,n, tx, ty):
    """ compute number fluctuations"""
    ### compute histogram of the number of bins in each cell
    hist, xedges, yedges = np.histogram2d(x,y,bins=n, range = [[0,1],[0,1]])
    histx, xedges, yedges = np.histogram2d(x,y,bins=n, range = [[0,1],[0,1]], weights = tx)
    histy, xedges, yedges = np.histogram2d(x,y,bins=n, range = [[0,1],[0,1]], weights = ty)
    ### compute the order parameter
    k = 0
    order_parameter = 0
    for i in range(int(n)):
        for j in range(int(n)):
            h = hist[i,j]
            hx = histx[i,j]
            hy = histy[i,j]
            if hist[i,j] > 0:
                k = k  + 1
                order_parameter += ((hx)**2 + (hy)**2)/h**2 
    order_parameter /= k
    return order_parameter

##################################################################

def analyse_data(charfile, headerfile, nfil, nskip):
    """ analyse the simulation data"""
    ### open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
    ### skip the first snapshots
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    ### generate output array
    polar_order = np.zeros((nsteps,nfrac))
    nematic_order = np.zeros((nsteps, nfrac))
    t = np.zeros((nsteps))
    ### loop over the file
    for step in range(nsteps):
        ### print progress to screen
        print 'Current Step / All Steps', step, '/', nsteps
        ### read in coordinates
        xs, ys, lx, ly, tstep, natoms = read_char.read_snapshot(hfile, ifile)
        ### compute tangent vectors
        phi, tx, ty, tx2, ty2 = compute_orientation(xs,ys,1.0,1.0,nfil)
        # compute the number fluctuations for all bin widths
        #   and add results to output array
        t[step] = tstep
        for j in range(nfrac):
            Pj = compute_order(xs,ys,frac[j], tx, ty)
            polar_order[step,j] = Pj
            Sj = compute_order(xs,ys,frac[j], tx2, ty2)
            nematic_order[step,j] = Sj
    ### close input files
    ifile.close()
    hfile.close()
    ### compute the time averages
    P_av = np.average(polar_order, axis = 0)
    P_std = np.average(polar_order, axis = 0)/np.sqrt(nsteps)
    S_av = np.average(nematic_order, axis = 0)
    S_std = np.average(nematic_order, axis = 0)/np.sqrt(nsteps)
    ### compute an array with the average size of the box in each bin
    w_bins = lx / frac

    return w_bins, P_av, P_std, S_av, S_std, polar_order, nematic_order, t

##################################################################

def write_output(w, P, Pstd, S, Sstd, polar_order, nematic_order, t, ofname):
    """ create figures and tables of the simulation results"""
    ### generate output folder
    os.system('mkdir ' + ofname)
    ### generate lin-log plot of the avergaed order parameters
    # polar order
    fig = plt.figure()
    ax = plt.subplot(111)
    # add measured values
    ax.errorbar(w, P, yerr = Pstd, ls = '', marker = 'o')
    # modify the axes
    ax.set_xscale('log')
    ax.set_xlabel(r'$w [\sigma]$')
    ax.set_ylabel(r'$P$')
    # save and show the results
    plt.savefig(ofname + '/polar_order.png')
    plt.close()
    # nematic order
    fig = plt.figure()
    ax = plt.subplot(111)
    # add measured values
    ax.errorbar(w, S, yerr = Sstd, ls = '', marker = 'o')
    # modify the axes
    ax.set_xscale('log')
    ax.set_xlabel(r'$w [\sigma]$')
    ax.set_ylabel(r'$S$')
    # save and show the results
    plt.savefig(ofname + '/nematic_order.png')
    plt.close()
    ### generate plots of the global order parameter over time
    fig = plt.figure()
    ax = plt.subplot(111)
    # add the data
    ax.plot(t,polar_order[:,-1], label = 'polar order')
    ax.plot(t,nematic_order[:,-1], label = 'nematic order')
    # add the legend
    legend = ax.legend()
    # add the labels
    ax.set_xlabel(r'$t$ [timesteps]')
    ax.set_ylabel(r'order parameter')
    # save and close
    plt.savefig(ofname + '/time_evolution.png')
    plt.close()
    ### generate a table with the results
    # polar order parameter
    ofile = open(ofname + '/polar_order_parameter.data', 'w')
    ofile.write('Polar Order Parameter\n\n')
    ofile.write('\tbin_width\n')
    ofile.write('t')
    for j in range(len(w)):
        ofile.write('\t' + str(w[j]))
    ofile.write('\n')
    for i in range(len(t)):
        ofile.write(str(t[i]))
        for j in range(len(w)):
            ofile.write('\t' + str(polar_order[i,j]))
        ofile.write('\n')
    ofile.close()
    # nematic order parameter
    ofile = open(ofname + '/nematic_order_parameter.data', 'w')
    ofile.write('Nematic Order Parameter\n\n')
    ofile.write('\tbin_width\n')
    ofile.write('t')
    for j in range(len(w)):
        ofile.write('\t' + str(w[j]))
    ofile.write('\n')
    for i in range(len(t)):
        ofile.write(str(t[i]))
        for j in range(len(w)):
            ofile.write('\t' + str(nematic_order[i,j]))
        ofile.write('\n')
    ofile.close()
    return

##################################################################

def main():
    """ main function, called when script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip = read_settings()
    ### loop over trajectory and analyse the data
    w_bins, P_av, P_std, S_av, S_std, polar_order, nematic_order, t = analyse_data(charfile, headerfile, nfil, nskip)
    ### write results to files
    write_output(w_bins, P_av, P_std, S_av, S_std, polar_order, nematic_order, t, ofname)
      
    return

##################################################################

if __name__ == '__main__':
    main()
