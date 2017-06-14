#!/usr/local/bin/python2.7

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import codecs
import read_char
import os
import rdf_siteWrapper

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
    # number of bins
    line = ifile.readline()
    line = line.split()
    nbin = int(line[-1])
    # maximum cutoff
    line = ifile.readline()
    line = line.split()
    rmax = float(line[-1])
    # close file
    ifile.close()
    # return input values
    return charfile, headerfile, ofname, nfil, nskip, nbin, rmax

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

def compute_com_single(x,y,lx,ly):
    """ compute the absolute of the spiral number of a filament,
        correct pbc on the fly"""
    nbeads = len(x)
    # correct pbcs
    for i in range(1,nbeads):
        dx = x[i] - x[i-1]
        dy = y[i] - y[i-1]
        x[i] = x[i-1] + neigh_min(dx,lx)
        y[i] = y[i-1] + neigh_min(dy,ly)
    # compute comx and comy
    comx = np.average(x)
    comy = np.average(y)
    # put comx and comy into simulation box
    comx /= lx
    comx = comx - math.floor(comx)
    comx *= lx
    comy /= ly
    comy = comy - math.floor(comy)
    comy *= ly
    return comx, comy

############################################################################

def gen_linked_list(x,y,lx,ly,dcrit,nmol):
    """ generate a linked list"""
    # determine the number of cells in each direction
    nsegx = int(lx/dcrit)
    nsegy = int(ly/dcrit)

    # allocate head and llist
    ncells = nsegx*nsegy
    head = np.zeros((ncells), dtype = np.int32)
    llist = np.zeros((nmol), dtype = np.int32)

    # fill list and head
    for i in range(nmol):
        segx = int(x[i]/lx*nsegx)
        segy = int(y[i]/ly*nsegy)
        cell = segx*nsegy + segy
        llist[i] = head[cell]
        head[cell] = i

    return nsegx,nsegy,head,llist

##################################################################

def compute_rdf_single(x, y, nbin, rmax, lx, ly, nmol):
    """ compute the radial distribution function of the center of mass"""
    ### allocate array to store the results
    rdf = np.zeros((nbin), dtype = np.float64)
    ### generate a linked list of the com
    nsegx, nsegy, head, llist = gen_linked_list(x,y,lx,ly,rmax, nmol)
    ### loop over the linked list
    mol = np.arange(nmol, dtype = np.int32)
    rdf_site = rdf_siteWrapper.Rdf_site()
    rdf_site.compute(nsegx, nsegy, nmol, head, llist, mol, x, y, rdf, lx, ly, rmax, nbin)
    """
    for i in range(nsegx):
        for j in range(nsegy):
            # store header of current cell
            sv1 = head[i*nsegy + j]
            # loop over neighboring cells
            for a in range(3):
                i2 = (i-1+a)%nsegx
                for b in range(3):
                    j2 = (j-1+b)%nsegy
                    # store header of neighbor cell
                    sv2 = head[i2*nsegy + j2]
                    
                    # restore head values at for each new cell
                    val1 = sv1
                    val2 = sv2
                    while val1 != 0:
                        x1 = x[val1]/lx
                        y1 = y[val1]/ly
                        while val2 != 0:
                            if val1 != val2:
                                x2 = x[val2]/lx
                                y2 = y[val2]/ly

                                dx = x2-x1
                                dx = dx - math.floor(dx + 0.5)
                                dx = dx*lx
                                dy = y2-y1
                                dy = dy - math.floor(dy + 0.5)
                                dy = dy*ly
                                rsq = dx**2 + dy**2
                                if rsq < rmax**2:
                                    r = np.sqrt(rsq)
                                    seg = int(r/rmax*nbin)
                                    rdf[seg] += 1
                            val2 = llist[val2]
                        val1 = llist[val1]
                        val2 = sv2
    """
    return rdf

##################################################################

def compute_rdf_com(charfile, headerfile, nfil, nskip, nbin, rmax):
    """ compute the radial distribution function
        of the center of mass of all filaments"""
    ### open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
    ### skip the first snapshots
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - nskip
    ### allocate arrays, including output: radial distribution function
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    # center of mass coordinates
    comx = np.zeros((nmol), dtype = np.float64)
    comy = np.zeros((nmol), dtype = np.float64)
    # time
    time = np.zeros((nsteps))
    # rdf
    rdf = np.zeros((nbin, nsteps), dtype = np.float64)
    # distance array
    r = np.linspace(0.,rmax,nbin+1)
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
        ### compute com for each filament, correct pbc on the fly
        for i in range(nmol):
            comx[i], comy[i] = compute_com_single(x[i*nfil:(i+1)*nfil-1], y[i*nfil:(i+1)*nfil-1], lx, ly)
        ### compute the radial distribution function for the current snapshot
        rdfi = compute_rdf_single(comx, comy, nbin, rmax, lx, ly, nmol)
        ### copy rdfi to rdf
        for i in range(nbin):
            rdf[i,step] = rdfi[i]
    ### close the input files
    hfile.close()
    ifile.close()
    ### normalize the rdf
    V = np.zeros(nbin)
    for i in range(nbin):
        V[i] = np.pi*(r[i+1]**2 - r[i]**2)
    rho_av = nmol/lx/ly
    for i in range(nsteps):
        rdf[:,i] *= 1./V/rho_av/nmol
    return nsteps, time, r, rdf


##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip, nbin, rmax = read_settings()
    ### loop over the entire trajectory and compute the center of mass, correct for pbc
    nsteps, time, r, rdf = compute_rdf_com(charfile, headerfile, nfil, nskip, nbin, rmax)
    ### compute the average radial distributation function
    rdf_av = np.average(rdf, axis = 1)

    ### write results to file and generage a plot
    # generate folder structure
    os.system('mkdir ' + ofname)
    os.system('mkdir ' + ofname + '/time')
    # write data of the averaged rdf
    ofile = open(ofname + '/rdf_av.data', 'w')
    ofile.write('Averaged RDF over entire time\n\n')
    ofile.write('r_min\tr_max\trdf\n')
    for i in range(nbin):
        ofile.write(str(r[i]) + '\t' + str(r[i+1]) + '\t' + str(rdf_av[i]) + '\n')
    ofile.close()
    # gen figure of the averaged rdf
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(0.5*(r[:-1]+r[1:]), rdf_av)
    ax.set_xlabel(r'r [$\sigma]')
    ax.set_ylabel(r'g(r)')
    ax.set_title('Radial Distribution Function')
    plt.savefig(ofname + '/rdf_av.png')
    plt.close()
    # write data of the non-averaged rdf
    ofile = open(ofname + '/rdf_time.data', 'w')
    ofile.write('RDF of individual snapshots\n\n')
    ofile.write('nbin = ' + str(nbin) + '\n')
    ofile.write('nsteps = ' + str(nsteps) + '\n')
    ofile.write('\t\tTimesteps\n')
    ofile.write('rmin\trmax\t')
    for i in range(nsteps):
        ofile.write(str(time[i]) + '\t')
    ofile.write('\n')
    for i in range(nbin):
        ofile.write(str(r[i]) + '\t' + str(r[i+1]) + '\t')
        for j in range(nsteps):
            ofile.write(str(rdf[i,j]) + '\t')
        ofile.write('\n')
    ofile.close()
    # gen figures of the non-averaged rdf
    for i in range(nsteps):
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(0.5*(r[:-1]+r[1:]), rdf[:,i], label = 'tstep = ' + str(time[i]))
        ax.legend(loc = 'lower right')
        ax.set_xlabel(r'r [$\sigma]')
        ax.set_ylabel(r'g(r)')
        ax.set_title('Radial Distribution Function')
        plt.savefig(ofname + '/time/rdf_' + str(time[i]) + '.png')
        plt.close()
    return
    
##################################################################

if __name__ == '__main__':
    main()
    
