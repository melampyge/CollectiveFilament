#!/usr/local/bin/python2.7

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mplcolors
import sys
import math
import os
import codecs
from scipy import optimize
import read_char

try:
    ifname = sys.argv[1]

except:
    print 'Usage: ' + sys.argv[0] + '     parameter file'
    exit()
    
##################################################################

def read_settings():
    """ read in the settings from the parameter file"""
    # open file for reading
    ifile = open(ifname, 'r')
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
    # length of the filaments
    line = ifile.readline()
    line = line.split()
    nfil = int(line[-1])
    # critical distance for mol/mol neighbor search
    line = ifile.readline()
    line = line.split()
    dcrit = float(line[-1])
    # critical overlap for mol/mol neibhor search
    line = ifile.readline()
    line = line.split()
    lcrit = float(line[-1])
    # critical angle for mol/mol neighbor search
    line = ifile.readline()
    line = line.split()
    pcrit = float(line[-1])
    pcrit *= np.pi/180.0
    # critical distance for cluster/cluster neighbor search
    # close file
    ifile.close()
    # return input values
    return charfile, headerfile, nfil, dcrit, lcrit, pcrit

############################################################################

def gen_mol_info(natoms, nfil):
    """ generate information about molecular ID"""
    mol = np.zeros((natoms), dtype = int)
    nmol = natoms/nfil
    k = 0
    for i in range(nmol):
        for j in range(nfil):
            mol[k] = i
            k = k + 1
    return mol,nmol

############################################################################

def nearest_neighbor(x1,x2,lx):
    """ compute vector to nearest neighbor"""
    dx1 = x1 - x2
    dx2 = x1 - x2 + lx
    dx3 = x1 - x2 - lx
    if dx1**2 < dx2**2 and dx1**2 < dx3**2:
        return dx1
    if dx2**2 < dx3**2:
        return dx2
    return dx3

############################################################################

def compute_orientation(x,y,lx,ly,nfil):
    """ compute orientation of all beads from bond vectores"""
    # number of molecules
    natoms = len(x)
    nmol = natoms/nfil
    # allocate aray for results
    phi = np.zeros((natoms))
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
            dx = nearest_neighbor(x1,x2,lx)
            dy = nearest_neighbor(y1,y2,ly)
            # compute angle using atan2
            pi = math.atan2(dy,dx)
            phi[k] = pi
            # increment k
            k = k + 1
    return phi

############################################################################

def gen_linked_list(x,y,lx,ly,dcrit):
    """ generate a linked list"""
    # determine the number of cells in each direction
    nsegx = int(lx/dcrit)
    nsegy = int(ly/dcrit)

    # allocate head and llist
    ncells = nsegx*nsegy
    natoms = len(x)
    head = np.zeros((ncells), dtype = int)
    llist = np.zeros((natoms), dtype = int)

    # fill list and head
    for i in range(natoms):
        segx = int(x[i]/lx*nsegx)
        segy = int(y[i]/ly*nsegy)
        cell = segx*nsegy + segy
        llist[i] = head[cell]
        head[cell] = i

    return nsegx,nsegy,head,llist

##################################################################

def fill_neigh_matrix(neighs,llist,head,nsegx,nsegy,x,y,phi,mol,lx,ly,dcrit,pcrit):
    """ count the number of neighboring beads of all molecules"""
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
                        p1 = phi[val1]
                        mol1 = mol[val1]
                        while val2 != 0:
                            if val1 != val2:
                                x2 = x[val2]/lx
                                y2 = y[val2]/ly
                                p2 = phi[val2]
                                mol2 = mol[val2]
                                if mol1!= mol2:
                                    dx = x2-x1
                                    dx = dx - math.floor(dx + 0.5)
                                    dx = dx*lx
                                    dy = y2-y1
                                    dy = dy - math.floor(dy + 0.5)
                                    dy = dy*ly
                                    dphi = p1 - p2
                                    if dx**2 + dy**2 < dcrit**2 and min([dphi**2, (dphi+2*np.pi)**2, (dphi-2*np.pi)**2])  < pcrit**2:
                                        neighs[mol1,mol2] += 1
                                        neighs[mol2,mol1] += 1
                            val2 = llist[val2]
                        val1 = llist[val1]
                        val2 = sv2
    return

##################################################################

def fill_neigh_cl_matrix((neighs,llist,head,nsegx,nsegy,x,y,phi,mol,lx,ly,dcrit,pcrit)):
    """ count the number of neighboring beads of all molecules"""
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
                        p1 = phi[val1]
                        mol1 = mol[val1]
                        while val2 != 0:
                            if val1 != val2:
                                x2 = x[val2]/lx
                                y2 = y[val2]/ly
                                p2 = phi[val2]
                                mol2 = mol[val2]
                                if mol1!= mol2:
                                    dx = x2-x1
                                    dx = dx - math.floor(dx + 0.5)
                                    dx = dx*lx
                                    dy = y2-y1
                                    dy = dy - math.floor(dy + 0.5)
                                    dy = dy*ly
                                    dphi = p1 - p2
                                    if dx**2 + dy**2 < dcrit**2 and min([dphi**2, (dphi+2*np.pi)**2, (dphi-2*np.pi)**2])  < pcrit**2:
                                        neighs[mol1,mol2] += 1
                                        neighs[mol2,mol1] += 1
                            val2 = llist[val2]
                        val1 = llist[val1]
                        val2 = sv2
    return

##################################################################

def recursion(neighs,cl,i,nmol,ncrit):
    """ recursiely find clusters"""
    for j in range(nmol):
        if neighs[i,j] > ncrit:
                if cl[j] == -1:
                    cl[j] = cl[i]
                    recursion(neighs,cl,j,nmol,ncrit)

##################################################################

def cluster_search(neighs,cl,lcrit,nfil,nmol):
    """ recursively search for clusters"""
    # define critical overlap
    ncrit = lcrit*nfil
    # loop over all colums of the matrix:
    for i in range(nmol):
        # assign a cluster id to the current molecule
        if cl[i] != -1:
            continue
        cl[i] = max(cl) + 1
        recursion(neighs,cl,i,nmol,ncrit)

    return
                    


##################################################################

def find_clusters(x,y,phi,mol,nmol,nfil,lx,ly,dcrit,lcrit,pcrit):
    """ search for clusters; two molecules are defined as being
        part of the same cluster if lcrit of their body length is
        within a distance of dcrit and when the difference in orientation
        is less than phi"""
    natoms = len(x)
    ### STEP ONE: COMBINE MOLECULES TO CLUSTERS BASED ON MOLECULAR NEIGHBORSHIP
    # allocate required arrays
    cl = np.zeros((nmol), dtype = int)
    cl -= 1
    clusters = np.zeros((natoms), dtype = int)
    neighs_mol = np.zeros((nmol,nmol), dtype = int)
    # generate a linked list
    nsegx, nsegy, head, llist = gen_linked_list(x,y,lx,ly,dcrit)
    # allocate array to store neighbors
    fill_neigh_matrix(neighs_mol,llist,head,nsegx,nsegy,x,y,phi,mol,lx,ly,dcrit,pcrit)
    # recursive search for clusters in neighbor matrix
    cluster_search(neighs_mol,cl,lcrit,nfil,nmol)
    # fill cluster results to per atom array
    k = 0
    for i in range(nmol):
        for j in range(nfil):
            clusters[k] = cl[i]
            k = k + 1
    ### return results
    return clusters

##################################################################

def gen_cluster_plot(x,y,clusters,tstep):
    """ generate a plot of the atoms color coded with the cluster
        to which they belong"""
    fig = plt.figure()
    ax = plt.subplot(111)
    #ax.scatter(x,y,s=5,c=clusters, linewidths = 0)
    nclusters = max(clusters)
    for i in range(nclusters):
        ax.plot(x[clusters==i],y[clusters==i],ls = '', marker = 'o', markeredgewidth = 0, markersize = 1, color = mplcolors.cnames.keys()[i%len(mplcolors.cnames)])
    ax.axis('equal')
    os.system('mkdir test')
    plt.savefig('test/cluster_' + str(tstep) + '.png')
    plt.close()
    return

##################################################################

def transform_cluster_data(cluster_atm, nfil):
    """ transform data: input (cluster_atm) is a per/atom array.
        Transform to array that contains all molecules of each
        cluster"""
    clusters = []
    # find the maximum number of clusters
    clmax = max(cluster_atm)
    # add emtpy lists to cluster
    for i in range(clmax+1):
        clusters.append([])
    # fill lists with mol IDs
    nmol = len(cluster_atm)/nfil
    for i in range(nmol):
        ci = cluster_atm[i*nfil]
        clusters[ci].append(i)
    return clusters

##################################################################

def recurse_cluster(h,h2,nx,ny,sx0,sy0):
    """ find the cluster recursively"""
    # write input file for cpp code
    ofile = open('h.data', 'w')
    ofile.write(str(nx) + ' ' + str(ny) + ' ' + str(sx0) + ' ' + str(sy0) + '\n')
    for i in range(nx):
        for j in range(ny):
            ofile.write(str(int(h[i,j])) + '\n')
    ofile.close()
    # start cpp code
    os.system('~/Applications/Scripts/clusters/recurse_cluster')
    # read in results
    ifile = open('h2.data', 'r')
    for i in range(nx):
        for j in range(ny):
            line = ifile.readline()
            line = line.split()
            h2[i,j] = int(line[0])
    ifile.close() 
 
    # cpp code does the following computation
    """
    for i in range(-1,2):
        for j in range(-1,2):
            sx = sx0 + i
            sy = sy0 + j
            if sx < 0 or sx >= nx:
                continue
            if sy < 0 or sy >= ny:
                continue
            if h[sx,sy] == 1 and h2[sx,sy] == 0:
                h2[sx,sy] = 1
                recurse_cluster(h,h2,nx,ny,sx,sy)
    """
    return

##################################################################

def explore_cluster(h,nx,ny,sx0,sy0):
    """ recursively search for an isolated cluster"""
    h2 = np.zeros((nx,ny), dtype = int)
    # set starting value for h2
    h2[sx0,sy0] = 1
    # fill h2 iteratively
    recurse_cluster(h,h2,nx,ny,sx0,sy0)
    return h2

##################################################################

def expand_cluster(x,y,lx,ly,clusters,i,nfil):
    """ expand the current cluster to periodic images"""
    nmol = len(clusters[i])
    xcl = np.zeros((9*nmol*nfil))
    ycl = np.zeros((9*nmol*nfil))
    # create copies of the atoms of all clusters
    l = 0
    for j in range(nmol):
        mi = clusters[i][j]
        for k in range(nfil):
            xcl[l] = x[mi*nfil + k]
            ycl[l] = y[mi*nfil + k]
            l = l + 1
    # add pbc copies
    for j in range(-1,2):
        for k in range(-1,2):
            if j == 0 and k == 0:
                continue
            for m in range(nmol*nfil):
                xcl[l] = xcl[m] + j*lx
                ycl[l] = ycl[m] + k*ly
                l = l + 1
    return xcl, ycl

##################################################################

def cluster_hist_1d(xcl,ycl,lx,ly,nmol,dcrit):
    """ try to detect cluster by 1d histograms"""
    # compute 1D histograms
    nxbins = int(3*lx/(dcrit))
    bx = 3*lx/nxbins
    hx, xedges = np.histogram(xcl, bins = nxbins, range = (-lx,2*lx))
    nybins = int(3*ly/(dcrit))
    by = 3*ly/nybins
    hy, yedges = np.histogram(ycl, bins = nybins, range = (-ly,2*ly))
    # check whether approach is feasible
    if 0 in hx and 0 in hy:  # 1d approach works
        # find boundaries
        xmin = -1
        xmax = -1
        ymin = -1
        ymax = -1
        jstop = -1
        for j in range(nxbins-1):
            if hx[j] == 0 and hx[j+1] > 0:
                xmin = xedges[j+1]
                jstop = j
                break
        for j in range(jstop,nxbins-1):
            if hx[j] > 0 and hx[j+1] == 0:
                xmax = xedges[j+1]
                break
        for j in range(nybins-1):
            if hy[j] == 0 and hy[j+1] > 0:
                ymin = yedges[j+1]
                jstop = j
                break
        for j in range(jstop,nybins-1):
            if hy[j] > 0 and hy[j+1] == 0:
                ymax = yedges[j+1]
                break
        # map all atom coordinates into relevant area
        for j in range(nmol*nfil):
            while xcl[j] < xmin:
                xcl[j] += lx
            while xcl[j] > xmax:
                xcl[j] -= lx
            while ycl[j] < ymin:
                ycl[j] += ly
            while ycl[j] > ymax:
                ycl[j] -= ly
        return 1
    return 0

#################################################################

def cluster_hist_2d(xcl,ycl,lx,ly,nmol,dcrit):
    """ try to cluster based on 2d histogramms"""
    nxbins = int(3*lx/(2*dcrit))
    bx = 3*lx/nxbins
    nybins = int(3*ly/(2*dcrit))
    by = 3*ly/nybins
    hxy, xedges, yedges = np.histogram2d(xcl, ycl, bins = [nxbins,nybins], range = [[-lx,2*lx], [-ly,2*ly]])
    # modify hxy histogram
    for i in range(nxbins):
        for j in range(nybins):
            if hxy[i,j] > 0:
                hxy[i,j] = 1
    # start from arbitrary point and fill cluster recursively
    segx = int((xcl[0] + lx)/bx)
    segy = int((ycl[0] + ly)/by)
    hxy2 =  explore_cluster(hxy,nxbins,nybins,segx,segy)
    # check whether cluster is isolated
    if np.sum(hxy2) < np.sum(hxy) - 0.5:
        for j in range(nmol*nfil):
            for k in range(9):
                xi = xcl[j + k*nmol*nfil]
                yi = ycl[j + k*nmol*nfil]
                # check whether xi and yi is cluster
                segx = int((xi + lx)/bx)
                segy = int((yi + ly)/by)
                if hxy2[segx,segy] == 1:
                    xcl[j] = xi
                    ycl[j] = yi
                    break
        return 1
    return 0

##################################################################

def correct_pbc_single(xcl,ycl,nmol,lx,ly,nfil):
    """ correct periodic boundary conditions of all molecules
        separately"""
    # loop over individual molecules in the cluster
    #   connect nearest neighbors
    l = 0
    for j in range(nmol):
        for k in range(nfil-1):
            x0 = xcl[l]
            y0 = ycl[l]
            x1 = xcl[l+1]
            y1 = ycl[l+1]
            dx = nearest_neighbor(x0,x1,lx)
            dy = nearest_neighbor(y1,y2,ly)
            xcl[l+1] = x0 + dx
            ycl[l+1] = y0 + dy
            l = l + 1
        l = l + 1
    # loop over all molecules in the cluster
    #   adjust com position
    for j in range(nmol):
        comx = np.average(xcl[j*nfil:(j+1)*nfil - 1])
        comy = np.average(ycl[j*nfil:(j+1)*nfil - 1])
        while comx < 0:
            comx += lx
            xcl[j*nfil:(j+1)*nfil - 1] += lx
        while comx > lx:
            comx -= lx
            xcl[j*nfil:(j+1)*nfil - 1] -= lx
        while comy < 0:
            comy += ly
            ycl[j*nfil:(j+1)*nfil - 1] += ly
        while comy > ly:
            comy -= ly
            ycl[j*nfil:(j+1)*nfil - 1] -= ly
    return

##################################################################

def adjust_com_cluster(xcl,ycl,lx,ly,nmol,nfil):
    """ move cluster such that com is in periodic box"""
    comx = np.average(xcl[0:nmol*nfil])
    comy = np.average(ycl[0:nmol*nfil])
    while comx < 0:
        comx += lx
        xcl += lx
    while comx > lx:
        comx -= lx
        xcl -= lx
    while comy < 0:
        comy += ly
        ycl += ly
    while comy > ly:
        comy -= ly
        ycl -= ly
    return

##################################################################

def correct_cluster_pbc(x,y,lx,ly,clusters,nfil,dcrit):
    """ create copy of x and y such that clusters are not separated
        by pbcs"""
    # number of clusters
    nclusters = len(clusters)
    # initialize output clusters and isolated array
    xcluster = np.copy(x)
    ycluster = np.copy(y)
    isolated = np.ones((nclusters), dtype = int)
    # loop over all clusters
    for i in range(nclusters):
        # allocate array to store cluster + copies
        nmol = len(clusters[i])
        # expand cluster to periodic images
        xcl, ycl = expand_cluster(x,y,lx,ly,clusters,i,nfil)
        # try simple approach with 1D histograms
        flag_1d = cluster_hist_1d(xcl,ycl,lx,ly,nmol,dcrit)
        # try approach with 2D histgrams if 1d failed
        if flag_1d == 0:
            flag_2d = cluster_hist_2d(xcl,ycl,lx,ly,nmol,dcrit)
        # move molecules if both approaches failed
        if flag_1d == 0 and flag_2d == 0:
            isolated[i] = 0
            correct_pbc_single(xcl,ycl,nmol,lx,ly,nfil)
        # adjust center of mass of entire cluster
        adjust_com_cluster(xcl,ycl,lx,ly,nmol,nfil)
        # copy coordinates to cluster array
        l = 0
        for j in range(nmol):
            mi = clusters[i][j]
            for k in range(nfil):
                xcluster[mi*nfil + k] = xcl[l]
                ycluster[mi*nfil + k] = ycl[l]
                l = l + 1
    return xcluster, ycluster, isolated

##################################################################

def make_rot_max(x,y,ncl,clusters,i,dfx_com,dfy_com):
    def rot_max(c):
        cxi = c[0]
        cyi = c[1]
        fr = 0
        for j in range(ncl):
            mi = clusters[i][j]
            for k in range(nfil - 1):
                x0 = x[mi*nfil + k]
                y0 = y[mi*nfil + k]
                x1 = x[mi*nfil + k + 1]
                y1 = y[mi*nfil + k + 1]
                dx = x0 - x1
                dy = y0 - y1
                # project bond force on translational force
                pr = dx*dfx_com + dy*dfy_com
                dx = dx*(1-pr)
                dy = dy*(1-pr)
                xc = 0.5*(x0 + x1)
                yc = 0.5*(y0 + y1)
                fr += ((xc-cxi)*dy - (yc-cyi)*dx)/np.sqrt((xc-cxi)**2 + (yc-cyi)**2)
        frinv = 1./math.fabs(fr)
        return frinv
    return rot_max

##################################################################

def characterize_cluster(x,y,phi,clusters,lx,ly,nfil):
    """ compute the external torque and net force that act on each cluster"""
    n = len(clusters)
    cl_mol = np.zeros((n), dtype = int)
    forcex = np.zeros(n)
    forcey = np.zeros(n)
    torque = np.zeros(n)
    comx = np.zeros(n)
    comy = np.zeros(n)
    rgyrsq = np.zeros(n)
    vx = np.zeros(n)
    vy = np.zeros(n)
    omega = np.zeros(n)
    ftotal = np.zeros(n)
    ftrans = np.zeros(n)
    frot = np.zeros(n)
    frot_max = np.zeros(n)
    xrot = np.zeros(n)
    yrot = np.zeros(n)
    etrans = np.zeros(n)
    erot = np.zeros(n)
    straightness = np.zeros(n)
    swirlicity = np.zeros(n)
    # loop over all clusters
    for i in range(n):
        # compute the center of mass of the cluster
        cxi = 0
        cyi = 0
        ncl = len(clusters[i])
        cl_mol[i] = ncl
        natoms = ncl*nfil
        for j in range(ncl):
            mi = clusters[i][j]
            for k in range(nfil):
                cxi += x[mi*nfil + k]
                cyi += y[mi*nfil + k]
        cxi /= natoms
        cyi /= natoms
        comx[i] = cxi
        comy[i] = cyi
        # compute the radius of gyration squared
        rgisq = 0
        for j in range(ncl):
            mi = clusters[i][j]
            for k in range(nfil): 
                dx = x[mi*nfil + k] - cxi
                dy = y[mi*nfil + k] - cyi
                rgisq += dx**2 + dy**2
        rgisq /= natoms
        rgyrsq[i] = rgisq
        # compute the net external force (equivalent to net velocity)
        fix = 0
        fiy = 0
        for j in range(ncl):
            mi = clusters[i][j]
            x0 = x[mi*nfil]
            y0 = y[mi*nfil]
            x1 = x[(mi+1)*nfil - 1]
            y1 = y[(mi+1)*nfil - 1]
            dx = x0 - x1
            dy = y0 - y1
            fix += dx
            fiy += dy
        forcex[i] = fix
        forcey[i] = fiy
        # compute the net torque (equivalent to net angular velocity)
        ti = 0
        for j in range(ncl):
            mi = clusters[i][j]
            for k in range(nfil - 1):
                x0 = x[mi*nfil + k]
                y0 = y[mi*nfil + k]
                x1 = x[mi*nfil + k + 1]
                y1 = y[mi*nfil + k + 1]
                dx = x0 - x1
                dy = y0 - y1
                xc = 0.5*(x0 + x1)
                yc = 0.5*(y0 + y1)
                ti = ti + (xc-cxi)*dy - (yc-cyi)*dx
        torque[i] = ti
        # compute pseudo velocities based on gamma
        gamma = 2.
        vxi = fix/gamma/natoms
        vyi = fiy/gamma/natoms
        wi = 2*ti/gamma/rgisq
        vx[i] = vxi
        vy[i] = vyi
        omega[i] = wi
        # compute kinetic rotation and translation energies based on m ( <- possibly nonsense)
        m = 1.0
        eti = 0.5*natoms*m*(vxi**2 + vyi**2)
        eri = 0.5*m*wi**2*natoms*rgisq
        etrans[i] = eti
        erot[i] = eri
        # compute the total external force
        fext = 0
        for j in range(ncl):
            mi - clusters[i][j]
            for k in range(nfil - 1):
                x0 = x[mi*nfil + k]
                y0 = y[mi*nfil + k]
                x1 = x[mi*nfil + k + 1]
                y1 = y[mi*nfil + k + 1]
                dx = x0 - x1
                dy = y0 - y1
                fext += np.sqrt(dx**2 + dy**2)
        ftotal[i] = fext
        # compute the total translational force
        ftr = np.sqrt(fix**2 + fiy**2)
        ftrans[i] = ftr
        # compute the total rotational force with respect to com
        fr = 0
        dfx_com = fix/ftr  # direction of the com force
        dfy_com = fiy/ftr
        for j in range(ncl):
            mi = clusters[i][j]
            for k in range(nfil - 1):
                # force along the bond
                x0 = x[mi*nfil + k]
                y0 = y[mi*nfil + k]
                x1 = x[mi*nfil + k + 1]
                y1 = y[mi*nfil + k + 1]
                dx = x0 - x1
                dy = y0 - y1
                # project bond force on translational force
                pr = dx*dfx_com + dy*dfy_com
                dx = dx*(1-pr)
                dy = dy*(1-pr)
                xc = 0.5*(x0 + x1)
                yc = 0.5*(y0 + y1)
                fr += ((xc-cxi)*dy - (yc-cyi)*dx)/np.sqrt((xc-cxi)**2 + (yc-cyi)**2)
        frot[i] = math.fabs(fr)
        # compute the total rotational force with reference point such
        #   that fr is maximized
        """
        res = optimize.minimize(make_rot_max(x,y,ncl,clusters,i,dfx_com,dfy_com), [cxi,cyi])
        xri = res.x[0]
        yri = res.x[1]
        fr_max = 1./res.fun
        xrot[i] = xri
        yrot[i] = yri
        frot_max[i] = fr_max
        """
        # compute the straightness and swirlicity
        straightness[i] = ftrans[i]/ftotal[i]
        swirlicity[i] = frot_max[i]/ftotal[i]
    return (cl_mol, forcex, forcey, torque,
            comx, comy, rgyrsq, vx, vy, omega,
            ftotal, ftrans, frot, frot_max, xrot, yrot,
            etrans, erot, straightness, swirlicity)

##################################################################

def gen_cluster_force_plot(x,y,clusters,
                           cl_mol, forcex, forcey, torque,
                           comx, comy, rgyrsq, vx, vy, omega,
                           ftotal, ftrans, frot, frot_max, xrot, yrot,
                           etrans, erot, straightness, swirlicity,tstep):
    """ create a plot for a few clusters, put torque and forces
        to legend """
    # loop over all clusters
    for i in range(len(clusters)):
        # check whether to create a plot for the current cluster
        #if len(clusters[i]) < 3:
        #    continue
        # generate a plot with the relevant coordinates
        xs = []
        ys = []
        for j in range(len(clusters[i])):
            mol = clusters[i][j]
            for k in range(nfil):
                xs.append(x[mol*nfil + k])
                ys.append(y[mol*nfil + k])
        xs = np.array(xs)
        ys = np.array(ys)

        # text for the label
        fx = forcex[i]
        fy = forcey[i]
        t = torque[i]
        rgsq = rgyrsq[i]
        wi = omega[i]
        eti = etrans[i]
        eri = erot[i]
        st = straightness[i]
        si_com = frot[i]/ftotal[i]
        label = ('straightness = ' + str(st) + '\n' +
                 'swirlicity = ' + str(si_com)
                )
        # generate the figure
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(xs, ys, ls = '', marker = 'o', markersize = 2, markeredgewidth = 0, color = 'r',
                label = label)
        ax.plot(comx[i],comy[i], ls = '', marker = 's', markersize = 5, color = 'k', label = 'com')
        #ax.plot(xrot[i],yrot[i], ls = '', marker = '^', markersize = 5, color = 'b', label = 'crot')
        ax.legend()
        ax.axis('equal')
        os.system('mkdir test > /dev/NULL')
        os.system('mkdir test/' + str(tstep) + ' > dev/NULL')
        plt.savefig('test/' + str(tstep) + '/cl_' + str(i) + '.png')
        plt.show()
        plt.close()
    return


##################################################################

def cluster_analysis(charfile, headerfile, nfil, dcrit, lcrit, pcrit):
    """ perform cluster analysis"""
    ### open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    ### generate information about molecule IDs
    mol, nmol = gen_mol_info(natoms, nfil)
    ### allocate arrays required for processing the data


    ### allocat arrays to store all data

    
    

    ### skip some snapshots for testing purpose
    nskip = nsteps - 1000
    read_char.skip_snapshots(hfile, ifile, nskip)
    
    ### loop over all steps
    for i in range(nsteps-nskip):
        ### print stats
        print 'Progress:',i,'/',nsteps
        ### read in the data
        xs,ys,lx,ly,tstep,natoms = read_char.read_snapshot(hfile, ifile)
        x = xs*lx
        y = ys*ly
        ### compute orientation of each bead, incorporate pbc
        phi = compute_orientation(x,y,lx,ly,nfil)
        ### find clusters
        clusters_atm = find_clusters(x,y,phi,mol,nmol,nfil,lx,ly,dcrit,lcrit,pcrit)



        # generate a plot of the clusters
        gen_cluster_plot(x,y,clusters_atm,str(tstep))


        ### Transform data representation of the clusters
        clusters = transform_cluster_data(clusters_atm, nfil)

        ### Compute coordinates of clusters with corrections for pbcs
        xcluster, ycluster, isolated = correct_cluster_pbc(x,y,lx,ly,clusters,nfil)

        ### examine properties of the cluster
        (cl_mol, forcex, forcey, torque,
         comx, comy, rgyrsq, vx, vy, omega,
         ftotal, ftrans, frot, frot_max, xrot, yrot,
         etrans, erot, straightness, swirlicity) = characterize_cluster(xcluster,ycluster,phi,clusters,lx,ly,nfil)

        ### generate plot which shows the 10 largest clusters and their external forces
        gen_cluster_force_plot(xcluster,ycluster,clusters,
                               cl_mol, forcex, forcey, torque,
                               comx, comy, rgyrsq, vx, vy, omega,
                               ftotal, ftrans, frot, frot_max, xrot, yrot,
                               etrans, erot, straightness, swirlicity, tstep)

    ifile.close()
    # still need to decide what exactly will be returned
    return 0

##################################################################

def main():
    """ main function, called when script is started"""
    ### read parameters from input file
    charfile, headerfile, nfil, dcrit, lcrit, pcrit = read_settings()
    # analyse the clusters of the simulation
    results = cluster_analysis(charfile, headerfile, nfil, dcrit, lcrit, pcrit)
    return

##################################################################

if __name__ == '__main__':
    main()
