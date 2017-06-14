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
import performance_toolsWrapper
#import bottleneck as bot

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
    # output-folder name
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
    # critical distance for mol/mol neighbor search
    line = ifile.readline()
    line = line.split()
    dcrit = float(line[-1])
    # critical overlap for mol/mol neighbor search
    line = ifile.readline()
    line = line.split()
    lcrit = float(line[-1])
    # critical angle for mol/mol neighbor search
    line = ifile.readline()
    line = line.split()
    pcrit = float(line[-1])
    pcrit *= np.pi/180.0
    ccrit = math.cos(pcrit)
    # critical width for the histogram
    line = ifile.readline()
    line = line.split()
    wcrit = float(line[-1])
    # close file
    ifile.close()
    # return input values
    return charfile, headerfile, ofname, nfil, nskip,dcrit, lcrit, ccrit, wcrit

############################################################################

def gen_mol_info(natoms, nfil):
    """ generate information about molecular ID"""
    mol = np.zeros((natoms), dtype = np.int32)
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

def nearbyint(x):
    """ Round to the nearby integer"""
    if x >= 0:
        return math.floor(x+0.5)
    else:
        return math.floor(x-0.5)

############################################################################
    
def minimum_image_dist(x1,x2,lx):
    """ compute vector to nearest neighbor"""
    dx = x2 - x1 
    return dx-nearbyint(dx/lx)*lx

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
            tx[k] = dx / np.sqrt(dx**2 + dy**2)
            ty[k] = dy / np.sqrt(dx**2 + dy**2)
            # increment k
            k = k + 1
    return phi, tx, ty

############################################################################

def gen_linked_list(x,y,lx,ly,dcrit):
    """ generate a linked list"""
    # determine the number of cells in each direction
    nsegx = int(lx/dcrit)
    nsegy = int(ly/dcrit)

    # allocate head and llist
    ncells = nsegx*nsegy
    natoms = len(x)
    head = np.zeros((ncells), dtype = np.int32)
    llist = np.zeros((natoms), dtype = np.int32)

    # fill list and head
    for i in range(natoms):
        segx = int(x[i]/lx*nsegx)
        segy = int(y[i]/ly*nsegy)
        cell = segx*nsegy + segy
        llist[i] = head[cell]
        head[cell] = i

    return nsegx,nsegy,head,llist

############################################################################
    
def sample_new_cluster(cl_mol,maxs):
    """ decide whether to do a new cluster index sampling or not"""
    sample = True
    for siz in cl_mol:
        for maxn in maxs:
            if siz == maxn:
                sample = False
                return sample
    return sample
   
############################################################################
    
def compare_clusters(cluster_idx, old_cluster_idx, nfil, lx, ly, ofile, unique_id):
    """ compare the clusters in different time frames"""
        
#        ## compute the minimum image distance between center of masses with each cluster from the previous point in time
#        dx = np.zeros((len(old_cluster_idx)))
#        dy = np.zeros((len(old_cluster_idx)))            
#        for oj in range(len(old_cluster_idx)):
#            dx[oj] = minimum_image_dist(cluster_idx[j][2],old_cluster_idx[oj][2],lx)
#            dy[oj] = minimum_image_dist(cluster_idx[j][3],old_cluster_idx[oj][3],ly)
#        d2 = dx**2 + dy**2
#        
#        ## get the 5 closest clusters in terms of com distance
#        mins = bot.partsort(d2, 5)[:5]         # 5 closest clusters
#        
#        ## get the indices of the 5 closest clusters
#        closest_clusters = np.zeros(5, dtype=int)
#        for i, dist2 in enumerate(d2):
#            if dist2 == mins[0]:
#                closest_clusters[0] = i
#            elif dist2 == mins[1]:
#                closest_clusters[1] = i
#            elif dist2 == mins[2]:
#                closest_clusters[2] = i
#            elif dist2 == mins[3]:
#                closest_clusters[3] = i
#            elif dist2 == mins[4]:
#                closest_clusters[4] = i

#        print 'Index of the closest clusters are\t', closest_clusters, ' with the comparison cluster index being ', j
#        print 'Minimum distances are\t', d2[closest_clusters[0]], d2[closest_clusters[1]], d2[closest_clusters[2]], d2[closest_clusters[3]], d2[closest_clusters[4]]
#        print 'Position of the old cluster at the minimum distances are: '
#        print 'Old cluster 1: ', old_cluster_idx[closest_clusters[0]][2], old_cluster_idx[closest_clusters[0]][3], old_cluster_idx[closest_clusters[0]][9]
#        print 'Old cluster 2: ', old_cluster_idx[closest_clusters[1]][2], old_cluster_idx[closest_clusters[1]][3], old_cluster_idx[closest_clusters[1]][9]
#        print 'Old cluster 3: ', old_cluster_idx[closest_clusters[2]][2], old_cluster_idx[closest_clusters[2]][3], old_cluster_idx[closest_clusters[2]][9]
#        print 'Old cluster 4: ', old_cluster_idx[closest_clusters[3]][2], old_cluster_idx[closest_clusters[3]][3], old_cluster_idx[closest_clusters[3]][9]
#        print 'Old cluster 5: ', old_cluster_idx[closest_clusters[4]][2], old_cluster_idx[closest_clusters[4]][3], old_cluster_idx[closest_clusters[4]][9]
#        print 'Position of the cluster being compared is:\t', cluster_idx[j][2], cluster_idx[j][3]
        
#        ## among the closest 5 clusters, choose the cluster with the largest number of intersecting/overlapping filaments
#        len_intersect_clusters = np.zeros(5)
#        for i, idx in enumerate(closest_clusters):
#            len_intersect_clusters[i] = len(np.intersect1d(cluster_idx[j][9], old_cluster_idx[int(idx)][9]))
#        maxi_intersect = max(len_intersect_clusters)
#        maxi_intersect_idx = int(np.argmax(len_intersect_clusters))
#        guess_idx = closest_clusters[maxi_intersect_idx]

#        ## choose the cluster with the number of overlapping filaments larger than the size
#        if maxi_intersect < old_cluster_idx[guess_idx][1]/2:
#            len_intersect_clusters = np.zeros(len(old_cluster_idx))
#            for oj in range(len(old_cluster_idx)):
#                len_intersect_clusters[oj] = len(np.intersect1d(cluster_idx[j][9], old_cluster_idx[oj][9]))
#            maxi_intersect = max(len_intersect_clusters)
#            maxi_intersect_idx = int(np.argmax(len_intersect_clusters))  
#            if maxi_intersect > 5:
#                guess_idx = maxi_intersect_idx 
#            else:
#                jcnt += 1
#                guess_idx = -1
                
    visited = np.ones(len(old_cluster_idx), dtype=int)*(-1)

    ## for each cluster with identities stored in cluster_idx
    for j in range(len(cluster_idx)):                
                
        ## choose the cluster with the maximum number of overlapping filaments 
        len_intersect_clusters = np.zeros(len(old_cluster_idx))
        for oj in range(len(old_cluster_idx)):
            len_intersect_clusters[oj] = len(np.intersect1d(cluster_idx[j][9], old_cluster_idx[oj][9]))
        if len(len_intersect_clusters) != 0:
            maxi_intersect = max(len_intersect_clusters)
            maxi_intersect_idx = int(np.argmax(len_intersect_clusters))  
        else:
            maxi_intersect = 0

        ## if the maximum intersection size is larger than the given threshold size        
        if maxi_intersect > 5:
            
            ## maximum intersecting id from the previous time frame is the guess for now
            guess_idx = maxi_intersect_idx 
            
            ## check if the found cluster id is already occupied or not
            ## if it is not occupied, occupy it
            if visited[guess_idx] == -1:
                visited[guess_idx] = j
                
            ## if it is occupied, compare with the occupying cluster id
            else:
                len_intersect_new = len(np.intersect1d(cluster_idx[j][9], old_cluster_idx[guess_idx][9]))
                len_intersect_old = len(np.intersect1d(cluster_idx[visited[guess_idx]][9], old_cluster_idx[guess_idx][9]))
                
                ## in case of clash of cluster ids, choose the id with larger overlapping filaments and destroy the other one
                if len_intersect_new > len_intersect_old:
                    ## assign a new unique id for the clashing cluster
                    unique_id += 1
                    cluster_idx[visited[guess_idx]][0] = unique_id
                    ## change the visiting cluster to the new cluster
                    visited[guess_idx] = j
                    
                ## if the older one is larger, make a unique cluster id for the next one and leave the older one intact
                else:
                    unique_id += 1
                    guess_idx = -1
                    
        ## if the maximum intersection size is smaller than the threshold, then find a unique cluster id
        else:
            unique_id += 1
            guess_idx = -1     
        
        ## update the cluster if unique id flag is not on
        if guess_idx != -1:
#            print 'Cluster being tested is: ', j, ' with positions: ', cluster_idx[j][2], cluster_idx[j][3]
#            print 'and identities: ', cluster_idx[j][9]
#            print 'Cluster identity is changed to: ', guess_idx, ' with positions in the previous step: ', old_cluster_idx[guess_idx][2], old_cluster_idx[guess_idx][3]
#            print 'and identities: ', old_cluster_idx[guess_idx][9]
#            print 'Maximum number of overlapping filaments is: ', maxi_intersect
            cluster_idx[j][0] = old_cluster_idx[guess_idx][0]
        else:
            cluster_idx[j][0] = unique_id
    
    ## write the data after the all the update is conducted!    
    for j in range(len(cluster_idx)):
        ofile.write(str(cluster_idx[j][0]) + '\t' + str(cluster_idx[j][1]) + '\t'
        + str(cluster_idx[j][2]) + '\t' + str(cluster_idx[j][3]) + '\t'
        + str(cluster_idx[j][4]) + '\t' + str(cluster_idx[j][5]) + '\t'
        + str(cluster_idx[j][6]) + '\t' + str(cluster_idx[j][7]) + '\t'
        + str(cluster_idx[j][8]) + '\t')        
        for kj in range(len(cluster_idx[j][9])):
            ofile.write(str(cluster_idx[j][9][kj]) + '\t')
        ofile.write('\n')

    return unique_id
    
############################################################################

def find_clusters(x,y,tx,ty,mol,nmol,nfil,lx,ly,dcrit,lcrit,ccrit):
    """ search for clusters; two molecules are defined as being
        part of the same cluster if lcrit of their body length is
        within a distance of dcrit and when the difference in orientation
        is less than phi"""
    natoms = len(x)
    # allocate required arrays
    cl = np.zeros((nmol), dtype = np.int32)
    cl -= 1
    clusters = np.zeros((natoms), dtype = np.int32)
    neighs_mol = np.zeros((nmol,nmol), dtype = np.int32)
    # generate a linked list
    nsegx, nsegy, head, llist = gen_linked_list(x,y,lx,ly,dcrit)
    # find neighborhood information
    neighs_molf = neighs_mol.ravel()
    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.fill_neigh_matrix(neighs_molf,llist,head,nsegx,nsegy,x,y,tx,ty,mol,nmol,lx,ly,dcrit,ccrit)
    # recursive search for clusters in neighbor matrix
    performance_tools.cluster_search(neighs_molf,cl,lcrit,nfil,nmol)
    # fill cluster results to per atom array
    k = 0
    for i in range(nmol):
        for j in range(nfil):
            clusters[k] = cl[i]
            k = k + 1
    ### return results
    return clusters

##################################################################

def gen_cluster_plot(x,y,clusters,tstep,phi,lx,ly,ofname):
    """ generate a plot of the atoms color coded with the cluster they belong"""
    
    ofile = open(ofname + '/beads_' + tstep + '.txt', 'w')
    for i in range(len(x)):
        ofile.write(str(x[i]) + '\t' + str(y[i]) + '\t' + str(phi[i]) + '\t' + str(clusters[i]) + '\n')
    ofile.close()
        
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
    os.system('~/Applications/Scripts/clusters_collective/recurse_cluster')
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
    h2 = np.zeros((nx,ny), dtype = np.int32)
    # set starting value for h2
    h2[sx0,sy0] = 1
    # fill h2 iteratively,
    performance_tools = performance_toolsWrapper.Performance_tools()
    performance_tools.recurse_cluster(h,h2,nx,ny,sx0,sy0)
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

def cluster_hist_1d(xcl,ycl,lx,ly,nmol,dcrit,nfil):
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

def cluster_hist_2d(xcl,ycl,lx,ly,nmol,dcrit,nfil,ofname):
    """ try to cluster based on 2d histogramms"""
    nxbins = int(3*lx/(2*dcrit))
    bx = 3*lx/nxbins
    nybins = int(3*ly/(2*dcrit))
    by = 3*ly/nybins
    hxy, xedges, yedges = np.histogram2d(xcl, ycl, bins = [nxbins,nybins], range = [[-lx,2*lx], [-ly,2*ly]])
    hxy *= 2
    hxy = hxy.astype(np.int32)
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
    if 8*np.sum(hxy2) < np.sum(hxy):
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
            dy = nearest_neighbor(y0,y1,ly)
            xcl[l+1] = x0 - dx
            ycl[l+1] = y0 - dy
            l = l + 1
        l = l + 1
    # loop over all molecules in the cluster
    #   adjust com position
    for j in range(nmol):
        comx = np.average(xcl[j*nfil:(j+1)*nfil - 1])
        comy = np.average(ycl[j*nfil:(j+1)*nfil - 1])
        xcl[j*nfil:(j+1)*nfil - 1] += -math.floor(comx/lx)*lx
        ycl[j*nfil:(j+1)*nfil - 1] += -math.floor(comy/ly)*ly
    return

##################################################################

def adjust_com_cluster(xcl,ycl,lx,ly,nmol,nfil):
    """ move cluster such that com is in periodic box"""
    comx = np.average(xcl[0:nmol*nfil])
    comy = np.average(ycl[0:nmol*nfil])
    xcl[0:nmol*nfil] += -math.floor(comx/lx)*lx
    ycl[0:nmol*nfil] += -math.floor(comy/ly)*ly

    return

##################################################################

def correct_cluster_pbc(x,y,lx,ly,clusters,nfil,dcrit,ofname):
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
        flag_1d = cluster_hist_1d(xcl,ycl,lx,ly,nmol,dcrit,nfil)
        if flag_1d == 0:
            isolated[i] = 0
            correct_pbc_single(xcl,ycl,nmol,lx,ly,nfil)
#        # try approach with 2D histgrams if 1d failed
#        if flag_1d == 0:
#            flag_2d = cluster_hist_2d(xcl,ycl,lx,ly,nmol,dcrit,nfil,ofname)
#        # move molecules if both approaches failed
#        if flag_1d == 0 and flag_2d == 0:
#            isolated[i] = 0
#            correct_pbc_single(xcl,ycl,nmol,lx,ly,nfil)
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

def make_rot_max(x,y,ncl,clusters,i,dfx_com,dfy_com,nfil):
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

def cluster_size(clusters):
    """ compute the size of all clusters"""
    n = len(clusters)
    cl_mol = np.zeros((n), dtype = int)
    for i in range(n):
        ncl = len(clusters[i])
        cl_mol[i] = ncl
    return cl_mol

##################################################################

def characterize_cluster(x,y,phi,clusters,nfil):
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
    #frot_max = np.zeros(n)
    xrot = np.zeros(n)
    yrot = np.zeros(n)
    etrans = np.zeros(n)
    erot = np.zeros(n)
    straightness = np.zeros(n)
    swirlicity = np.zeros(n)
    asphericity = np.zeros(n)
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
        rgidx = 0
        rgidy = 0
        rgidxdy = 0
        for j in range(ncl):
            mi = clusters[i][j]
            for k in range(nfil): 
                dx = x[mi*nfil + k] - cxi
                dy = y[mi*nfil + k] - cyi
                rgisq += dx**2 + dy**2
                rgidx += dx**2
                rgidy += dy**2
                rgidxdy += dx*dy
        rgisq /= natoms
        rgidx /= natoms
        rgidy /= natoms
        rgidxdy /= natoms
        det = rgidx*rgidy - rgidxdy**2
        trace = rgidx + rgidy
        #avg_rgi = 0.5*(rgidx+rgidy)
        rgyrsq[i] = rgisq
        asphericity[i] = 1 - 4*det/trace**2
        #asphericity[i] = 0.5*( (rgidx-avg_rgi)**2 + (rgidy-avg_rgi)**2 )/avg_rgi**2
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
        fext = 0
        for j in range(ncl):
            mi = clusters[i][j]
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
        swirlicity[i] = frot[i]/ftotal[i]
        #swirlicity[i] = frot_max[i]/ftotal[i]
    return (cl_mol, forcex, forcey, torque,
            comx, comy, rgyrsq, vx, vy, omega,
            ftotal, ftrans, frot, xrot, yrot,
            etrans, erot, straightness, swirlicity, asphericity)



##################################################################

def compute_rotation(vx,vy,rho, bx,by,lx,ly):
    """ compute the rotation of a vector field"""
    # allocate output array
    rotation = np.zeros((bx,by))
    duy_dx = np.zeros((bx,by))
    dux_dy = np.zeros((bx,by))
    # bin width
    wx = lx/bx
    wy = ly/by
    # number of included elements
    included = 0
    # compute rotation in each bin
    for i in range(bx):
        for j in range(by):
            # compute velocity gradients using finite differences
            vyxp = vy[(i+1)%bx,j]
            vyxm = vy[i-1,j]
            vxyp = vx[i,(j+1)%by]
            vxym = vx[i,j-1]
            if vyxp*vyxm*vxyp*vxym == 0:
                continue
            included = included + 1
            duy_dx[i,j] = (vyxp - vyxm)/(2*wx)
            dux_dy[i,j] = (vxyp - vxym)/(2*wy)
            rotation[i,j] = duy_dx[i,j] - dux_dy[i,j]
    return duy_dx, dux_dy, rotation, included

##################################################################

def characterize_particular_rotation(x, y, phi, tx, ty, ind, clusters_atm, nfil, width, tstep):
    """ characterize the rotation of a cluster based on the integral over the surface"""
    ### compute the local density and forces based on a histogram with width dcrit    
    xcl = x[clusters_atm == ind]
    ycl = y[clusters_atm == ind]
    txcl = tx[clusters_atm == ind]
    tycl = ty[clusters_atm == ind]
    natoms = len(xcl)
    xmin = min(xcl) - 2*width
    xmax = max(xcl) + 2*width
    ymin = min(ycl) - 2*width
    ymax = max(ycl) + 2*width
    lx = xmax - xmin
    ly = ymax - ymin
    nx = int(lx/width) + 1
    ny = int(ly/width) + 1
    xmax = nx*width + xmin
    ymax = ny*width + ymin
    lx = xmax - xmin
    ly = ymax - ymin
    ### generate a histogram of the density and the forces
    rho = np.zeros((nx,ny))
    fx_b = np.zeros((nx,ny))
    fy_b = np.zeros((nx,ny))
    for j in range(natoms):
        # get coordinates
        xi = xcl[j]
        yi = ycl[j]
        # get current bin
        segx = int((xi-xmin)/lx*nx)
        segy = int((yi-ymin)/ly*ny)
        # add data to bin
        rho[segx,segy] += 1
        fx_b[segx,segy] += txcl[j]
        fy_b[segx,segy] += tycl[j]
    ### divide by density to get polarity
    for ix in range(nx):
        for iy in range(ny):
            if rho[ix,iy] > 0:
                fx_b[ix,iy] /= rho[ix,iy]
                fy_b[ix,iy] /= rho[ix,iy]
    ### compute the rotation
    duy_dx, dux_dy, forticity, included = compute_rotation(fx_b, fy_b, rho, nx, ny, lx, ly)
    enstrophy_value = 0.
    if included > 0:
        #swirlicity = np.sum(forticity) / included
        enstrophy_value = 0.5*np.sum(forticity**2) / included

    return enstrophy_value
    
##################################################################


def cluster_analysis(charfile, headerfile, ofname, nfil, nskip, dcrit, lcrit, ccrit, wcrit):
    """ perform cluster analysis"""
    ### open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    ### get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    nsteps = nsteps - nskip
    if nsteps <= 0:
        print 'Error, nsteps < 0'
        exit()
        
    ### generate information about molecule IDs
    mol, nmol = gen_mol_info(natoms, nfil)

    ### allocate arrays to store the cluster statistics
    ncl = np.zeros((nsteps), dtype = np.int32)
    #clsize = [] # unfortunately a list    
    t = np.zeros((nsteps), dtype = np.int32)
    clmax = np.zeros((nsteps), dtype = np.int32)

    ### skip some initial snapshots
    read_char.skip_snapshots(hfile, ifile, nskip)
    
    first_time = True
    old_cluster_idx = []
    
    ### loop over all steps
    unique_id = 0    
    for i in range(nsteps):
        ### print stats
         
        ### read in the data
        xs,ys,lx,ly,tstep,natoms = read_char.read_snapshot(hfile, ifile)
        x = xs*lx
        y = ys*ly
        
#        ### check if the file is already examined or not -- this is wrong!
#        if os.path.exists(ofname + '/beads_' + str(tstep) + '.txt'):
#            print 'This step is passed'
#            continue
        
        print 'Progress:',i,'/',nsteps

        ### compute orientation of each bead, incorporate pbc
        phi,tx,ty = compute_orientation(x,y,lx,ly,nfil)
        
        ### find clusters
        clusters_atm = find_clusters(x,y,tx,ty,mol,nmol,nfil,lx,ly,dcrit,lcrit,ccrit)

        ### write data of cluster id of atoms to generate a plot of the clusters with post processsing
        ofile = open(ofname + '/beads_' + str(tstep) + '.txt', 'w')
        for ii in range(len(x)):
            ofile.write(str(x[ii]) + '\t' + str(y[ii]) + '\t' + str(phi[ii]) + '\t' + str(clusters_atm[ii]) + '\n')
        ofile.close()

        ### Transform data representation of the clusters
        clusters = transform_cluster_data(clusters_atm, nfil)

        ### compute cluster size
        cl_mol = cluster_size(clusters)
        
        ### store cluster sizes in form of number of filaments
        t[i] = tstep
        ofile = open(ofname + '/cluster_sizes_' + str(tstep) + '.txt', 'w')
        ncli = len(cl_mol)
        ncl[i] = ncli
        for j in range(ncli):
            #clsize.append(cl_mol[j])
            ofile.write(str(cl_mol[j]) + '\n')
        clmax[i] = max(cl_mol)
        ofile.close()        
        #print tstep, clmax[i]
        
        ### Compute coordinates of clusters with corrections for pbcs      
        xcluster, ycluster, isolated = correct_cluster_pbc(x,y,lx,ly,clusters,nfil,dcrit,ofname)
         
        ### examine properties of the cluster
        (cl_mol, forcex, forcey, torque,
         comx, comy, rgyrsq, vx, vy, omega,
         ftotal, ftrans, frot, xrot, yrot,
         etrans, erot, straightness, swirlicity, asphericity) = characterize_cluster(xcluster,ycluster,phi,clusters,nfil)        
         
        ## parametrize the clusters with more than 10 filaments
        ofile = open(ofname + '/cluster_evolution_' + str(tstep) + '.txt', 'w') 
        cluster_idx = []
        for j in range(ncl[i]):
            if cl_mol[j] >= 10:
                unique_id += 1
                enstrophy = characterize_particular_rotation(xcluster, ycluster, phi, tx, ty, j, clusters_atm, nfil, wcrit, tstep)
                cluster_idx.append([unique_id, cl_mol[j], comx[j], comy[j], rgyrsq[j], asphericity[j], straightness[j], swirlicity[j], enstrophy, clusters[j]])

        ## form the old cluster ids when running for the first time
        if first_time:
            first_time = False
            old_cluster_idx = cluster_idx
            for j in range(len(cluster_idx)):
                ofile.write(str(cluster_idx[j][0]) + '\t' + str(cluster_idx[j][1]) + '\t'
                + str(cluster_idx[j][2]) + '\t' + str(cluster_idx[j][3]) + '\t' 
                + str(cluster_idx[j][4]) + '\t' + str(cluster_idx[j][5]) + '\t'
                + str(cluster_idx[j][6]) + '\t' + str(cluster_idx[j][7]) + '\t'
                + str(cluster_idx[j][8]) + '\t')
                for kj in range(len(cluster_idx[j][9])):
                    ofile.write(str(cluster_idx[j][9][kj]) + '\t')
                ofile.write('\n')
            ofile.close()
            continue
        
        ## compare clusters with the previous time step and change ids accordingly
        unique_id = compare_clusters(cluster_idx, old_cluster_idx, nfil, lx, ly, ofile, unique_id)
        ofile.close()    
        old_cluster_idx = cluster_idx
           
#    # transform lists to arrays
#    clsize = np.array(clsize)
    
    # close input files
    ifile.close()
    hfile.close()
    
    # return list with cluster sizes, maybe more in the future
    return t, ncl, clmax            
                
                
#        ## Characterize the top 5 clusters which are still in the largest top 10 clusters in the next few steps
#        maxs = -bot.partsort(-cl_mol, 10)[:10]         # 10 largest clusters
#        maxs = sorted(maxs, reverse=True)   
#        
#        if its_first_time:
#            its_first_time = False
#            maxs = maxs[:5] 
#            maxs_idx = []
#            maxs_idx = [idx for idx, siz in enumerate(cl_mol) for maxn in maxs if siz == maxn]
#            max_clusters_idx = np.asarray(maxs_idx)
#            max_clusters = np.array( [clusters[maxs_idx[0]], clusters[maxs_idx[1]], 
#                                      clusters[maxs_idx[2]], clusters[maxs_idx[3]], clusters[maxs[4]]] )
#            new_cluster = False
#        else:    
#            new_cluster = sample_new_cluster(cl_mol,maxs)
#        
#        print max_clusters, tstep
#        print max_clusters_idx, tstep
#        if new_cluster:
#            
#            maxs = maxs[:5] 
#            maxs_idx = []
#            maxs_idx = [idx for idx, siz in enumerate(cl_mol) for maxn in maxs if siz == maxn]
#            max_clusters_idx = np.asarray(maxs_idx)
#            max_clusters = np.array( [len(clusters[maxs_idx[0]]), len(clusters[maxs_idx[1]]), 
#                                      len(clusters[maxs_idx[2]]), len(clusters[maxs_idx[3]]), len(clusters[maxs_idx[4]])] )
#
#        ### Keep the identity of the maximum cluster in the next steps if it is in the largest 10 clusters in the next step     
#        
#        ### Compute coordinates of clusters with corrections for pbcs      
#        xcluster, ycluster, isolated = correct_cluster_pbc(x,y,lx,ly,max_clusters,nfil,dcrit,ofname)
#
#        ### examine properties of the cluster
#        (max_cl_mol, forcex, forcey, torque,
#         comx, comy, rgyrsq, vx, vy, omega,
#         ftotal, ftrans, frot, xrot, yrot,
#         etrans, erot, straightness, swirlicity) = characterize_cluster(xcluster,ycluster,phi,max_clusters,nfil)        
#
#        ### characterize rotation
#        #characterize_rotation(xcluster, ycluster, phi, tx, ty, clusters, clusters_atm, nfil, wcrit, ofname,tstep)
#                           
#        ### generate plot which shows the 10 largest clusters and their external forces and store data about them
#        gen_cluster_force_plot(xcluster,ycluster,tx,ty,max_clusters,max_clusters_idx,clusters_atm,
#                               max_cl_mol, forcex, forcey, torque,
#                               comx, comy, rgyrsq, vx, vy, omega,
#                               ftotal, ftrans, frot, xrot, yrot,
#                               etrans, erot, straightness, swirlicity, dcrit,tstep,nfil,ofname)

##################################################################

def main():
    """ main function, called when script is started"""
    ### read parameters from input file
    charfile, headerfile, ofname, nfil, nskip, dcrit, lcrit, ccrit, wcrit = read_settings()
    ### analyse the clusters of the simulation
    if os.path.exists(ofname) != True:
        os.system('mkdir ' + ofname)
    t, ncl, clmax = cluster_analysis(charfile, headerfile, ofname, nfil, nskip, dcrit, lcrit, ccrit, wcrit)
    ### write results to files
    # plot the time evolution of the largest cluster
    ax = plt.subplot(111)
    ax.plot(t,clmax)
    ax.set_xlabel('timestep')
    ax.set_ylabel('maximum cluster size')
    plt.savefig(ofname + '/clmax.png')
    plt.close()
#    # plot a histogram of the cluster size
#    ax = plt.subplot(111)
#    ax.hist(clsize, bins = np.sqrt(len(clsize)-1))
#    ax.set_xlabel('cluster size [mol]')
#    ax.set_ylabel('hist cluster size')
#    plt.savefig(ofname + '/histogram.png')
#    plt.close()
    return

##################################################################

if __name__ == '__main__':
    main()
