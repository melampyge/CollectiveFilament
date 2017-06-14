import ctypes
lib = ctypes.cdll.LoadLibrary('/usr/users/iff_th2/duman/RolfData/Scripts/clusters_collective/libperformance_tools.so')

class Performance_tools(object):
    def __init__(self):
        self.obj = lib.Performance_tools_new()

    ########################################################################
        
    def fill_neigh_matrix(self,neighs_molf,llist,head,nsegx,nsegy,x,y,tx,ty,mol,nmol,lx,ly,dcrit,ccrit):

        neighs_molfc = neighs_molf.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        headc = head.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        llistc = llist.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        xc = x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        yc = y.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        txc = tx.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        tyc = ty.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        molc = mol.ctypes.data_as(ctypes.POINTER(ctypes.c_int))


        nsegxc = ctypes.c_int(nsegx)
        nsegyc = ctypes.c_int(nsegy)
        nmolc = ctypes.c_int(nmol)
        lxc = ctypes.c_double(lx)
        lyc = ctypes.c_double(ly)
        dcritc = ctypes.c_double(dcrit)
        ccritc = ctypes.c_double(ccrit)
        
        lib.fill_neigh_matrix(self.obj, neighs_molfc, llistc, headc,
                              nsegxc, nsegyc, xc, yc, txc, tyc, molc, nmolc, lxc,
                              lyc, dcritc, ccritc)
        return

    
    ########################################################################

    def cluster_search(self, neighs, cl, lcrit, nfil, nmol):

        neighsc = neighs.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        clc = cl.ctypes.data_as(ctypes.POINTER(ctypes.c_int))

        lcritc = ctypes.c_double(lcrit)
        nfilc = ctypes.c_int(nfil)
        nmolc = ctypes.c_int(nmol)

        lib.cluster_search(self.obj, neighsc, clc, lcritc, nfilc, nmolc)
        return

    ########################################################################

    def recurse_cluster(self, h, h2, nx, ny, sx0, sy0):

        # flatten arrays
        hf = h.ravel()
        h2f = h2.ravel()
        # transform arrays to ctypes
        hc = hf.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        h2c = h2f.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform other entries to ctypes
        nxc = ctypes.c_int(nx)
        nyc = ctypes.c_int(ny)
        sx0c = ctypes.c_int(sx0)
        sy0c = ctypes.c_int(sy0)

        lib.recurse_cluster(self.obj, hc, h2c, nxc, nyc, sx0c, sy0c)
        return
        

        
