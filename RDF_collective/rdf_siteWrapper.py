import ctypes
lib = ctypes.cdll.LoadLibrary('/usr/users/iff_th2/isele/Applications/Scripts/RDF_collective/librdf_site.so')

class Rdf_site(object):
    def __init__(self):
        self.obj = lib.Rdf_site_new()

    def compute(self,nsegx,nsegy,natoms, head, llist,
                mol, x, y, rdf, lx, ly, rmax, nbin):

        headc = head.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        llistc = llist.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        molc = mol.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        xc = x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        yc = y.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        rdfc = rdf.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        nsegxc = ctypes.c_int(nsegx)
        nsegyc = ctypes.c_int(nsegy)
        natomsc = ctypes.c_int(natoms)
        lxc = ctypes.c_double(lx)
        lyc = ctypes.c_double(ly)
        rmaxc = ctypes.c_double(rmax)
        nbinc = ctypes.c_int(nbin)
        
        lib.compute(self.obj, nsegxc, nsegyc, natomsc, headc, llistc,
                    molc, xc, yc, rdfc, lxc, lyc, rmaxc, nbinc)
