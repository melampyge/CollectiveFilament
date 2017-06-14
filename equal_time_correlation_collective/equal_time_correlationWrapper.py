import ctypes
lib = ctypes.cdll.LoadLibrary('/usr/users/iff_th2/duman/RolfData/Scripts/equal_time_correlation_collective/libequal_time_correlation.so')

class Equal_time_correlation(object):
    def __init__(self):
        self.obj = lib.Equal_time_correlation_new()

    def compute(self,nsegx,nsegy,natoms, head, llist,
                mol, x, y, etcorr, lx, ly, rmax, nbin, tx, ty, counter):

        headc = head.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        llistc = llist.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        molc = mol.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        xc = x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        yc = y.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        etcorrc = etcorr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        txc = tx.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        tyc = ty.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        counterc = counter.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        nsegxc = ctypes.c_int(nsegx)
        nsegyc = ctypes.c_int(nsegy)
        natomsc = ctypes.c_int(natoms)
        lxc = ctypes.c_double(lx)
        lyc = ctypes.c_double(ly)
        rmaxc = ctypes.c_double(rmax)
        nbinc = ctypes.c_int(nbin)
        
        lib.compute(self.obj, nsegxc, nsegyc, natomsc, headc, llistc,
                    molc, xc, yc, etcorrc, lxc, lyc, rmaxc, nbinc, txc, tyc, counterc)
