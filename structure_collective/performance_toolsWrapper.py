
import ctypes

lib = ctypes.cdll.LoadLibrary('/usr/users/iff_th2/duman/RolfData/Scripts/structure_collective/libperformance_tools.so')

class Performance_tools(object):
    
    def __init__(self):
        self.obj = lib.Performance_tools_new()
        
    ########################################################################

    def compute_static_structure(self, S, x, y, nsteps, natoms, box_size):

        # transfrom arrays to ctypes
        Sc = S.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        xc = S.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        yc = S.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        natomsc = ctypes.c_int(natoms)
        box_sizec = ctypes.c_double(box_size)

        lib.compute_static_structure(self.obj, Sc, xc, yc, nstepsc, natomsc, box_sizec)
        
        return
        
