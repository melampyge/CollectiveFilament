import ctypes
lib = ctypes.cdll.LoadLibrary('/usr/users/iff_th2/isele/Applications/Scripts/rod_swimmer_compute_orientations/libperformance_tools.so')

class Performance_tools(object):
    def __init__(self):
        self.obj = lib.Performance_tools_new()


    ############################################################################

    def compute_angle(self, ex, ey, phi, dp, nsteps):

        # transfrom arrays to ctypes
        exc = ex.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        eyc = ey.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        phic = phi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        dpc = dp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        # transform integers to ctypes
        nstepsc = ctypes.c_int(nsteps)

        lib.compute_angle(self.obj, exc, eyc, phic, dpc, nsteps)
        return

