import ctypes
lib = ctypes.cdll.LoadLibrary('/usr/users/iff_th2/duman/RolfData/Scripts/MSR_collective/libperformance_tools.so')

class Performance_tools(object):
    def __init__(self):
        self.obj = lib.Performance_tools_new()

    ############################################################################

    def compute_MSR(self, t_msr, time, phi, msr_tmp, counter, logval, nsteps, nmsr, limit, nevery):

        # transfrom arrays to ctypes
        t_msrc = t_msr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        phic = phi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        msr_tmpc = msr_tmp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        counterc = counter.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        logvalc = logval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        nmsrc = ctypes.c_int(nmsr)
        limitc = ctypes.c_int(limit)
        neveryc = ctypes.c_int(nevery)

        lib.compute_msr(self.obj, t_msrc, timec, phic, msr_tmpc, counterc, logvalc, nstepsc, nmsrc, limitc, neveryc)
        return

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
