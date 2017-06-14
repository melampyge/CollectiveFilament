import ctypes
lib = ctypes.cdll.LoadLibrary('/usr/users/iff_th2/isele/Applications/Scripts/sphere_swimmer_compute_orientations/libperformance_tools.so')

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


    ############################################################################

    def compute_MSR_loglog(self, t_msr, time, phi, msr, counter, logval, nsteps, nmsr, limit):

        # transfrom arrays to ctypes
        t_msrc = t_msr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        phic = phi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        msrc = msr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        counterc = counter.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        logvalc = logval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        nmsrc = ctypes.c_int(nmsr)
        limitc = ctypes.c_int(limit)

        lib.compute_msr_loglog(self.obj, t_msrc, timec, phic, msrc, counterc, logvalc, nstepsc, nmsrc, limitc)
        return

    ############################################################################

    def compute_correlation(self, t_ee, time, tx, ty, c_ee, counter, linval, nsteps, limit):

        # transfrom arrays to ctypes
        t_eec = t_ee.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        txc = tx.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        tyc = ty.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_eec = c_ee.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        counterc = counter.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        linvalc = linval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        limitc = ctypes.c_int(limit)

        lib.compute_correlation(self.obj, t_eec, timec, txc, tyc, c_eec, counterc, linvalc, nstepsc, limitc)
        return

    ############################################################################

    def compute_correlation_with_errorbars(self, t_ee, time, tx, ty, c_ee, std, blockdata, block_std, block_uncert, linval, nsteps, limit):

        # transfrom arrays to ctypes
        t_eec = t_ee.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        txc = tx.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        tyc = ty.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        c_eec = c_ee.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        stdc = std.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        blockdatac = blockdata.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        block_stdc = block_std.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        block_uncertc = block_uncert.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        linvalc = linval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        limitc = ctypes.c_int(limit)

        lib.compute_correlation_with_errorbars(self.obj, t_eec, timec, txc, tyc, c_eec, stdc, blockdatac, block_stdc, block_uncertc, linvalc, nstepsc, limitc)
        return

