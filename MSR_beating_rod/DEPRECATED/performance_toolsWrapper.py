import ctypes
lib = ctypes.cdll.LoadLibrary('/usr/users/iff_th2/isele/Applications/Scripts/MSR_beating_rod/libperformance_tools.so')

class Performance_tools(object):
    def __init__(self):
        self.obj = lib.Performance_tools_new()

    ############################################################################

    def compute_MSR(self, t_msr, time, phi, msr, counter, linval, nsteps, nmsr, limit, dt):

        # transfrom arrays to ctypes
        t_msrc = t_msr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        phic = phi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        msrc = msr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        counterc = counter.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        linvalc = linval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        nmsrc = ctypes.c_int(nmsr)
        limitc = ctypes.c_int(limit)
        dtc = ctypes.c_int(dt)

        lib.compute_msr(self.obj, t_msrc, timec, phic, msrc, counterc, linvalc, nstepsc, nmsrc, limitc, dtc)
        return

    ############################################################################

    def compute_MSR_with_derivatives(self, t_msr, time, phi, dphi, d2phi, msr, dmsr, d2msr, d2msr_t1, d2msr_t2, counter, linval, nsteps, nmsr, limit, dt):

        # transfrom arrays to ctypes
        t_msrc = t_msr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        phic = phi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        dphic = dphi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        d2phic = d2phi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        msrc = msr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        dmsrc = dmsr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        d2msrc = d2msr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        d2msr_t1c = d2msr_t1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        d2msr_t2c = d2msr_t2.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        counterc = counter.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        linvalc = linval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        nmsrc = ctypes.c_int(nmsr)
        limitc = ctypes.c_int(limit)
        dtc = ctypes.c_int(dt)

        lib.compute_msr_with_derivatives(self.obj, t_msrc, timec, phic, dphic, d2phic, msrc, dmsrc, d2msrc, d2msr_t1c, d2msr_t2c, counterc, linvalc, nstepsc, nmsrc, limitc, dtc)
        return

    ############################################################################

    def compute_autocorrelation(self, t_ac, time, x, ac, counter, linval, nsteps, nac, limit, dt):

        # transfrom arrays to ctypes
        t_acc = t_ac.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        xc = x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        acc = ac.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        counterc = counter.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        linvalc = linval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        nacc = ctypes.c_int(nac)
        limitc = ctypes.c_int(limit)
        dtc = ctypes.c_int(dt)

        lib.compute_autocorrelation(self.obj, t_acc, timec, xc, acc, counterc, linvalc, nstepsc, nacc, limitc, dtc)
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

    ############################################################################

    def compute_autocorrelation_with_errorbars(self, t_ac, time, x, ac, std, blockdata, block_std, block_uncert, linval, nsteps, nac, limit, dt):

        # transfrom arrays to ctypes
        t_acc = t_ac.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        xc = x.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        acc = ac.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        stdc = std.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        blockdatac = blockdata.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        block_stdc = block_std.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        block_uncertc = block_uncert.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        linvalc = linval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        nacc = ctypes.c_int(nac)
        limitc = ctypes.c_int(limit)
        dtc = ctypes.c_int(dt)

        lib.compute_autocorrelation_with_errorbars(self.obj, t_acc, timec, xc, acc, stdc, blockdatac, block_stdc, block_uncertc, linvalc, nstepsc, nacc, limitc, dtc)
        return

    ############################################################################

    def compute_MSR_with_errorbars(self, t_msr, time, phi, msr, std, blockdata, block_std, block_uncert, linval, nsteps, nmsr, limit, dt):

        # transfrom arrays to ctypes
        t_msrc = t_msr.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        phic = phi.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        msrc = msr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        stdc = std.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        blockdatac = blockdata.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        block_stdc = block_std.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        block_uncertc = block_uncert.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        linvalc = linval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        nmsrc = ctypes.c_int(nmsr)
        limitc = ctypes.c_int(limit)
        dtc = ctypes.c_int(dt)

        lib.compute_msr_with_errorbars(self.obj, t_msrc, timec, phic, msrc, stdc, blockdatac, block_stdc, block_uncertc, linvalc, nstepsc, nmsrc, limitc, dtc)
        return
