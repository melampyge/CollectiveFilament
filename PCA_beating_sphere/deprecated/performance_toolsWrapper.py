import ctypes
lib = ctypes.cdll.LoadLibrary('/usr/users/iff_th2/isele/Applications/Scripts/PCA_beating_rod/libperformance_tools.so')

class Performance_tools(object):
    def __init__(self):
        self.obj = lib.Performance_tools_new()

    ############################################################################

    def compute_curvature(self, curv, xall, yall, nsteps, natoms):

        # flatten arrays
        curvf = curv.ravel()
        xallf = xall.ravel()
        yallf = yall.ravel()
        # transfrom arrays to ctypes
        curvc = curvf.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        xallc = xallf.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        yallc = yallf.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        natomsc = ctypes.c_int(natoms)

        lib.compute_curvature(self.obj, curvc, xallc, yallc, nstepsc, natomsc)
        return
        

    ############################################################################

    def compute_crosscorrelation(self, t_cc, time, a1, a2, cc1, cc2, counter, linval, nsteps, ncc, limit, dt):

        # transfrom arrays to ctypes
        t_ccc = t_cc.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        a1c = a1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        a2c = a2.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        cc1c = cc1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        cc2c = cc2.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        counterc = counter.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        linvalc = linval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        nccc = ctypes.c_int(ncc)
        limitc = ctypes.c_int(limit)
        dtc = ctypes.c_int(dt)

        lib.compute_crosscorrelation(self.obj, t_ccc, timec, a1c, a2c, cc1c, cc2c, counterc, linvalc, nstepsc, nccc, limitc, dtc)
        return

    ############################################################################

    def compute_crosscorrelation_with_errorbars(self, t_cc, time, x1, x2, cc, std, blockdata, block_std, block_uncert, linval, nsteps, ncc, limit, dt):

        # transfrom arrays to ctypes
        t_ccc = t_cc.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        x1c = x1.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        x2c = x2.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        ccc = cc.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        stdc = std.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        blockdatac = blockdata.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        block_stdc = block_std.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        block_uncertc = block_uncert.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        linvalc = linval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        nccc = ctypes.c_int(ncc)
        limitc = ctypes.c_int(limit)
        dtc = ctypes.c_int(dt)

        lib.compute_crosscorrelation_with_errorbars(self.obj, t_ccc, timec, x1c, x2c, ccc, stdc, blockdatac, block_stdc, block_uncertc, linvalc, nstepsc, nccc, limitc, dtc)
        return
