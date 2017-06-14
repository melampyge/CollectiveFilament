import ctypes
lib = ctypes.cdll.LoadLibrary('/usr/users/iff_th2/isele/Applications/Scripts/MSD_beating_rod/libperformance_tools.so')

class Performance_tools(object):
    def __init__(self):
        self.obj = lib.Performance_tools_new()

    ############################################################################

    def compute_MSD(self, t_msd, time, comx, comy, msd, counter, linval, nsteps, nmsd, limit, dt):

        # transfrom arrays to ctypes
        t_msdc = t_msd.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        comxc = comx.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        comyc = comy.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        msdc = msd.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        counterc = counter.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        linvalc = linval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        nmsdc = ctypes.c_int(nmsd)
        limitc = ctypes.c_int(limit)
        dtc = ctypes.c_int(dt)

        lib.compute_msd(self.obj, t_msdc, timec, comxc, comyc, msdc, counterc, linvalc, nstepsc, nmsdc, limitc, dtc)
        return

    ############################################################################

    def compute_MSD_with_errorbars(self, t_msd, time, comx, comy, msd, std, blockdata, block_std, block_uncert, linval, nsteps, nmsd, limit, dt):

        # transfrom arrays to ctypes
        t_msdc = t_msd.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        comxc = comx.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        comyc = comy.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        msdc = msd.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        stdc = std.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        blockdatac = blockdata.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        block_stdc = block_std.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        block_uncertc = block_uncert.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        linvalc = linval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        nmsdc = ctypes.c_int(nmsd)
        limitc = ctypes.c_int(limit)
        dtc = ctypes.c_int(dt)

        lib.compute_msd_with_errorbars(self.obj, t_msdc, timec, comxc, comyc, msdc, stdc, blockdatac, block_stdc, block_uncertc, linvalc, nstepsc, nmsdc, limitc, dtc)
        return

