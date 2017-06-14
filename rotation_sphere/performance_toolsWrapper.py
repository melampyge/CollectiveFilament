import ctypes
lib = ctypes.cdll.LoadLibrary('/usr/users/iff_th2/isele/Applications/Scripts/rotation_sphere/libperformance_tools.so')

class Performance_tools(object):
    def __init__(self):
        self.obj = lib.Performance_tools_new()

    ############################################################################

    def compute_ttac_with_errorbars(self, t_ac, time, tx, ty, ac, std, blockdata, block_std, block_uncert, linval, nsteps, nac, limit, dt):

        # transfrom arrays to ctypes
        t_acc = t_ac.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        txc = tx.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        tyc = ty.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
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

        lib.compute_ttac_with_errorbars(self.obj, t_acc, timec, txc, tyc, acc, stdc, blockdatac, block_stdc, block_uncertc, linvalc, nstepsc, nacc, limitc, dtc)
        return
