import ctypes
lib = ctypes.cdll.LoadLibrary('/usr/users/iff_th2/duman/RolfData/Scripts/MSD_collective/libperformance_tools.so')

class Performance_tools(object):
    def __init__(self):
        self.obj = lib.Performance_tools_new()

    ########################################################################

    def correct_pbc(self, comx, dx, lx, nsteps, nmol):
        # flatten arrays if required
        comxf = comx.ravel()
        # transfrom arrays to ctypes
        comxc = comxf.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        dxc = dx.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        # transform integers and doubles to ctypes
        nstepsc = ctypes.c_int(nsteps)
        lxc = ctypes.c_double(lx)
        nmolc = ctypes.c_int(nmol)

        lib.correct_pbc(self.obj, comxc, dxc, lxc, nstepsc, nmolc)
        return
        
    ########################################################################

    def compute_MSD(self, t_msd, time, comx, comy, msd_tmp, counter, logval, nsteps, nmsd, limit, nevery):

        # transfrom arrays to ctypes
        t_msdc = t_msd.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        timec = time.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        comxc = comx.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        comyc = comy.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        msd_tmpc = msd_tmp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
        counterc = counter.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        logvalc = logval.ctypes.data_as(ctypes.POINTER(ctypes.c_int))
        # transform integers and double to ctypes
        nstepsc = ctypes.c_int(nsteps)
        nmsdc = ctypes.c_int(nmsd)
        limitc = ctypes.c_int(limit)
        neveryc = ctypes.c_int(nevery)

        lib.compute_msd(self.obj, t_msdc, timec, comxc, comyc, msd_tmpc, counterc, logvalc, nstepsc, nmsdc, limitc, neveryc)
        return
        
