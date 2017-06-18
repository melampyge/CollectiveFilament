
""" Helper functions for miscelleanous operations"""

##############################################################################

import read_write
import numpy as np

##############################################################################

class Simulation:
    """ container for general simulation information"""
    
    def __init__(self, dt, ti, lx, ly, totalStep, nsamp, nbpf, \
                 nbeads, bl, kT, fmc, kbend):
        
        self.dt = dt*nsamp
        self.ti = ti
        self.lx = lx
        self.ly = ly
        self.nsteps = (totalStep - ti)/nsamp
        self.nbpf = nbpf
        self.nbeads = nbeads
        self.bl = bl
        self.kT = kT
        self.fp = fmc
        self.kappa = kbend
        self.nfils = self.nbeads/self.nbeads
        self.length = bl*(self.nbpf-1)
        self.area = self.lx*self.ly
        self.pe = self.fp*self.length**2/self.kT
        self.xil = self.kappa/self.kT/self.length
        self.flex = self.pe/self.xil
        self.tau_D = self.length**2*(self.nbpf+1)/4./self.kT
        self.tau_A = 0.
        if self.fp != 0.:
            self.tau_A = (self.nbpf+1)/self.fp
        
        return

##############################################################################

class Subplots:
    """ plot structure"""
    
    totcnt = -1             # Total number of subplots 
    
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction
        
    def addSubplot(self):
        """ add a subplot in the grid structure"""
        
        ## increase the number of subplots in the figure
        
        self.totcnt += 1
        
        ## get indices of the subplot in the figure
        
        self.nx = self.totcnt%(self.tot)
        self.ny = self.totcnt/(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])

##############################################################################

def gen_folders(base, fl, d, k, f):
    """ generate paths of folders to the analysis files per parameter"""
    
    datafolder = base + 'density_' + str(d) + '/kappa_' + str(k) + \
        '/fp_' + str(f) + '/'
    datafile = datafolder + 'init_info.txt'
    analysisfile = datafolder + fl
    
    return datafile, analysisfile

##############################################################################

def collect_data(basefolder, analysisfilepath, read_fnc, param_choice):
    """ collect the analysis data -SINGULAR- as a function of both of the parameters"""
    
    density = 0.2
    kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 300.0, 400.0]
    fp = [0.0, 0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
    data = {}
    sims = {}
    for k in kappa:
        for f in fp:
            key = (k, f)
            datafile, analysisfile = gen_folders(basefolder, analysisfilepath, 
                                                   density, k, f)
            sims[key] = read_write.read_sim_info(datafile)
            data[key] = read_fnc(analysisfile)
                
    if param_choice == 'kappa':
        xp, yp = accumulate_data(kappa, fp, 'kappa', data, sims)
    elif param_choice == 'fp':
        xp, yp = accumulate_data(fp, kappa, 'fp', data, sims)

        
    return xp, yp, sims

##############################################################################

def collect_multiple_data(basefolder, analysisfilepath, read_fnc, fix_choice, fix_value):
    """ collect the analysis data -MULTIPLE- as a function of both of the parameters"""
    
    density = 0.2
    kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 300.0, 400.0]
    fp = [0.0, 0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    #fp = [0.0, 0.0112, 0.024, 0.08, 1.2, 2.4, 7.0]
    
    ### param is the value to plot as a function of
    
    if fix_choice == 'kappa':
        param = fp
    elif fix_choice == 'fp':
        param = kappa
    
    xp = {}
    yp = {}
    ystdp = {}
    sims = {}
    for p in param:
        if fix_choice == 'kappa':
            datafile, analysisfile = gen_folders(basefolder, analysisfilepath, 
                                                 density, fix_value, p)   
            key = sims
        elif fix_choice == 'fp':
            datafile, analysisfile = gen_folders(basefolder, analysisfilepath, 
                                                 density, p, fix_value)            
        sim = read_write.read_sim_info(datafile)
        if fix_choice == 'kappa':
            key = sim.pe
        elif fix_choice == 'fp':
            key = sim.xil
        sims[key] = sim
        data = read_fnc(analysisfile)
        xp[key] = data[0]
        yp[key] = data[1]
        ystdp[key] = data[2]
    
    return xp, yp, sims                      
    #return xp, yp, ystdp, sims

##############################################################################

def accumulate_data(param, other_param, param_choice, data, sims):
    """ acccumulate data -SINGULAR- as a function of the chosen parameter"""
    
    xp = {}
    yp = {}
    for p in param:
        x = []
        y = []
        for r in other_param: 
            if param_choice == 'kappa':
                key = (p, r) 
            elif param_choice == 'fp':
                key = (r, p)
            x.append(sims[key].pe)
            y.append(data[key])
        new_key = sims[key].xil
        xp[new_key] = x
        yp[new_key] = y
    
    return xp, yp

##############################################################################
