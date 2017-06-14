
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import codecs
import read_char
import h5py
import os

try:
    infilename = sys.argv[1]
except:
    print 'Usage: ' + sys.argv[0] + '       parameter file'
    exit()
    
##################################################################

def read_settings():
    """ read in the settings from the parameter file"""
    
    ## open file for reading
    ifile = open(infilename, 'r')
    
    ## skip comment line
    ifile.readline()
    
    ## char-file name
    line = ifile.readline()
    line = line.split()
    charfile = line[-1]
    
    ## header-file name
    line = ifile.readline()
    line = line.split()
    headerfile = line[-1]
    
    ## output folder
    line = ifile.readline()
    line = line.split()
    ofname = line[-1]
    
    ## length of the filaments
    line = ifile.readline()
    line = line.split()
    nfil = int(line[-1])
    
    ## number of snapshots to skip
    line = ifile.readline()
    line = line.split()
    nskip = int(line[-1])
    
    ## close file
    ifile.close()
    
    return charfile, headerfile, ofname, nfil, nskip

##################################################################

class Subplots:
    """ Arrange subplot grid structure (square box is assumed)"""
    
    totcnt = -1             # Total number of subplots 
    
    ## constructor
    def __init__(self, f, l, s, b, t):
        self.fig = f        # Figure axes handle
        self.length = l     # Length of the subplot box 
        self.sep = s        # Separation distance between subplots 
        self.beg = b        # Beginning (offset) in the figure box
        self.tot = t        # Total number of subplots in the x direction
        
    ## add a subplot in the grid structure
    def addSubplot(self):
        
        # Increase the number of subplots in the figure
        self.totcnt += 1
        
        # Indices of the subplot in the figure
        self.nx = self.totcnt%(self.tot)
        self.ny = self.totcnt/(self.tot)
        
        self.xbeg = self.beg + self.nx*self.length + self.nx*self.sep
        self.ybeg = self.beg + self.ny*self.length + self.ny*self.sep
        return self.fig.add_axes([self.xbeg,self.ybeg,self.length,self.length])
        
##################################################################

def nearest_neighbor(x1,x3,lx):
    """ compute vector to nearest neighbor"""
    
    dx1 = x3 - x1
    dx2 = x3 - x1 + lx
    dx3 = x3 - x1 - lx
    if dx1**2 < dx2**2 and dx1**2 < dx3**2:
        return dx1
    if dx2**2 < dx3**2:
        return dx2
    return dx3
    
##################################################################
    
def compute_velocity(x1,y1,x3,y3,lx,ly,tstep1,tstep3,natoms):
    """ compute the velocity of each bead, consider pbc"""
    
    ## define target arrays
    vx = np.zeros(natoms)
    vy = np.zeros(natoms)
    
    ## loop over all beads
    for i in range(natoms):
        xi = x1[i]
        yi = y1[i]
        xj = x3[i]
        yj = y3[i]
        dx = nearest_neighbor(xi,xj,lx)
        dy = nearest_neighbor(yi,yj,ly)
        vx[i] = dx/(tstep3-tstep1)
        vy[i] = dy/(tstep3-tstep1)
    return vx, vy

##################################################################

def bin_data(x,y,vx,vy,natoms,lx,ly,bx,by):
    """ create binned densities and velocities"""
    
    ## allocate output arrays
    rho = np.zeros((bx,by))
    vx_b = np.zeros((bx,by))
    vy_b = np.zeros((bx,by))
    px_b = np.zeros((bx,by))
    py_b = np.zeros((bx,by))
    
    ## fill all atoms into bins
    for i in range(natoms):
        ## get coordinates
        xi = x[i]
        yi = y[i]
        ## get current bin
        segx = int(xi/lx*bx)
        segy = int(yi/ly*by)
        ## add data to bin
        rho[segx,segy] += 1
        px_b[segx,segy] += vx[i]
        py_b[segx,segy] += vy[i]
        
    ## transform moments to velocities
    for i in range(bx):
        for j in range(by):
            if rho[i,j] > 1:
                vx_b[i,j] = px_b[i,j]/rho[i,j]
                vy_b[i,j] = py_b[i,j]/rho[i,j]
                
    ## transform number counts to densities
    wx = lx/bx
    wy = ly/by
    rho /= wx*wy
    
    ## generate an array that contains the edges
    xedges = np.linspace(0., 1., bx+1)*lx
    yedges = np.linspace(0., 1., by+1)*ly
    
    return xedges, yedges, rho, vx_b, vy_b

##################################################################

def compute_rotation(vx,vy,bx,by,lx,ly):
    """ compute the rotation of a vector field"""
    
    ## allocate output array
    rotation = np.zeros((bx,by))
    
    ## bin width
    wx = lx/bx
    wy = ly/by
    
    ## compute rotation in each bin
    for i in range(bx):
        for j in range(by):
            
            ## compute velocity gradients using finite differences
            duy_dx = (vy[(i+1)%bx,j] - vy[i-1,j])/(2*wx)
            dux_dy = (vx[i,(j+1)%by] - vx[i,j-1])/(2*wy)
            rotation[i,j] = duy_dx - dux_dy
            
    return rotation
    
##################################################################    
    
def analyse_spirals(charfile, headerfile, nfil, nskip):
    """ analyse the spiral structure"""
    
    ## open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    
    ## get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
    
    ## skip initial snapshots and correct nsteps value
    read_char.skip_snapshots(hfile, ifile, nskip)
    nsteps = nsteps - 2 - nskip

    ## read in the first two steps
    xs, ys, lx, ly, tstep2, natoms = read_char.read_snapshot(hfile, ifile)
    x2 = xs*lx
    y2 = ys*ly
    xs, ys, lx, ly, tstep3, natoms = read_char.read_snapshot(hfile, ifile)
    x3 = xs*lx
    y3 = ys*ly
    
    ## bin width and number of bins
    bwx = 1.
    bwy = 1.
    bx = int(np.ceil(lx/bwx))
    by = int(np.ceil(ly/bwy))

    ## allocate arrays
    time = np.zeros((nsteps), dtype = int)
    rho_hist = np.zeros((nsteps,bx,by))
    vx_hist = np.zeros((nsteps,bx,by))
    vy_hist = np.zeros((nsteps,bx,by))
    vorticity = np.zeros((nsteps,bx,by))
    
    ## loop over all steps
    for i in range(nsteps):
        
        ## print some information
        print 'Process:', i+1, '/', nsteps
        
        ## move previous data to other arrays
        x1 = np.copy(x2)
        y1 = np.copy(y2)
        tstep1 = tstep2
        x2 = np.copy(x3)
        y2 = np.copy(y3)
        tstep2 = tstep3
        
        ## read in the new coordinates
        xs,ys,lx,ly,tstep3,natoms = read_char.read_snapshot(hfile, ifile)
        x3 = xs*lx
        y3 = ys*ly
        
        ## approximate the velocity of the particles (finite difference scheme)
        vx,vy = compute_velocity(x1,y1,x3,y3,lx,ly,tstep1,tstep3,natoms)

        ## bin densities and velocities
        xedges, yedges, rho_hist[i], vx_hist[i], vy_hist[i] = \
            bin_data(x2,y2,vx,vy,natoms,lx,ly,bx,by)
        
        ## compute vorticity
        vorticity[i] = compute_rotation(vx_hist[i],vy_hist[i],bx,by,lx,ly)

        ## add current step to time array
        time[i] = tstep2
            
                    
    ## close the input files
    hfile.close()
    ifile.close()

    return nsteps, lx, ly, time, xedges, yedges, rho_hist, vx_hist, vy_hist, vorticity  

##################################################################

def save_data(ofile, nsteps, lx, ly, time, xedges, yedges, rho_hist, vx_hist, vy_hist, vorticity):
    """ save data in hdf5 format"""   
    
    edges_grp = ofile.create_group('edges')
    edges_grp.create_dataset('x', data=xedges, compression='gzip')
    edges_grp.create_dataset('y', data=yedges, compression='gzip')
    
    ofile.create_dataset('time', data=time, compression='gzip')
    
    tables_grp = ofile.create_group('tables')
    tables_grp.create_dataset('rho', data=rho_hist, compression='gzip')
    tables_grp.create_dataset('vx', data=vx_hist, compression='gzip')
    tables_grp.create_dataset('vy', data=vy_hist, compression='gzip')
    tables_grp.create_dataset('vorticity', data=vorticity, compression='gzip')
    
    box_grp = ofile.create_group('box')
    box_grp.create_dataset('x', data=lx)
    box_grp.create_dataset('y', data=ly)
    
    ofile.create_dataset('nsteps', data=nsteps)
    
    return
    
##################################################################

def load_data(infile):
    """ save data in hdf5 format"""   
    
    f = h5py.File(infile, 'r')
    
    edges_grp = f['edges']
    xedges = np.asarray(edges_grp['x'], dtype=float)
    yedges = np.asarray(edges_grp['y'], dtype=float)

    time = np.asarray(f['time'])

    tables_grp = f['tables']
    rho_hist = np.asarray(tables_grp['rho'], dtype=float)
    vx_hist = np.asarray(tables_grp['vx'], dtype=float)
    vy_hist = np.asarray(tables_grp['vy'], dtype=float)
    vorticity = np.asarray(tables_grp['vorticity'], dtype=float)  
    
    box_grp = f['box']
    lx = box_grp['x'][...]
    ly = box_grp['y'][...]
    
    nsteps = len(time)
    
    #nsteps = f['nsteps'][...]
    f.close()

    return nsteps, lx, ly, time, xedges, yedges, rho_hist, vx_hist, vy_hist, vorticity   

##################################################################

def plot_data(sfolder, nsteps, lx, ly, time, xedges, yedges, rho_hist, vx_hist, vy_hist, vorticity):
    """ plot data"""   
    
    xedges *= 2.
    yedges *= 2.
    xgrid, ygrid = np.meshgrid(xedges, yedges)
    
    fig = plt.figure()
    downlim = 0
    uplim = lx*2.
    ax_len = 0.5                         # Length of one subplot square box
    ax_b = 0.1                            # Beginning/offset of the subplot in the box
    ax_sep = 0.15                         # Separation length between two subplots
    total_subplots_in_x = 2               # Total number of subplots in the x direction  
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x)  
    
    print np.shape(xgrid)
    print np.shape(vx_hist[0])
    print np.shape(vorticity)
    print np.shape(vorticity[0])
    for i in range(nsteps):
        
        print i, ' of ', nsteps
#        vxgrid, vygrid = np.meshgrid(np.transpose(vx_hist[i]), np.transpose(vy_hist[i]))

        ax = subp.addSubplot()

        line = ax.pcolor(xgrid, ygrid, vorticity[i], cmap='jet', alpha=0.7)
#        ax.quiver(xedges, yedges, vx_hist[i], vy_hist[i])
#        q = ax.quiver(xgrid, ygrid, vxgrid, vygrid)
        ax.set_xlim((downlim,uplim))
        ax.set_ylim((downlim,uplim))
        ax.xaxis.set_ticks( np.linspace(downlim, uplim, num=4, endpoint=True) )
        ax.yaxis.set_ticks( np.linspace(downlim, uplim, num=4, endpoint=True) )
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_xlabel('x/b', fontsize=25)
        ax.set_ylabel('y/b', fontsize=25)
        #ax.set_title('Vorticity', fontsize=25)
        
        cax = plt.axes([subp.xbeg+ax_len+0.01, subp.ybeg, 0.02, ax_len])
        plt.colorbar(line, cax=cax)
        #cax.set_yticks([])
        #cax.set_xticks([])
        cax.tick_params(labelsize=5)
        cax.set_title('w',fontsize=10)
        
        timestep = time[i]
        
        plt.savefig(sfolder + '/frame-'+'{0:05d}'.format(int(timestep))+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
        plt.clf()
    
    return
    
##################################################################

def extended_plot_data(sfolder, nsteps, lx, ly, time, xedges, yedges, \
    rho_hist, vx_hist, vy_hist, vorticity, nfil, charfile, headerfile, nskip):
    """ plot data"""   
    
        ## open files for reading
    ifile = codecs.open(charfile, 'r', 'UTF-8')
    hfile = open(headerfile, 'r')
    
    ## get general information from header file
    natoms, nsteps = read_char.read_first(hfile)
    nmol = natoms/nfil
    
    ## skip initial snapshots and correct nsteps value
    read_char.skip_snapshots(hfile, ifile, nskip)

    ## read in the first two steps
    xs, ys, lx, ly, tstep2, natoms = read_char.read_snapshot(hfile, ifile)
    x2 = xs*lx
    y2 = ys*ly
    xs, ys, lx, ly, tstep3, natoms = read_char.read_snapshot(hfile, ifile)
    x3 = xs*lx
    y3 = ys*ly
    
    xedges *= 2.
    yedges *= 2.
    xgrid, ygrid = np.meshgrid(xedges, yedges)
    
    fig = plt.figure()
    downlim = 0
    uplim = lx*2.
    ax_len = 0.5                         # Length of one subplot square box
    ax_b = 0.1                            # Beginning/offset of the subplot in the box
    ax_sep = 0.15                         # Separation length between two subplots
    total_subplots_in_x = 2               # Total number of subplots in the x direction  
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x)  
    
    print np.shape(xgrid)
    print np.shape(vx_hist[0])
    print np.shape(vorticity)
    print np.shape(vorticity[0])
    
    ## loop over all steps
    for i in range(nsteps):
        
        ## print some information
        print 'Process:', i+1, '/', nsteps
        
        ## move previous data to other arrays
        x1 = np.copy(x2)
        y1 = np.copy(y2)
        tstep1 = tstep2
        x2 = np.copy(x3)
        y2 = np.copy(y3)
        tstep2 = tstep3
        
        ## read in the new coordinates
        xs,ys,lx,ly,tstep3,natoms = read_char.read_snapshot(hfile, ifile)
        x3 = xs*lx
        y3 = ys*ly
        
        ## approximate the velocity of the particles (finite difference scheme)
        vx,vy = compute_velocity(x1,y1,x3,y3,lx,ly,tstep1,tstep3,natoms)
           
        ax = subp.addSubplot()
        
        print tstep1, tstep2, time[i]
        
        line = ax.pcolor(xgrid, ygrid, np.transpose(vorticity[i]), cmap='cool', alpha=0.7)
#        ax.quiver(x1*2.,y1*2.,vx*2.,vy*2.)
        ax.set_xlim((downlim,uplim))
        ax.set_ylim((downlim,uplim))
        ax.xaxis.set_ticks( np.linspace(downlim, uplim, num=4, endpoint=True) )
        ax.yaxis.set_ticks( np.linspace(downlim, uplim, num=4, endpoint=True) )
        ax.tick_params(axis='both', which='major', labelsize=12)
        ax.set_xlabel('x/r', fontsize=25)
        ax.set_ylabel('y/r', fontsize=25)
        #ax.set_title('Vorticity', fontsize=25)
        
        cax = plt.axes([subp.xbeg+ax_len+0.01, subp.ybeg, 0.02, ax_len])
        plt.colorbar(line, cax=cax)
        #cax.set_yticks([])
        #cax.set_xticks([])
        cax.tick_params(labelsize=5)
        cax.set_title('w',fontsize=10)
        
        timestep = time[i]
        
        plt.savefig(sfolder + '/frame-'+'{0:05d}'.format(int(timestep))+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
        plt.clf()           
                    
    ## close the input files
    hfile.close()
    ifile.close()
    
    return    
        
##################################################################

def main():
    """ main function, called when the script is started"""
    
    ## read parameters from input file
    charfile, headerfile, ofname, nfil, nskip = read_settings()

#    ## analyse the spiral structure 
#    nsteps, lx, ly, time, xedges, yedges, rho_hist, vx_hist, vy_hist, vorticity = \
#        analyse_spirals(charfile, headerfile, nfil, nskip)
#        
#    ## save data
#    os.system('mkdir -p ' + ofname)
#    ofile = h5py.File(ofname + '/data.hdf5', 'w')
#    save_data(ofile, nsteps, lx, ly, time, xedges, yedges, rho_hist, vx_hist, vy_hist, vorticity)
#    ofile.close()
    
    ## load data 
    infile = ofname + '/data.hdf5'
    nsteps, lx, ly, time, xedges, yedges, rho_hist, vx_hist, vy_hist, vorticity \
        = load_data(infile)
        
        
    ## plot data
    savebase = "/usr/users/iff_th2/duman/RolfData/long_filaments"
    sfolder = savebase + "/" + sys.argv[2] + "/" + ofname
    os.system('mkdir -p ' + sfolder)
    extended_plot_data(sfolder, nsteps, lx, ly, time, xedges, yedges, \
        rho_hist, vx_hist, vy_hist, vorticity, nfil, charfile, headerfile, nskip)



if __name__ == '__main__':
    main()
        