
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import codecs
import read_char
import h5py
import os
import argparse

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

def load_data(infile, nstep):
    """ save data in hdf5 format"""   
    
    f = h5py.File(infile, 'r')
    
    edges_grp = f['edges']
    xedges = np.asarray(edges_grp['x'][nstep], dtype=float)
    yedges = np.asarray(edges_grp['y'][nstep], dtype=float)

    time = np.asarray(f['time'][nstep])

    tables_grp = f['tables']
    rho_hist = np.asarray(tables_grp['rho'][nstep], dtype=float)
    vx_hist = np.asarray(tables_grp['vx'][nstep], dtype=float)
    vy_hist = np.asarray(tables_grp['vy'][nstep], dtype=float)
    vorticity = np.asarray(tables_grp['vorticity'][nstep], dtype=float)  
    
    box_grp = f['box']
    lx = box_grp['x'][...]
    ly = box_grp['y'][...]
        
    #nsteps = f['nsteps'][...]
    f.close()

    return lx, ly, time, xedges, yedges, rho_hist, vx_hist, vy_hist, vorticity   

##################################################################

def plot_data(sfolder, lx, ly, time, xedges, yedges, rho_hist, vx_hist, vy_hist, vorticity, save):
    """ plot data"""   
    
    xedges *= 2.
    yedges *= 2.
    print xedges
    xgrid, ygrid = np.meshgrid(xedges, yedges)
    
    fig = plt.figure()
    downlim = 0
    uplim = lx*2.
    ax_len = 0.5                         # Length of one subplot square box
    ax_b = 0.1                            # Beginning/offset of the subplot in the box
    ax_sep = 0.15                         # Separation length between two subplots
    total_subplots_in_x = 2               # Total number of subplots in the x direction  
    subp = Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x)  
    
    
    vxgrid, vygrid = np.meshgrid(np.transpose(vx_hist), np.transpose(vy_hist))

    ax = subp.addSubplot()

    #line = ax.pcolor(xgrid, ygrid, np.transpose(vorticity[i]), cmap='jet', alpha=0.7)
    #ax.quiver(xedges, yedges, vx_hist[i], vy_hist[i])
    q = ax.quiver(xgrid, ygrid, vxgrid, vygrid)
#        ax.set_xlim((downlim,uplim))
#        ax.set_ylim((downlim,uplim))
#        ax.xaxis.set_ticks( np.linspace(downlim, uplim, num=4, endpoint=True) )
#        ax.yaxis.set_ticks( np.linspace(downlim, uplim, num=4, endpoint=True) )
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.set_xlabel('x/r', fontsize=25)
    ax.set_ylabel('y/r', fontsize=25)
    #ax.set_title('Vorticity', fontsize=25)
    
#        cax = plt.axes([subp.xbeg+ax_len+0.01, subp.ybeg, 0.02, ax_len])
#        plt.colorbar(line, cax=cax)
#        #cax.set_yticks([])
#        #cax.set_xticks([])
#        cax.tick_params(labelsize=5)
#        cax.set_title('w',fontsize=10)
        
    if save:
        plt.savefig(sfolder + '/frame-'+'{0:05d}'.format(int(time))+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
        plt.savefig(sfolder + '/frame-'+'{0:05d}'.format(int(time))+'.eps',dpi=200,bbox_inches='tight',pad_inches=0.08)
    else:
        plt.savefig(sfolder + '/frame-'+'{0:05d}'.format(int(time))+'.png',dpi=200,bbox_inches='tight',pad_inches=0.08)
        
    plt.clf()
    
    return
        
##################################################################

def main():
    """ main function, called when the script is started"""
    
    ## read parameters from command line
    parser = argparse.ArgumentParser()    
    parser.add_argument("-d", "--density", help="Density")
    parser.add_argument("-k", "--kappa", help="Bending rigidity")
    parser.add_argument("-t", "--time", help="Timestep")    
    parser.add_argument("-s", "--save", help="Save options", action="store_true")      
    args = parser.parse_args()

    ## load data 
    database = "/local/duman/SIMULATIONS/long_filaments"
    datafolder = database + "/density_" + args.density + "/kappa_" + args.kappa
    infile =  datafolder + "/VORTEX/data.hdf5"
    lx, ly, time, xedges, yedges, rho_hist, vx_hist, vy_hist, vorticity \
        = load_data(infile, int(args.time))
        
    ## plot data
    savebase = "~/RolfData/many_filaments_5"        
    if args.save:
        sfolder = savebase + "/plots/VORTEX" 
        os.system('mkdir -p ' + sfolder)
    else:
        sfolder = "~/Desktop"
    plot_data(sfolder, lx, ly, time, xedges, yedges, rho_hist, vx_hist, vy_hist, vorticity, args.save)



if __name__ == '__main__':
    main()
        