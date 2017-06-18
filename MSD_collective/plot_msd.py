
""" Plot number fluctuations of filaments,
    by fixing persistence length,
    and changing Pe"""

### example command line arguments: 

##############################################################################

import sys
sys.path.append('../Utility')

import argparse
import numpy as np
import os
import matplotlib as mpl
mpl.use('Agg', warn=False)
import matplotlib.pyplot as plt
import misc_tools
import read_write
import pandas as pd

##############################################################################

def get_args():
    """ get context arguments"""

    parser = argparse.ArgumentParser()
    parser.add_argument("-fl", "--folder", nargs="?", \
                        const='/local/duman/SIMULATIONS/many_polymers_5/density_0.2/', \
                        help="Folder containing data, as in /local/duman/SIMULATIONS/many_polymers_5/density_0.2/") 
    parser.add_argument("-af", "--analysisfile", nargs="?", const="NUMBER_FLUC/number_fluc.data", \
                        help="Address of the analysis file as in NUMBER_FLUC/number_fluc.data")     
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/RolfData/PLOTS/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/RolfData/PLOTS/") 
    parser.add_argument("-sf", "--savefolder", nargs="?", const="NUMBER_FLUC", \
                        help="Specific folder for saving, as in NUMBER_FLUC")      
    args = parser.parse_args()
    
    return args

##############################################################################

def read_msd_data(path):
    """ read specific analysis data"""
    
    if os.path.exists(path):
        data = pd.read_csv(path, sep='\t', skiprows=1, header=0)
        x = data['Timestep']
        y = data['MSD'] 
    else:
        x = 0.
        y = 0.
    
    return x, y

##############################################################################

def e2e_theoretical(xil):
    """ yield the analytical expression for the end-to-end vector,
    with Kraktky-Porod model"""
    
    return 2.*xil - 2.*(xil)**(2)*(1-np.exp(-1/xil))

##############################################################################

def plot_data(xp, yp, sims, savebase, savefolder, param_choice):
    """ plot the data,
    NOTE THAT xp and yp are dictionaries with keys as chosen parameter,
    and values as the x and y axis of the plot"""
    
    ### set general plot properties

    sim = sims[sims.keys()[0]]
    #downlim = -1
    #uplim = sim.lx/4.
    num_ticks = 5
    ax_len = 1.0                          # Length of one subplot square box
    ax_b = 0.0                            # Beginning/offset of the subplot in the box
    ax_sep = 0.0                          # Separation length between two subplots
    total_subplots_in_x = 1               # Total number of subplots    
    fig = plt.figure()
    subp = misc_tools.Subplots(fig, ax_len, ax_sep, ax_b, total_subplots_in_x) 
    ax0 = subp.addSubplot()
    
    ### save properties
    
    base = savebase + savefolder + '/'
    os.system("mkdir -p " + base)  
    savepath = base + '/' + savefolder + '_' + param_choice + '.pdf'
    
    keys = xp.keys()
    keys = np.sort(keys)
    colors = plt.cm.rainbow(np.linspace(0, 1, len(keys)))   
    
    for j, key in enumerate(keys):
        
        x = np.array(xp[key])
        y = np.array(yp[key])
    
        label = r'$\xi_{p}/L=$' + str(key)
        line0 = ax0.loglog(x/sim.tau_D, y/sim.length**2, \
                         linewidth=2.0, label=label, color=colors[j])
#        line1 = ax0.plot(x, yth, \
#                         linewidth=2.0, label='_nolegend_', color=colors[j])        
    
#    ax0.set_xscale('log')
#    ax0.set_yscale('log')
    
    ### title
    
#    ax0.set_title("$t/\\tau_{D}$ = " + "{0:.2f}".format(time/sim.tau_D) + \
#        ", $t/\\tau_{A}$ = " + "{0:.2f}".format(time/sim.tau_A), fontsize=30)
    
    ### labels

    ax0.set_xlabel(r'$t/\tau_{D}$', fontsize=40)
    ax0.set_ylabel(r'$\Delta r^{2}/L^2$', fontsize=40)

    ### limits

    #ax0.set_xlim((0.4, 1.05))
    #ax0.set_ylim((0.4, 1.05))
    
    ### ticks
    
    #ax0.xaxis.set_ticks(np.linspace(0, 15, num_ticks, endpoint=True))
    #ax0.yaxis.set_ticks(np.linspace(0, uplim, num_ticks, endpoint=True))
    ax0.tick_params(axis='both', which='major', labelsize=30)
    
    ### legend

    ax0.legend(bbox_to_anchor=(1.005, 0.,0.65, 1.), loc=2, borderaxespad=0., \
        prop={'size': 20}, mode="expand", frameon=False)
    
    ### save 
    
    plt.savefig(savepath, dpi=300, bbox_inches='tight', pad_inches=0.08)        
    fig.clf()                          
        
    return

##############################################################################

def main():
    
    args = get_args()
    fix_choice = 'fp'          # plot with fp as the legend
    param_choice = 'kappa'
    fix_value = 7.0
    x, y, sims = misc_tools.collect_multiple_data_2D(args.folder, args.analysisfile, 
                                                  read_msd_data, 
                                                  fix_choice, fix_value)
    plot_data(x, y, sims, args.savebase, args.savefolder, param_choice)

    return

##############################################################################

if __name__ == '__main__':
    main()
    
##############################################################################    
