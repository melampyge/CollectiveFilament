
""" Plot shape parameter with respect to time,
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
    parser.add_argument("-af", "--analysisfile", nargs="?", const="RGYR/rgyr.data", \
                        help="Address of the analysis file as in RGYR/rgyr.data")     
    parser.add_argument("-sb", "--savebase", nargs="?", \
                        const = "/usr/users/iff_th2/duman/RolfData/PLOTS/", \
                        help="Folder to save the data, as in /usr/users/iff_th2/duman/RolfData/PLOTS/") 
    parser.add_argument("-sf", "--savefolder", nargs="?", const="RGYR", \
                        help="Specific folder for saving, as in RGYR")      
    args = parser.parse_args()
    
    return args

##############################################################################

def read_rgyr_data(path):
    """ read specific analysis data"""
    
    if os.path.exists(path):
        data = np.transpose(np.loadtxt(path, dtype=float))
        time = data[0]
        lamda_1 = data[1] 
        lamda_12 = data[2]
        lamda_2 = data[3]
        eig_1 = data[4]
        eig_2 = data[5]
    else:
        time = lamda_1 = lamda_12 = lamda_2 = eig_1 = eig_2 = 0.
    
    return time, lamda_1, lamda_12, lamda_2, eig_1, eig_2
    
##############################################################################

def unpack_data(data):
    """ unpack multidimensional array data into scalars"""

    return data[0], data[1], data[2], data[3], data[4], data[5]

##############################################################################

def compute_shape(x, y):
    """ compute the shape parameter"""
            
    return y/x

##############################################################################

def get_shape_param(data):
    """ calculate the shape parameter"""
    
    x = {}
    y = {}
    for key in data.keys():
        t, lamda_1, lamda_12, lamda_2, eig_1, eig_2 = unpack_data(data[key])
        x[key] = t
        y[key] = compute_shape(eig_1, eig_2)
        
    return x, y

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
    
        label = r'$Pe=$' + str(key)
        line0 = ax0.loglog(x/sim.tau_D, y, \
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
    ax0.set_ylabel(r'$\lambda_{1}/\lambda_{2}$', fontsize=40)

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
    fix_choice = 'kappa'          # plot with fp as the legend
    param_choice = 'fp'
    fix_value = 200.0
    data, sims = misc_tools.collect_multiple_data_MultiD(args.folder, args.analysisfile, 
                                                         read_rgyr_data, 
                                                         fix_choice, fix_value)
    x, y = get_shape_param(data)
    plot_data(x, y, sims, args.savebase, args.savefolder, param_choice)

    return

##############################################################################

if __name__ == '__main__':
    main()
    
##############################################################################    
