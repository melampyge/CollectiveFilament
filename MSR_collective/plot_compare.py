
## Load needed libraries and necessary files
import argparse
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os.path
import glob
import pandas as pd
from string import atof
import sys
from scipy.optimize import curve_fit

##########################################################################

#
# Function definitions
#

##########################################################################
 
def loadSimData(datafile):
    """ Load initial simulation data"""
    
    global dt, ti, Lx, Ly, nsamp, N, M, L, B, totalStep, Fmc, Kbend, kT, \
        dtSamp, T, box_area, nt, body_length, Pe, persistence 

    datafile = open(datafile,"r")
    for line in datafile:
        A = line.split()
        if A[0] == "dt":                    # Time interval between MD steps
            dt = float(A[-1])
        elif A[0] == "ti":                  # Beginning time for data acquisition
            ti = float(A[-1])
        elif A[0] == "Lx":                  # Box size in x
            Lx = float(A[-1])            
        elif A[0] == "Ly":                  # Box size in y
            Ly = float(A[-1])
        elif A[0] == "totalStep":           # Total MD steps
            totalStep = float(A[-1])
        elif A[0] == "nsamp":               # Data sampling frequency
            nsamp = float(A[-1])
        elif A[0] == "nfil":                # Number of particles per polymer
            N = float(A[-1])
        elif A[0] == "L":                   # Number of particles
            L = float(A[-1])
        elif A[0] == "B":                   # Bond length between particles of a body
            B = float(A[-1])
        elif A[0] == "kT":                  # Boltzmann constant*Temperature
            kT = float(A[-1])
        elif A[0] == "Fmc":                 # Self propulsion force constant
            Fmc = float(A[-1])     
        elif A[0] == "Kbend":               # Bending constant
            Kbend = float(A[-1])
    
    Lx /= B
    Ly /= B
    M = L/N
    dtSamp = dt*nsamp
    T = totalStep - ti
    nt = T/nsamp
    box_area = Lx*Ly
    body_length = B*N
    Pe = Fmc*body_length**2/kT
    persistence = Kbend/(kT*body_length)
    
    return

##########################################################################

def analyticFunc(t,a):
    """ Theoretical MSR"""
    return 2*a*t
    
##########################################################################
    
def linearFunc(x,a):
    """ Linear function"""
    return a * x
    
##########################################################################    
    
def quadFunc(x,a):
    """ Quadratic function"""
    return a * x**2  

##########################################################################    

def read_data(datafile):
    """ Read the data with black magic"""
    
    data = pd.read_csv(datafile, sep='\t', skiprows=1, header=0) 
    x = data['Timestep']*dt
    y = data['MSR']

    return x, y
   
##########################################################################    
    
def fit_data(x,y):
    """ Fit the y data"""
    
    ## Curve fitting
    popt, pcov = curve_fit(analyticFunc, x, y)
    
    ## Fitted function
    yfit = analyticFunc(x,popt[0])
    
    return yfit, popt[0]
    
##########################################################################       

def save_data(f, out_fname, x, y, yfit):
    """ Save the data for plotting it better in the future"""
    
    ## create the folder in which saving is gonna be performed
    out_fname = '/usr/users/iff_th2/duman/RolfData/GraphData/' + out_fname
    if os.path.exists(out_fname) == False:
        os.mkdir(out_fname)
    out_fname += '/'
        
    ## parse the filename
    fp = f.split('/')[-1]
    ka = f.split('/')[-2]
    de = f.split('/')[-3]
    
    ## save the data into the designated file
    save_path = out_fname + 'density_' + de + \
        '_kappa_' + ka + '_fp_' + fp + '.data'
    save_file = open(save_path, 'w')
    
    for jj in range(len(x)):
        save_file.write(str(x[jj]) + '\t\t' + str(y[jj]) + '\t\t' + str(yfit[jj]) +'\n')
        
    save_file.close()    
 
    
##########################################################################

def save_fit_vars(f, out_fname, *args):
    """ Save fit variables for plotting it better in the future"""
    
    ## create the folder in which saving is gonna be performed
    out_fname = '/usr/users/iff_th2/duman/RolfData/GraphData/' + out_fname
    if os.path.exists(out_fname) == False:
        os.mkdir(out_fname)
    out_fname += '/'
        
    ## parse the filename
    fp = f.split('/')[-1]
    ka = f.split('/')[-2]
    de = f.split('/')[-3]
    
    ## save the data into the designated file
    save_path = out_fname + 'density_' + de + \
        '_kappa_' + ka + '_fp_' + fp + '.data'
    save_file = open(save_path, 'w')
    
    if len(args) == 1:
        save_file.write(args)
    elif len(args) == 2:
        save_file.write(args[0] + '\t\t' + args[1])
    elif len(args) == 3:
        save_file.write(args[0] + '\t\t' + args[1] + '\t\t' + args[2])
    elif len(args) == 4:
        save_file.write(args[0] + '\t\t' + args[1] + '\t\t' + args[2] + '\t\t' + args[3])
    else:
        print "More than 4 parameters are there in the fit here : ", f
        
    save_file.close()    
 
    
##########################################################################
    
def format_data(datafile, density, kappa, fp):
    """ Format the data"""

    ## Continue indexing the data 
    folders = []
    for d in density:
        for k in kappa:
            for f in fp:
                folders.append( datafile + '/density_' + str(d) + \
                    '/kappa_' + str(k) + '/fp_' + str(f) + '/' )
                    
    return folders    
       
##########################################################################
   
def main(): 
    
    ## Argument parsing (command line options)
    parser = argparse.ArgumentParser()
    parser.add_argument("datafolder", help="Folder containing data")
    parser.add_argument("analysis", help="Type of analysis")
    parser.add_argument("savefile", help="Folder in which figures should be saved")
    args = parser.parse_args()
    
    ## Index the data
    density = [0.2]
    kappa = [1.25, 2.5, 5.0, 25.0, 62.5, 125.0, 200.0, 400.0]
    fp = [0.0, 0.0048, 0.0112, 0.024, 0.08, 0.24, 0.8, 1.2, 2.4, 7.0]
    
    ## Format the data
    datafolder = "/local/duman/SIMULATIONS/" + args.datafolder
    folders = format_data(datafolder, density, kappa, fp)    
    
    ## Labels to be used in plotting 
    labels={'Pe':r"Pe",
       'xi_L':r"$\xi_p/L$",
       'density':r'\phi'}   
       
    ## Save the data
    x = []          # x data
    y = []          # y data
    yfit = []       # fit data
    yf1 = -1        # fit variable data
    yf2 = -1        # fit variable dataa
    yf3 = -1        # fit variable data
    yf4 = -1        # fit variable data
    for f in folders:
        
        ## Load initial simulation data
        initfile = datafolder + "init_info.txt"
        loadSimData(initfile)
        
        ## Load analysis data
        datafile = f + args.analysis
        if os.path.exists(datafile):
            
            ## Use black magic to read the text file
            data = read_data(datafile)
            if len(data) == 2:
                x = data[0]
                y = data[1]
            else:
                print "More than 2 arrays (corresponding to x and y) are here : ", datafile 
            
            
            ## Perform fitting and save the fit
            fit_data = fit_data(x,y)
            yfit = fit_data[0]
            if len(fit_data) > 1:
                if len(fit_data) == 2:
                    yf1 = fit_data[1]
                elif len(fit_data) == 3:
                    yf2 = fit_data[2]
                elif len(fit_data) == 4:
                    yf3 = fit_data[3]
                elif len(fit_data) == 5:
                    yf4 = fit_data[4]                    
                else:
                    print "More than 5 parameters are there in the fit here : ", datafile
                    
            ## Save data for plotting     
            save_data(f, 'MSR', x, y, yfit)
            save_fit_vars(f, 'ROT_DIFF/MSR', yf1)
            
        else:
            print "Analysis file could not be loaded for : ", datafile




##########################################################################            

if __name__ == '__main__':
    main()

