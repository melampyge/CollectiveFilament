#!/usr/local/bin/python2.7

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


#####################################################################
### define / read in some global variables
gamma = 2.0  # viscosity
kT = 1.0     # thermal energy

ifile = open('scale_params.data')
line = ifile.readline()
line = line.split()
L = float(line[-1])    # polymer length
line = ifile.readline()
line = line.split()
dt = float(line[-1])     # simulation timestep

#####################################################################

def read_data():
    """ read in the data"""
    # basic variables
    t = [0]
    MSD = [0]
    # read in the results
    ifile = open('msd_logscale.data')
    ifile.readline()
    for line in ifile:
        line = line.split()
        ti = float(line[0])
        d = float(line[1])
        if ti != t[-1]:
            t.append(ti)
            MSD.append(d)
    ifile.close()
    # check for problems
    mes_ok = 1
    if t[-1] == 0:
        mes_ok = 0
    # transfrom to arrays
    t = np.array(t)
    MSD = np.array(MSD)
    # scale the results
    MSD /= L**2
    ttrans = gamma*L**3/6./kT
    t *= dt/ttrans
    return mes_ok, t, MSD

#####################################################################

def find_dr(t,MSD):
    """ compute the rotational diffusion coefficient"""
    # compute log versions of the data
    tlog = np.log10(t)
    Mlog = np.log10(MSD)
    kmin = MSD[0] / t[0]
    kmax = MSD[-1] / t[-1]
    n = 6
    k = np.logspace(np.log10(kmin), np.log10(kmax),n)
    s = []
    for i in range(n):
        s.append(np.log10(k[i]*t))
    # find the linear regime of the data
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(tlog,Mlog)
    for i in range(n):
        ax.plot(tlog,s[i],lw = 0.5, color = '0.5')

    ax.grid(which = 'major')
    ax.grid(which = 'minor')
    plt.ion()
    plt.show()
    while True:
        try:
            # read in the input and output values
            inpt = raw_input('lower and upper bounds:  ')
            inpt = inpt.split()
            tmin = float(inpt[0])
            tmax = float(inpt[1])
            # compute the slope of the selected regime
            res = stats.linregress(tlog[(tlog >= tmin) & (tlog <= tmax)],Mlog[(tlog >= tmin ) & ( tlog <= tmax) ])
            lslope = res[0]
            inpt = raw_input('slope = ' + str(lslope) + '. Accept? (0 / 1):  ')
            inpt = int(inpt)
            if inpt == 1:
                break
        except:
            pass
    plt.close()
    # perform linear fit
    res  = stats.linregress(t[(tlog >= tmin) & (tlog <= tmax)],MSD[(tlog >= tmin ) & ( tlog <= tmax) ])
    slope = res[0]
    dt = 0.25*slope
    print 'Dt = ' + str(dt)
    return dt, lslope

#####################################################################

def  write_results(t,MSD,dt,mes_ok,lslope):
     ofile = open('msd_analyse_results_dt.data', 'w')
     ofile.write('Results from analyzing the MSD\n\n')
     ofile.write('dt = ' + str(dt) + '\n')
     ofile.write('mes_ok = ' + str(mes_ok) + '\n')
     ofile.write('lslope = ' + str(lslope) + '\n')
     ofile.write('\n')
     ofile.write('t\tMSD\n')
     for i in range(len(t)):
         ofile.write(str(t[i]) + '\t' + str(MSD[i]) + '\n')
     ofile.close()
     return
#####################################################################

def main():
    """ main function"""
    # read in the MSD data
    mes_ok, t, MSD = read_data()
    # fit curves
    dt, lslope = find_dr(t[1:],MSD[1:])
    # write results to file
    write_results(t,MSD,dt,mes_ok,lslope)
    return

#####################################################################

if __name__ == '__main__':
    main()
