#!/usr/local/bin/python2.7

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


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
ttrans = gamma*L**3/6./kT

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
    t *= dt/ttrans
    MSD /= L**2
    return mes_ok, t, MSD

#####################################################################

def get_fit(t):
    """ determine the values for dt, v0, and tr from theory"""
    # read in dt
    ifile = open('../../Pe_0.0/OUTPUT/msd_analyse_results_dt.data')
    ifile.readline()
    ifile.readline()
    line = ifile.readline()
    line = line.split()
    dt = float(line[-1])
    ifile.close()
    # read in tr
    ifile = open('msr_analyse_results.data')
    ifile.readline()
    ifile.readline()
    line = ifile.readline()
    line = line.split()
    tr = 1./float(line[-1])
    ifile.close()
    # read in re2e
    ifile = open('re2e_av.data')
    ifile.readline()
    ifile.readline()
    line = ifile.readline()
    line = line.split()
    re2e = float(line[-1])
    ifile.close()
    # read in fp
    ifile = open('../input.lammps')
    while True:
        line = ifile.readline()
        line = line.split()
        if len(line) >= 2:
            if line[1] == 'fp':
                fp = float(line[3])
                break
    ifile.close()
    # predict v0
    v0 = re2e*fp/gamma
    v0 *= ttrans/L
    # compute MSD
    Mfit = MSD_func(t,dt,tr,v0)
    return dt, tr, v0, Mfit

#####################################################################

def MSD_func(t,dt,tr,v0):
    """ computes the MSD according to theoretical prediction"""
    d = 2
    msd = 2*d*dt*t + 2*tr**2*v0**2*(t/tr + np.exp(-t/tr) - 1)
    return msd

#####################################################################

def MSD_func_log(t,dt,tr,v0):
    """ computes the log of the MSD according to theoretical prediction"""
    d = 2
    msd = 2*d*dt*t + 2*tr**2*v0**2*(t/tr + np.exp(-t/tr) - 1)
    msd = np.log(msd)
    return msd

#####################################################################

def fit_curves(t, MSD, dt0, tr0, v00):
    """ fit measured MSDs to theoretical prediction"""
    # starting array
    p0 = [dt0,tr0,v00]
    # select timespan for which to perform the fit
    tmax = t[-1]
    tfit = []
    Mfit = []
    for i in range(1,len(t)):
        if t[i] < tmax/20:
            tfit.append(t[i])
            Mfit.append(MSD[i])
    tfit = np.array(tfit)
    Mfit = np.array(Mfit)
    # perform the fit
    Mlog = np.log(Mfit)
    popt,pcov = optimize.curve_fit(MSD_func_log,tfit,Mlog, p0 = p0)
    dt = popt[0]
    tr = popt[1]
    v0 = popt[2]

    Mfit = MSD_func(t,dt,tr,v0)
    return dt,tr,v0,Mfit

#####################################################################

def gen_plot(t, MSD, MSD_fit, MSD_fit2, mes_ok):
    """ generate a plot of the data"""
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(t,MSD_fit)
    ax.plot(t,MSD_fit2, color = 'k')
    ax.plot(t,MSD, color = 'r', label = 'mes_ok = ' + str(mes_ok))
    ax.legend(loc = 'upper left')
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.show()
    plt.close()
    return

#####################################################################

def  write_results(dt, tr, v0, dt2, tr2, v02, mes_ok, t, MSD, MSD_fit, MSD_fit2):
     ofile = open('msd_analyse_results.data', 'w')
     ofile.write('Results from analyzing the MSD\n\n')
     ofile.write('dt = ' + str(dt) + '\n')
     ofile.write('tr = ' + str(tr) + '\n')
     ofile.write('v0 = ' + str(v0) + '\n')
     ofile.write('dt2 = ' + str(dt2) + '\n')
     ofile.write('tr2 = ' + str(tr2) + '\n')
     ofile.write('v02 = ' + str(v02) + '\n')
     ofile.write('mes_ok = ' + str(mes_ok) + '\n')
     ofile.write('\n')
     ofile.write('t\tMSD\tMSD_fit\tMSD_fit2\n')
     for i in range(len(t)):
         ofile.write(str(t[i]) + '\t' + str(MSD[i]) + '\t' + str(MSD_fit[i]) + str(MSD_fit2[i]) + '\n')
     ofile.close()
     return
#####################################################################

def main():
    """ main function"""
    # read in the MSD data
    mes_ok, t, MSD = read_data()
    # read in theoretical predictions
    dt, tr, v0, MSD_fit = get_fit(t)
    # fit curves
    if mes_ok == 1:
        dt2, tr2, v02, MSD_fit2 = fit_curves(t,MSD,dt,tr,v0)
    else:
        dt2 = dt
        tr2 = tr
        v02 = v0
        MSD_fit2 = MSD_fit
    # print fitted values
    print 'dt, tr, v0, dt2, tr2, v02', dt, tr, v0, dt2, tr2, v02
    # create plots with fitted and measured curves
    gen_plot(t, MSD, MSD_fit, MSD_fit2, mes_ok)
    # write results to file
    write_results(dt, tr, v0, dt2, tr2, v02, mes_ok, t, MSD, MSD_fit, MSD_fit2)
    return

#####################################################################

if __name__ == '__main__':
    main()
