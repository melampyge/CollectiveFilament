#!/usr/local/bin/python2.7

import matplotlib as mpl
mpl.use('Agg')
import sys
import numpy as np
import math
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

def read_coilicity():
    """ read in the coilcity"""
    t = []
    c = []
    ifile = open('coilicity.data')
    ifile.readline()
    ifile.readline()
    for line in ifile:
        line = line.split()
        try:
            t.append(float(line[0]))
            c.append(float(line[2]))
        except:
            pass
    ifile.close()
    t = np.array(t)
    c = np.array(c)
    # transform time and coility units
    ttrans = gamma*L**3/6./kT
    t *= dt/ttrans
    c *= L/2/np.pi
    return t,c

#####################################################################

def read_cacf():
    """ read in the coilicity autocorrelation function"""
    tacf = []
    cacf = []
    ifile = open('coil2_acf.data', 'r')
    ifile.readline()
    for line in ifile:
        line = line.split()
        tacf.append(float(line[0]))
        cacf.append(float(line[1]))
    ifile.close()
    tacf = np.array(tacf)
    cacf = np.array(cacf)
    # transform time units
    ttrans = gamma*L**3/6./kT
    tacf *= dt/ttrans
    return tacf,cacf

#####################################################################

def compute_moments(c):
    """ compute the coil moments"""
    n = len(c)
    cav = np.average(c)
    cav_std = np.std(c)/np.sqrt(n)
    csq = np.average(c**2)
    csq_std = np.std(c**2)/np.sqrt(n)
    curt = stats.kurtosis(c, fisher = False)
    # compute mirrored statistics
    cm = -np.copy(c)
    cboth = np.append(c,cm)
    curt2 = stats.kurtosis(cboth, fisher = False)
    return cav, cav_std, csq, csq_std,curt,curt2

#####################################################################

def compute_thalf(tacf,cacf):
    """ check where the autocorrelation function drops below 0.5"""
    n = len(tacf)
    thalf = -1
    for i in range(n):
        if cacf[i] < 0.5:
            thalf = tacf[i]
            break
    plt.plot(tacf, cacf)
    plt.savefig('coilicity_acf.png')
    plt.close()
    return thalf

#####################################################################

def main():
    """ main function"""
    # read in the coilicity
    t,c = read_coilicity()
    # read in the time autocorrelation function
    tacf, cacf = read_cacf()
    # compute the moments and standard deviations
    cav, cav_std, csq, csq_std,curt,curt2 = compute_moments(c)
    # compute the moments for only the second part of the array
    n = len(c)
    cavh, cav_stdh, csqh, csq_stdh, curth, curt2h = compute_moments(c[n/2:])
    # compute the time where the acf drops below 0.5
    thalf = compute_thalf(tacf,cacf)
    # write results to file
    ofile = open('coil_phase.data', 'w')
    ofile.write('Information required to identify coil phase\n\n')
    ofile.write('cav\tcav_std\tcsq\tcsq_std\tthalf\tcurt\tcurt2\n')
    ofile.write(str(cav) + '\t' + str(cav_std) + '\t' + str(csq) + '\t' + str(csq_std) + '\t' + str(thalf) + '\t' + str(curt) + '\t' + str(curt2) + '\n')
    ofile.close()
    ofile = open('coil_phaseh.data', 'w')
    ofile.write('Information required to identify coil phase\n\n')
    ofile.write('cav\tcav_std\tcsq\tcsq_std\tthalf\tcurt\tcurt2\n')
    ofile.write(str(cavh) + '\t' + str(cav_stdh) + '\t' + str(csqh) + '\t' + str(csq_stdh) + '\t' + str(thalf) + '\t' + str(curth) + '\t' + str(curt2h) + '\n')
    ofile.close()

    return
        
#####################################################################

if __name__ == '__main__':
    main()
