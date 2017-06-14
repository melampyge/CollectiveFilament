#!/usr/local/bin/python2.7


import numpy as np
import sys
import math
import os
import performance_toolsWrapper
from scipy import optimize

try:
    tanfilename = sys.argv[1]
    centerfilename = sys.argv[2]
except:
    print 'Usage: ' + sys.argv[0] + '       tangent vector file        center position file'
    exit()


##################################################################

def read_data(infilename):
    """ read in the time and tangent vectors"""
    t = []
    tx = []
    ty = []
    ifile = open(infilename)
    for i in range(3):
        ifile.readline()
    for line in ifile:
        line = line.split()
        t.append(int(line[0]))
        tx.append(float(line[1]))
        ty.append(float(line[2]))
    ifile.close()
    t = np.array(t, dtype = np.int32)
    tx = np.array(tx, dtype = np.float64)
    ty = np.array(ty, dtype = np.float64)
    return t, tx, ty

##################################################################

def calc_R(x, y, xc, yc):
    """ calculate the distance of each 2D points from the center (xc, yc)"""
    return np.sqrt((x-xc)**2 + (y-yc)**2)

##################################################################

def func(x,y):
    """ calculate the algebraic distance between the data points and the
        mean circle centered at c=(xc, yc)"""
    def func_2(center):
        xc = center[0]
        yc = center[1]
        Ri = calc_R(x,y,xc,yc)
        return Ri - Ri.mean()
    return func_2

##################################################################

def fit_arc(x,y):
    """ fit a circular arc to given x and y coordinates"""
    # compute the starting points for the optimization
    xc = np.average(x)
    yc = np.average(y)
    x0 = [xc, yc]

    xout, ier = optimize.leastsq(func(x,y), x0)
    
    xc_2, yc_2 = xout
    Ri_2 = calc_R(x,y,xc_2, yc_2)
    R_2 = Ri_2.mean()
    residu_2 = sum((Ri_2 - R_2)**2)
    return R_2

##################################################################

def compute_radius(x,y,tradius):
    """ compute the swimming radius using a circular fit"""
    n = len(x) - tradius
    radius = np.zeros((n), dtype = np.float64)
    for i in range(n):
        radius[i] = fit_arc(x[i:i+tradius], y[i:i+tradius])
    return radius

##################################################################

def main():
    """ main function, called when the script is started"""
    ### read parameters from input file
    print 'reading in the data'
    t, tx, ty = read_data(tanfilename)
    t, xc, yc = read_data(centerfilename)
    
    print 'computing the angular velocity omega'
    dtx = np.gradient(tx)
    dty = np.gradient(ty)
    omega = dty*tx - dtx*ty
    omega_av = np.average(np.fabs(omega))
    omega_std = np.std(np.fabs(omega)) / np.sqrt(len(omega))

    print 'computing the average radius'
    # use tfrac time of an entire period to commpute the radius
    tfrac = 0.25
    tradius = int(tfrac*2*np.pi/omega_av)
    radius = compute_radius(xc, yc, tradius)
    radius_av = np.average(radius)
    radius_std = np.std(radius) / np.sqrt(len(radius))
    
    print 'writing results to file'
    ofile = open('angular_velocity.data', 'w')
    ofile.write('Angular velocity of the rigid head group\n\n')
    ofile.write('average_omega = ' + str(omega_av) + ' +- ' + str(omega_std) + '\n')
    ofile.write('t\tomega\n')
    for i in range(len(t)):
        ofile.write(str(t[i]) + '\t' + str(omega[i]) + '\n')
    ofile.close()

    ofile = open('radius.data', 'w')
    ofile.write('Swimming radius of the rigid head group\n\n')
    ofile.write('average_radius = ' + str(radius_av) + ' +- ' + str(radius_std) + '\n')
    ofile.write('t\tradius\n')
    for i in range(len(radius)):
        ofile.write(str(t[i]) + '\t' + str(radius[i]) + '\n')
    ofile.close()

    

##################################################################

if __name__ == '__main__':
    main()
    
