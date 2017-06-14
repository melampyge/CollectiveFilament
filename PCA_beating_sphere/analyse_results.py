#!usr/bin/python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import sys
import math
import os
from scipy import linalg
from scipy import stats
from scipy import optimize

try:
    ampfile = sys.argv[1]
    eigenfile = sys.argv[2]
    corrfile = sys.argv[3]
    curvfile = sys.argv[4]
except:
    print 'Usage: ' + sys.argv[0] + '        amplitude evolution file           eigenvalue file        correlaton file     curvature file'
    exit()


##################################################################

def read_amplitudes(ampfile):
    """ read in the amplitudes"""
    t = []
    amp = []
    ifile = open(ampfile)
    for i in range(3):
        ifile.readline()
    for line in ifile:
        line = line.split()
        t.append(float(line[0]))
        amp.append([])
        for i in range(1, len(line)):
            amp[-1].append(float(line[i]))
    ifile.close()

    t = np.array(t, dtype = np.int32)
    amp = np.array(amp, dtype = np.float64)

    return t, amp

##################################################################

def read_correlation(corrfile):
    """ read in the correlation functions"""
    tc = []
    ac1 = []
    ac1_std = []
    ac2 = []
    ac2_std = []
    cc12 = []
    cc12_std = []
    ifile = open(corrfile)
    for i in range(3):
        ifile.readline()
    for line in ifile:
        line = line.split()
        tc.append(int(line[0]))
        ac1.append(float(line[1]))
        ac1_std.append(float(line[2]))
        ac2.append(float(line[3]))
        ac2_std.append(float(line[4]))
        cc12.append(float(line[5]))
        cc12_std.append(float(line[6]))
    ifile.close()
    tc = np.array(tc, dtype = np.int32)
    ac1 = np.array(ac1, dtype = np.float64)
    ac1_std = np.array(ac1_std, dtype = np.float64)
    ac2 = np.array(ac2, dtype = np.float64)
    ac2_std = np.array(ac2_std, dtype = np.float64)
    cc12 = np.array(cc12, dtype = np.float64)
    cc12_std = np.array(cc12_std, dtype = np.float64)
    return tc, ac1, ac1_std, ac2, ac2_std, cc12, cc12_std

##################################################################

def read_eigen(eigenfile):
    """ read in the eigenvalues and eigenvectors"""
    ifile = open(eigenfile)
    ifile.readline()
    ifile.readline()
    line = ifile.readline()
    line = line.split()
    n = int(line[2])
    eigval = np.zeros((n), dtype = np.float64)
    eigvec = np.zeros((n,n), dtype = np.float64)
    # read in eigenvalues
    line = ifile.readline()
    line = line.split()
    for i in range(n):
        eigval[i] = float(line[i+1])
    # read in eigenvectors
    for i in range(n):
        line = ifile.readline()
        line = line.split()
        for j in range(n):
            eigvec[i,j] = float(line[j+1])
    # close file
    ifile.close()
    return eigval, eigvec

##################################################################

def read_curvature(curvfile):
    """ read in the average curvature"""
    ifile = open(curvfile)
    for i in range(3):
        ifile.readline()
    curv_av = []
    curv_std = []
    for line in ifile:
        line = line.split()
        curv_av.append(float(line[1]))
        curv_std.append(float(line[2]))
    ifile.close()
    curv_av = np.array(curv_av, dtype = np.float64)
    curv_std = np.array(curv_std, dtype = np.float64)
    return curv_av, curv_std

##################################################################

def cart2pol(x,y):
    """ transform cartesian coordinates to polar coordinates"""
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y,x)
    return rho,  phi

##################################################################

def kurtosis2d(x,y):
    """ compute the 2D kurtosis of data x and y"""
    # combine to a single data set
    data = np.array([x,y])
    # average
    av = np.average(data, axis = 1)
    # compute the covariance matrix
    S = np.cov(data)
    # compute the inverse of S
    Sm = linalg.inv(S)


    # compute the kurtosis
    n = len(x)
    b2 = 0
    for i in range(n):
        d = data[:,i] - av
        b2 += np.dot(d, np.dot(Sm, d))**2
    b2 /= n
  
    return b2

##################################################################

def correlation_with_errorbars(x,y):
    """ compute the correlation of two signals including the std of the
        correlation using the blocking method"""
    xav = np.average(x)
    yav = np.average(y)
    n = len(x)
    blockdata = np.zeros((n))
    block_std = []
    block_uncert = []
    for i in range(n):
        blockdata[i] = (x[i] - xav)*(y[i] - yav)
    # compute the average of the blockdata
    correlation = np.average(blockdata)
    # compute the std using the blocking method
    while n >= 2:
        # compute std for current block
        c0 = 0
        for i in range(n):
            c0 += (blockdata[i] - correlation)**2
        c0 /= (n-1)
        # compute the standard deviation and fill it to blocking array
        c0 = np.sqrt(c0/(n-1))
        block_std.append(c0)
        block_uncert.append(c0/np.sqrt(2*(n-1)))
        # perform block operation
        for i in range(n/2):
            blockdata[i] = 0.5*(blockdata[2*i] + blockdata[2*i+1])
        n = n/2;

    # set initial value for the std
    std = 0.0
    ## find the plateau value and fill it to stdl, set to -1 if no plateau is found
    nignore = 3;
    n = len(block_std) - nignore
    smin_old = block_std[0] - 2*block_uncert[0]
    smax_old = block_std[0] + 2*block_uncert[0]
    for i in range(1,n): 
        smin = block_std[i] - 2*block_uncert[i]
        smax = block_std[i] + 2*block_uncert[i]
        if (smax_old > smin):
	        std = 0.5*(block_std[i] + block_std[i-1])
	        break
        smin_old = smin
        smax_old = smax
        

    return correlation, std

##################################################################

def reconstruct_shape(curv_av):
    """ reconstruct the principal modes from the average curvature"""
    n = len(curv_av)
    # allocate arrays
    xs = np.zeros((n+2), dtype = np.float64)
    ys = np.zeros((n+2), dtype = np.float64)
    # set starting value
    xs[1] = 1
    # loop over all curvatures
    for i in range(n):
        phi = curv_av[i]
        cos = np.cos(phi)  # tangent displacement
        sin = np.sin(phi)  # normal displacement
        dx = xs[i+1] - xs[i]
        dy = ys[i+1] - ys[i]
        dx2 = cos*dx - sin*dy
        dy2 = sin*dx + cos*dy
        xs[i+2] = xs[i+1] + dx2
        ys[i+2] = ys[i+1] + dy2
    return xs, ys

##################################################################

def reconstruct_principal_modes(eigvec, eigval):
    """ reconstruct the principal modes from the eigenvalues and amplitudes"""
    angle=[]
    N=eigvec.shape[0]
    xs=np.zeros(N)
    ys=np.zeros(N)
    angle=eigval*np.cumsum(eigvec)
    C=np.cos(angle)
    S=np.sin(angle)
    xs=np.cumsum(C)
    ys=np.cumsum(S)
    return xs, ys

##################################################################

def cross_correlate(time, x1, x2, dt, limit):
    """ compute the autocorrelation function of x"""
    ### compute required values
    nsteps = len(x1)
    limit = int(limit*nsteps)
    ncc = int(limit/dt)
    ### allocate output array to store the results
    t_cc = np.zeros((ncc), dtype = np.int32)
    cc = np.zeros((ncc), dtype = np.float64)
    linval = np.zeros((ncc), dtype = np.int32)
    std = np.zeros((ncc), dtype = np.float64)
    blockdata = np.zeros((nsteps), dtype = np.float64)
    # size of the array for the blocking method
    nblock = int(math.log(nsteps, 2)) + 2
    block_std = np.zeros((nblock), dtype = np.float64)
    block_uncert = np.zeros((nblock), dtype = np.float64)
    ### loop over all molecules, submit MSD computation, read and average
    ###   the results
    performance_tools = performance_toolsWrapper.Performance_tools()
    # submit MSR computation
    performance_tools.compute_crosscorrelation_with_errorbars(t_cc, time, x1, x2, cc, std, blockdata, block_std, block_uncert, linval, nsteps, ncc, limit, dt)
    ### return time and msr
    return t_cc, cc, std

##################################################################

def zero_passing_time(t,hyac):
    """ find where the autocorrelation function passes through zero
        for the first time"""
    n = len(hyac)
    idx = 0
    for i in range(n-1):
        if hyac[i]*hyac[i+1] < 0:
            idx = i
            break
    # find zero passing time from linear interpolation
    tm = t[idx]
    tp = t[idx+1]
    hm = hyac[idx]
    hp = hyac[idx+1]
    
    t0 = tm + (tp - tm) / (hp - hm) * (0.0 - hm)
    return t0

##################################################################

def find_extrema(t,hyac,std):
    """ find where the autocorrelation function passes through zero
        for the first time"""
    n = len(hyac)
    te = []
    he = []
    stde = []
    for i in range(1, n-1):
        if (hyac[i-1] - hyac[i]) * ( hyac[i] - hyac[i+1])  < 0:
            te.append(t[i])
            he.append(hyac[i])
            stde.append(std[i])
    te = np.array(te, dtype = np.int32)
    he = np.array(he, dtype = np.float64)
    stde = np.array(stde, dtype =  np.float64)
    # loop over all extrema and check whether they are valid
    nextr = len(te)
    valid = 0
    for i in range(nextr):
        hi = he[i]
        stdi = stde[i]
        print 'hi, stdi, i', hi, stdi, i
        if valid % 2 == 0:    # next item should be a minimum
            if hi + 2*stdi > 0:
                break
        else:                 # next item should be a maximum
            if hi - 2*stdi < 0:
                break
        valid = valid + 1
    return te, he, stde, valid

##################################################################

def exponential_decay_a0_fixed(a0):
    def func(t,tau_decay, tau_beat):
        y = a0*np.cos(2*np.pi*t/tau_beat)*np.exp(-t/tau_decay)
        return y
    return func

##################################################################

def exponential_decay(t,tau_decay, tau_beat, a0):
    y = a0*np.cos(2*np.pi*t/tau_beat)*np.exp(-t/tau_decay)
    return y

##################################################################

def compute_decay_and_frequency(tc, ac, tce, ace, v_ace):
    """ compute the decay time and frequency from the autocorrelation functions"""
    # check whether there is at least one significant minimum
    if v_ace < 1 :
        return 1.0, 1.0, 0.0
    # compute first guesses for all parameters
    a0 = ac[0]
    a1 = math.fabs(ace[0])
    if a1 > a0:
        a1 = 0.9999999*a0
    t_decay_0 = tce[0]/math.log(a0/a1)
    t_beat_0 = 2*tce[0]
    # fit damped cosine function to the results
    try:
        popt, pcov = optimize.curve_fit(exponential_decay_a0_fixed(a0), tc, ac, p0 = [t_decay_0, t_beat_0], maxfev = 2000)
        t_decay = popt[0]
        t_beat = popt[1]
    except:
        t_decay = t_decay_0
        t_beat = t_beat_0
    
    """   alternative fitting in which a0 is also optimized
    popt, pcov =optimize.curve_fit(exponential_decay, tc, ac, p0 = [t_decay_0, t_beat_0, a0])
    t_decay = popt[0]
    t_beat = popt[1]
    a0 = popt[2]
    """
    return t_decay, t_beat, a0

##################################################################

def multivariate_histogram(x,y,fname):
    """ compute multivariate histograms"""
    ### generate a histogram of amp1 and amp1
    fig = plt.figure()

    nullfmt = NullFormatter()
    
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02

    rect_hist2d = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    ax = plt.axes(rect_hist2d)
    ax1 = plt.axes(rect_histx)
    ax2 = plt.axes(rect_histy)

    ax1.xaxis.set_major_formatter(nullfmt)
    ax2.yaxis.set_major_formatter(nullfmt)

    nbins = 50
    ax.hist2d(x, y, bins = nbins, cmap = 'copper')
    ax1.hist(x, bins = nbins)
    ax2.hist(y, bins = nbins, orientation = 'horizontal')
    
    plt.savefig(fname)
    plt.close()

    return
    
##################################################################
    
def main():
    """ main function, called when the script is started"""
    ### read in the amplitutes
    print '   reading amplitudes'
    t, amp = read_amplitudes(ampfile)

    ### read in the correlations files
    tc, ac1, ac1_std, ac2, ac2_std, cc12, cc12_std = read_correlation(corrfile)
    
    ### read in the eigenvalues and functions
    eigval, eigvec = read_eigen(eigenfile)

    ### read in the average curvature
    curv_av, curv_std = read_curvature(curvfile)
        
    ### analysis of the amplitudes
    # average of of squared amplitudes
    amp_av = np.sqrt(np.average(amp**2, axis = 0))
    
    # compute the variance of the "circle" made by a1 and a2
    amp_var = np.var((amp[:,0]/amp_av[0])**2 + (amp[:,1]/amp_av[1])**2)

    # compute the multivariate curtosis of amp1 and amp2
    b2 = kurtosis2d(amp[:,0], amp[:,1])

    # compute the correlations of higher order amplitudes with a1 and a2
    namp = len(amp[0])
    c = np.zeros((namp), dtype = np.float64)
    c_std = np.zeros((namp), dtype = np.float64)
    for i in range(namp):
        c[i], c_std[i] = correlation_with_errorbars(amp[:,0]**2 + amp[:,1]**2, amp[:,i]**2)

    # transform data of a1 and a2 to polar coordinates
    rho, phi = cart2pol(amp[:,0]/amp_av[0], amp[:,1]/amp_av[1])

    ### reconstruct the shape of the eigenmodes from the eigenvalues
    x0, y0 = reconstruct_principal_modes(curv_av, 1)
    x1, y1 = reconstruct_principal_modes(eigvec[0], eigval[0])
    x2, y2 = reconstruct_principal_modes(eigvec[1], eigval[1])
    x3, y3 = reconstruct_principal_modes(eigvec[2], eigval[2])
    x4, y4 = reconstruct_principal_modes(eigvec[3], eigval[3])

    # zero passing time
    tc1_0 = zero_passing_time(tc,ac1)
    tc2_0 = zero_passing_time(tc,ac2)
    tc12_0 = zero_passing_time(tc[1:],cc12[1:]) # don't use first value

    # significant extrema
    tce1, ac1e, ac1e_std, ac1e_v = find_extrema(tc, ac1, ac1_std)
    tce2, ac2e, ac2e_std, ac2e_v = find_extrema(tc, ac2, ac2_std)
    tce12, cc12e, cc12e_std, cc12e_v = find_extrema(tc, cc12, cc12_std)

    # compute decay time and frequency from the autocorrelation functions
    try:
       ac1_tau, ac1_T, ac1_a0 = compute_decay_and_frequency(tc[tc <= tce1[ac1e_v]], ac1[tc <= tce1[ac1e_v]], tce1, ac1e, ac1e_v)
    except:
        try:
           ac1_tau, ac1_T, ac1_a0 = compute_decay_and_frequency(tc[tc <= tce1[ac1e_v-1]], ac1[tc <= tce1[ac1e_v-1]], tce1, ac1e, ac1e_v-1)
        except:
            ac1_tau = 1.0
            ac1_T = 1.0
            ac1_a0 = 0.0
    try:   
       ac2_tau, ac2_T, ac2_a0 = compute_decay_and_frequency(tc[tc <= tce2[ac2e_v]], ac2[tc <= tce2[ac2e_v]], tce2, ac2e, ac2e_v)
    except:
        try:
            ac2_tau, ac2_T, ac2_a0 = compute_decay_and_frequency(tc[tc <= tce2[ac2e_v-1]], ac2[tc <= tce2[ac2e_v-1]], tce2, ac2e, ac2e_v-1)
        except:
            ac2_tau = 1.0
            ac2_T = 1.0
            ac2_a0 = 0.0
    
    ### write down results to file
    ofile = open('results.data', 'w')
    ofile.write('Results from analysing the amplitudes\n\n')
    # average amplitudes
    ofile.write('amp_1_sq_av = ' + str(amp_av[0]) + '\n')
    ofile.write('amp_2_sq_av = ' + str(amp_av[1]) + '\n')
    ofile.write('amp_3_sq_av = ' + str(amp_av[2]) + '\n')
    # normalized variance of a1 and a2
    ofile.write('variance_a1_a2 = ' + str(amp_var) + '\n')
    # 2d-kurtosis of a1 and a2
    ofile.write('kurtosis_a1_a2 = ' + str(b2) + '\n')
    # a0, decay, and beating time
    ofile.write('amp1: a0, tau_decay, tau_beat = ' + str(ac1_a0) + ' ' + str(ac1_tau) + ' ' + str(ac1_T) + '\n')
    ofile.write('amp2: a0, tau_decay, tau_beat = ' + str(ac2_a0) + ' ' + str(ac2_tau) + ' ' + str(ac2_T) + '\n')
    ofile.close()

    
    ### generate plots
    # histogram of a1 and a2
    multivariate_histogram(amp[:,0], amp[:,1], 'hist_a1_a2.png')
    # histogram of a1**2 + a2**2 and a3
    multivariate_histogram((amp[:,0]**2 + amp[:,1]**2)*np.sign(amp[:,0]), amp[:,2], 'hist_a1a2_a3.png')
    # histogram of a1**2 + a2**2 and a4
    multivariate_histogram((amp[:,0]**2 + amp[:,1]**2)*np.sign(amp[:,0]), amp[:,3], 'hist_a1a2_a4.png')
    # histogram of phi
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.hist(phi, bins = 50)
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$p(\phi)$')
    plt.savefig('hist_phi.png')
    plt.close()
    # histogram of r
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.hist(rho, weights = 1./(2*np.pi*rho), bins = 100)
    #ax.hist(rho, weights = 2*np.pi*rho, bins = 100)
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$p_w(\rho)$')
    plt.savefig('hist_rho.png')
    plt.close()
    # amplitude correlations
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.errorbar(np.linspace(1,namp, namp), c, yerr = 2*c_std, ls = '', marker = 'o')
    ax.set_xlabel('amplitude')
    ax.set_ylabel('correlation with a1**2 + a2**2')
    plt.savefig('amplitude_correlations.png')
    plt.close()

    # cummulative sum of the eigenvalues
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(np.cumsum(eigval)/np.sum(eigval), ls = '', marker = 'o')
    plt.savefig('eigenvalues.png')
    plt.close()
    # reconstructed shape of the eigenmodes
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(x1, y1, label = 'mode 1')
    ax.plot(x2, y2, label = 'mode 2')
    ax.plot(x3, y3, label = 'mode 3')
    ax.plot(x4, y4, label = 'mode 4')
    ax.legend(loc = 'lower left')
    plt.savefig('eigenmodes.png')
    plt.close()

    # average shape of the molecule
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.plot(x0, y0)
    ax.set_aspect('equal')
    plt.savefig('average_shape.png')
    plt.close()

    # auto- and cross-correlation functions
    fig = plt.figure()
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(312)
    ax3 = plt.subplot(313)
    try:
        tcrit = tce1[0]*6
    except:
        tcrit = tc[-1]
    sel = (tc <= tcrit)
    ax1.errorbar(tc[sel], ac1[sel], yerr = 2*ac1_std[sel])
    ax1.plot(tc[sel], exponential_decay(tc[sel], ac1_tau, ac1_T, ac1_a0))
    ax2.errorbar(tc[sel], ac2[sel], yerr = 2*ac2_std[sel])
    ax2.plot(tc[sel], exponential_decay(tc[sel], ac2_tau, ac2_T, ac2_a0))
    ax3.errorbar(tc[sel], cc12[sel], yerr = 2*cc12_std[sel]) 
    plt.savefig('correlation_funcs.png')
    plt.close()
 
    return
    
    
##################################################################

if __name__ == '__main__':
    main()
    
