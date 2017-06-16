
""" Helper functions for reading and writing data"""

##############################################################################

import misc_tools

##############################################################################

def read_sim_info(filename):
    """ load general simulation info"""
    
    fl = open(filename, 'r')

    for line in fl:
        A = line.split()
        if A[0] == "dt":                    # Time interval between MD steps
            dt = float(A[-1])
        elif A[0] == "ti":                  # Beginning time for data acquisition
            ti = float(A[-1])
        elif A[0] == "Lx":                  # Box size in x
            lx = float(A[-1])            
        elif A[0] == "Ly":                  # Box size in y
            ly = float(A[-1])
        elif A[0] == "totalStep":           # Total MD steps
            totalStep = float(A[-1])
        elif A[0] == "nsamp":               # Data sampling frequency
            nsamp = float(A[-1])
        elif A[0] == "nfil":                # Number of particles per polymer
            nbpf = float(A[-1])
        elif A[0] == "L":                   # Number of particles
            nbeads = float(A[-1])
        elif A[0] == "B":                   # Bond length between particles of a body
            bl = float(A[-1])
        elif A[0] == "kT":                  # Boltzmann constant*Temperature
            kT = float(A[-1])
        elif A[0] == "Fmc":                 # Self propulsion force constant
            fmc = float(A[-1])     
        elif A[0] == "Kbend":               # Bending constant
            kbend = float(A[-1])
            
    sim = misc_tools.Simulation(dt, ti, lx, ly, totalStep, nsamp, nbpf, \
                                nbeads, bl, kT, fmc, kbend)
    
    fl.close()
    
    return sim

##############################################################################
            
