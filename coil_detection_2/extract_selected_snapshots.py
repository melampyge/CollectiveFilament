#!/usr/local/bin/python2.7

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import sys, math

extract = [0,25737,59517,169437,178016,200000]

ifile = open('run.xyz')
ofile = open('selected.xyz', 'w')

for i in range(200001):
    if i in extract:
        for line in range(102):
            line = ifile.readline()
            ofile.write(line)
    else:
        for line in range(102):
            ifile.readline()
ifile.close()
ofile.close()
