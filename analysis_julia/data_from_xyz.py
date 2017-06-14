# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 10:53:31 2015

@author: jaeger
"""

import numpy as np

#Read x and y positions from Lammps xyz-output
#A 2D array is returned with the data given in the following format: values[time][atom]
    #name= xyz-file
    #time= cut-off time for very long simulations. Set it equal 0 for the whole data to be used

def read_coordinates(name,time):
    n=extract_n(name)
    #for very long simulations this step should not be done and the number of time steps put in manually
    if time==0:
        f= open(name,'r')
        num_lines = sum(1 for line in f)
        sets=int(float(num_lines)/float(n+2))
        f.close()
    else:
        sets=time
    
    x_values=np.zeros((sets,n))
    y_values=np.zeros((sets,n))
    
    f= open(name,'r')
    step=-1
    atom=-1
    for i in xrange((n+2)*sets):
        st=f.readline()
        if st.split()[0] != str(n):    
            if st[0]=='A':
                step+=1
                atom=-1
            else:   
                atom+=1
                line=st.split()
                x_values[step][atom]=float(line[1])
                y_values[step][atom]=float(line[2])
    f.close()
    return [x_values,y_values,sets,n]
 
 
 
 
#Extract number of atoms from xyz-output 
 
def extract_n(name):
    f= open(name,'r')
    l=int(f.readline())
    f.close()
    return int(l)