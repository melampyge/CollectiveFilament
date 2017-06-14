# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 08:39:40 2015

@author: jaeger
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import linalg
from scipy import optimize


"""compute the force caused by the harmonic potential of the bond and the active force along the bond"""
def f_bond_active(x1,y1,x2,y2):
    # geometric features
    delx = x1 - x2
    dely = y1 - y2
    rsq = delx**2 + dely**2
    r = math.sqrt(rsq)
    # conservative force
    dr = r - r0
    fbond = -2.0*kspring*dr
    fbondx=delx/r*fbond
    fbondy=dely/r*fbond
    # propulsive force
    fpx = delx*fp
    fpy = dely*fp
    # resulting forces
    fx1 = fbondx + 0.5*fpx
    fy1 = fbondy + 0.5*fpy
    fx2 = -fbondx + 0.5*fpx
    fy2 = -fbondy + 0.5*fpy
    return fx1, fy1, fx2, fy2
    
""" compute the force caused by the angular interactions"""
def f_angle(x1,y1,x2,y2,x3,y3):
    # geometric feautres
    delx1 = x1 - x2
    dely1 = y1 - y2
    rsq1 = delx1*delx1 + dely1*dely1 
    r1 = math.sqrt(rsq1)
    delx2 = x3 - x2
    dely2 = y3 - y2
    rsq2 = delx2*delx2 + dely2*dely2 
    r2 = math.sqrt(rsq2)
    # calculate factors for forces on each of the atoms
    cos = (delx1*delx2 + dely1*dely2)/(r1*r2)
    if cos > 1.0:
         cos = 1.0;
    if cos < -1.0:
         cos = -1.0;
    a11 = kappa*cos / rsq1
    a12 = -kappa / (r1*r2)
    a22 = kappa*cos / rsq2
    # resulting forces
    fx1 = a11*delx1 + a12*delx2
    fy1 = a11*dely1 + a12*dely2
    fx3 = a22*delx2 + a12*delx1
    fy3 = a22*dely2 + a12*dely1
    fx2 = -fx1-fx3
    fy2 = -fy1-fy3
    return fx1, fy1, fx2, fy2, fx3, fy3
    
    
def compute_forces(x,v):
    forces=np.zeros((4*n))
    #force due to rod
    forces[2*n]-=rod*gamma/m*v
    #forces due to bonds
    count=2*n
    for i in range(n-1):
        fx1, fy1, fx2, fy2=f_bond_active(x[count+i*2],x[count+i*2+1],x[count+i*2+2],x[count+i*2+3])
        forces[count+i*2]+=fx1/m
        forces[count+i*2+1]+=fy1/m
        forces[count+i*2+2]+=fx2/m
        forces[count+i*2+3]+=fy2/m
    #forces due to angle interactions
    for i in range(n-2):
        fx1, fy1, fx2, fy2, fx3, fy3=f_angle(x[count+i*2],x[count+i*2+1],x[count+i*2+2],x[count+i*2+3],x[count+i*2+4],x[count+i*2+5])
        forces[count+i*2]+=fx1/m
        forces[count+i*2+1]+=fy1/m
        forces[count+i*2+2]+=fx2/m
        forces[count+i*2+3]+=fy2/m
        forces[count+i*2+4]+=fx3/m
        forces[count+i*2+5]+=fy3/m
    #forces due to drag
    for i in range(n):
        forces[count+i*2]-=gamma/m*(x[i*2]+v)
        forces[count+i*2+1]-=gamma/m*(x[i*2+1])
    return forces

'''   
def cpompute_jacobian():
    
    
    
def eigenvalues(jac):
'''

r0 = 1.0
kspring = 2000
fp = 0.55 #force per unit length
n=10
rod=4
kappa=200
m=1.0
gamma=2.0

x1=0.;x2=-0.999661;x3=-1.99937;x4=-2.99912;x5=-3.99891;x6=-4.99875;x7=-5.99863;x8=-6.99856;x9=-7.99854;x10=-8.99855;v=0.160688
x=np.zeros((4*n))
x[2*n]=x1
x[2*n+2]=x2
x[2*n+4]=x3
x[2*n+6]=x4
x[2*n+8]=x5
x[2*n+10]=x6
x[2*n+12]=x7
x[2*n+14]=x8
x[2*n+16]=x9
x[2*n+18]=x10
f=compute_forces(x,v)
print map(np.round,f)

   
'''  
r0 = 1.0
kspring = 2000
fp = 0.5 #force per unit length
n=100
rod=10
kappa=200
m=1.0
gamma=2.0
    
x1=0.;x2=-0.998826;x3=-1.99766;x4=-2.99651;x5=-3.99538;x6=-4.99425;x7=-5.99314;x8=-6.99204
x9=-7.99095;x10=-8.98988;x11=-9.98882;x12=-10.9878;x13=-11.9867;x14=-12.9857;x15=-13.9847
x16=-14.9837;x17=-15.9827;x18=-16.9817;x19=-17.9808;x20=-18.9798;x21=-19.9789;x22=-20.9779
x23=-21.977;x24=-22.9761;x25=-23.9752;x26=-24.9744;x27=-25.9735;x28=-26.9727;x29=-27.9718
x30=-28.971;x31=-29.9702;x32=-30.9694;x33=-31.9686;x34=-32.9678;x35=-33.967;x36=-34.9663
x37=-35.9656;x38=-36.9648;x39=-37.9641;x40=-38.9634;x41=-39.9627;x42=-40.962;x43=-41.9614
x44=-42.9607;x45=-43.9601;x46=-44.9595;x47=-45.9589;x48=-46.9583;x49=-47.9577;x50=-48.9571
x51=-49.9565;x52=-50.956;x53=-51.9554;x54=-52.9549;x55=-53.9544;x56=-54.9539;x57=-55.9534
x58=-56.9529;x59=-57.9525;x60=-58.952;x61=-59.9516;x62=-60.9511;x63=-61.9507;x64=-62.9503
x65=-63.9499;x66=-64.9496;x67=-65.9492;x68=-66.9488;x69=-67.9485;x70=-68.9482;x71=-69.9479
x72=-70.9476;x73=-71.9473;x74=-72.947;x75=-73.9467;x76=-74.9465;x77=-75.9462;x78=-76.946
x79=-77.9458;x80=-78.9456;x81=-79.9454;x82=-80.9452;x83=-81.9451;x84=-82.9449;x85=-83.9448
x86=-84.9447;x87=-85.9445;x88=-86.9444;x89=-87.9444;x90=-88.9443;x91=-89.9442;x92=-90.9442
x93=-91.9441;x94=-92.9441;x95=-93.9441;x96=-94.9441;x97=-95.9441;x98=-96.9441;x99=-97.9442
x100=-98.9442;v=0.224873
    
x=np.zeros((4*n))
x[2*n+0]=x1
x[2*n+2]=x2
x[2*n+4]=x3
x[2*n+6]=x4
x[2*n+8]=x5
x[2*n+10]=x6
x[2*n+12]=x7
x[2*n+14]=x8
x[2*n+16]=x9
x[2*n+18]=x10
x[2*n+20]=x11
x[2*n+22]=x12
x[2*n+24]=x13
x[2*n+26]=x14
x[2*n+28]=x15
x[2*n+30]=x16
x[2*n+32]=x17
x[2*n+34]=x18
x[2*n+36]=x19
x[2*n+38]=x20
x[2*n+40]=x21
x[2*n+42]=x22
x[2*n+44]=x23
x[2*n+46]=x24
x[2*n+48]=x25
x[2*n+50]=x26
x[2*n+52]=x27
x[2*n+54]=x28
x[2*n+56]=x29
x[2*n+58]=x30
x[2*n+60]=x31
x[2*n+62]=x32
x[2*n+64]=x33
x[2*n+66]=x34
x[2*n+68]=x35
x[2*n+70]=x36
x[2*n+72]=x37
x[2*n+74]=x38
x[2*n+76]=x39
x[2*n+78]=x40
x[2*n+80]=x41
x[2*n+82]=x42
x[2*n+84]=x43
x[2*n+86]=x44
x[2*n+88]=x45
x[2*n+90]=x46
x[2*n+92]=x47
x[2*n+94]=x48
x[2*n+96]=x49
x[2*n+98]=x50
x[2*n+100]=x51
x[2*n+102]=x52
x[2*n+104]=x53
x[2*n+106]=x54
x[2*n+108]=x55
x[2*n+110]=x56
x[2*n+112]=x57
x[2*n+114]=x58
x[2*n+116]=x59
x[2*n+118]=x60
x[2*n+120]=x61
x[2*n+122]=x62
x[2*n+124]=x63
x[2*n+126]=x64
x[2*n+128]=x65
x[2*n+130]=x66
x[2*n+132]=x67
x[2*n+134]=x68
x[2*n+136]=x69
x[2*n+138]=x70
x[2*n+140]=x71
x[2*n+142]=x72
x[2*n+144]=x73
x[2*n+146]=x74
x[2*n+148]=x75
x[2*n+150]=x76
x[2*n+152]=x77
x[2*n+154]=x78
x[2*n+156]=x79
x[2*n+158]=x80
x[2*n+160]=x81
x[2*n+162]=x82
x[2*n+164]=x83
x[2*n+166]=x84
x[2*n+168]=x85
x[2*n+170]=x86
x[2*n+172]=x87
x[2*n+174]=x88
x[2*n+176]=x89
x[2*n+178]=x90
x[2*n+180]=x91
x[2*n+182]=x92
x[2*n+184]=x93
x[2*n+186]=x94
x[2*n+188]=x95
x[2*n+190]=x96
x[2*n+192]=x97
x[2*n+194]=x98
x[2*n+196]=x99
x[2*n+198]=x100
f=compute_forces(x,v)
print map(np.round,f)'''
    
    
    

    