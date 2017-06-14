#!/usr/local/bin/python2.7

## compute the velocity in x and y dimension over time
## by using finite differences

import sys, numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

try:
    fname = sys.argv[1]
    npoints = int(sys.argv[2])
    d = float(sys.argv[3])
except:
    print 'Usage: ' + sys.argv[0] + '      infilename        npoints        framewidth'
    exit()


########################################################

def read_data():
    """ read in the positions of the COM"""
    t = []
    x = []
    y = []
    z = []
    ifile = open(fname, 'r')
    ifile.readline()
    ifile.readline()
    for line in ifile:
        line = line.split()
	t.append(float(line[0]))
        x.append(float(line[1]))
        y.append(float(line[2]))
        z.append(float(line[3]))
    ifile.close()
    t = np.array(t)
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    return t,x,y,z

############################################################

def run(dummy):
    # fill the data
    global counter
    counter = counter + 1
    if counter == 0:
	for i in range(npoints):
            xdata[i] = x[0]
	    ydata[i] = y[0]
    rep = counter % npoints
    xdata[rep] = x[counter]
    ydata[rep] = y[counter]
    # adjust the box windows
    xdmin = np.min(xdata)
    xdmax = np.max(xdata)
    ydmin = np.min(ydata)
    ydmax = np.max(ydata)
    xfmin, xfmax = ax.get_xlim()
    yfmin, yfmax = ax.get_ylim()


    # adjust axis if necessary
    if xdmax > xfmax or xdmin < xfmin or ydmax > yfmax or ydmin < yfmin:
        xc = 0.5*xdmax + 0.5*xdmin
        yc = 0.5*ydmax + 0.5*ydmin
        xfmin = xc - d
        xfmax = xc + d
        yfmin = yc - d
	yfmax = yc + d
	ax.set_xlim([xfmin, xfmax])
	ax.set_ylim([yfmin, yfmax])
        #ax.figure.canvas.draw()
    time_text.set_text(time_template%(t[counter]))
    
    #if counter >= npoints:
    line.set_data(xdata,ydata)

       
    # return the line
    return line, time_text

############################################################

# read in the center of mass positions
print '  reading the data ...'
t,x,y,z = read_data()
# create a figure
fig = plt.figure(figsize = (8,8))
ax = plt.subplot()
line, = ax.plot([],[], ls = '', marker = 'o')
time_template = 'time = %.0f'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
time_text.set_text('')
#ax.set_ylim(-1,1)
#ax.set_xlim(-1,1)
ax.grid()
counter = -1
xdata = np.zeros((npoints))
ydata = np.zeros((npoints))
    
ani = animation.FuncAnimation(fig, run, blit=False, interval=25, repeat = False)
plt.show()
plt.close()

