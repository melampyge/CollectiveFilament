
import numpy as np
import sys, math, os
import codecs

try:
    fname = sys.argv[1]
    ofname = sys.argv[2]
    hname = sys.argv[3]
    nsnap = int(sys.argv[4])
    last_T = int(sys.argv[5])
except:
    try:
        fname = sys.argv[1]
        ofname = sys.argv[2]
        hname = sys.argv[3]
        nsnap = 0
        last_T = 99950000
    except:
        print 'Usage: ' + sys.argv[0] + '      infilename         outfilename     header filename      nsnap (optional)     last timestep (optional)'
        exit()

##################################################################

def determine_backup():
    """ determine the backup configuration, find the last timestep read before"""
    
    last_tstep = -1
    cnt = -1
    if os.path.exists(hname):
        header_file = open(hname, 'r')
        for line in header_file:
            A = line.split()
            # if we are at timestep line, read it
            if cnt == 0:
                cnt = -1
                last_tstep = int(A[0])
                continue
            # find the timestep line
            if len(A) == 2:
                if A[1] == 'TIMESTEP':
                    cnt += 1
        header_file.close()
                                        
    return last_tstep


##################################################################

def compress_data(nsnap, first_timestep):
    """ compress the data into a two byte file"""
    # input file list including restart simulation runs
    files_to_try = []
    last_snap_in_file = {}
    # determine the number of atoms
    ifile = open(fname)
    ifile.readline()
    ifile.readline()
    ifile.readline()
    line = ifile.readline()
    line = line.split()
    natoms = int(line[0])
    ifile.close()
    # determine the number of snapshots
    if nsnap == 0:
        files_to_try = [fname]
        os.system('wc -l ' + fname + ' > tmp.txt')
        ifile = open('tmp.txt')
        line = ifile.readline()
        line = line.split()
        nlines = int(line[0])
        nsnap = nlines/(natoms+9)
        ifile.close()
        os.system('rm tmp.txt')
        last_snap_in_file[fname] = nsnap
    # if the total number of snapshots is given, so that the restart runs are being combined
    else:
        total_snap = 0
        for j in range(1,20):
            path_to_try = fname.split('.')[0][:3] + str(j) + '.dump'
            print path_to_try
            os.system('grep ' + str(last_T) + ' ' + path_to_try + ' | tail -1 > tmp2.txt')
            os.system('wc -l ' + path_to_try + ' > tmp.txt')
            ifile = open('tmp.txt')
            line = ifile.readline()
            line = line.split()
            nlines = int(line[0])
            snaps = nlines/(natoms+9)
            total_snap += snaps
            last_snap_in_file[path_to_try] = snaps
            ifile.close()
            os.system('rm tmp.txt')
            ifile = open('tmp2.txt')
            line = ifile.readline()
            ifile.close()
            files_to_try.append(path_to_try)
            os.system('rm tmp2.txt')            
            if len(line) != 0:
                break
    print files_to_try        
    # allocate an array to store the x and y coords and modified coords
    x = np.zeros((natoms))
    y = np.zeros((natoms))
    xi = np.zeros((natoms), dtype = int)
    yi = np.zeros((natoms), dtype = int)
    x1 = np.zeros((natoms), dtype = int)
    x2 = np.zeros((natoms), dtype = int)
    y1 = np.zeros((natoms), dtype = int)
    y2 = np.zeros((natoms), dtype = int)
    # open files, start to compress the data
    #hfile = open(hname,'a')
    #hfile.write('nsnapshots, natoms = ' + str(nsnap) + ' ' + str(natoms) + '\n')    
    ofile = codecs.open(ofname, 'a', 'UTF-8')
    snap_counter = 0
    last_timestep = -1      # last read timestep during reading of dump files in one run
    break_flag = False
    ##########
    ### COMBINE TIMESTEP MODULE WITH A BACKUP COMPRESSION MODULE
    #########
    # for all the dump files
    for fn in files_to_try:
        ifile = open(fn)
        hfile = open(hname,'a')        
        # run until the last timestep in each dump file
        for i in range(last_snap_in_file[fn]):
            # stop saving data if the flag is on
            if break_flag:
                print 'Run has ended with total number of snapshots / last timestep', snap_counter, last_T
                break
            # copy the header if the timestep is unique
            line1 = ifile.readline()
            line2 = ifile.readline()
            A = line2.split()
            timestep = int(A[0])     
            # if the current timestep is recorded before, don't save it again
            # last_timestep is the last timestep recorded in this run through the dump files
            # first_timestep is the first timestep of this run, because in the last runs other steps are added
            if timestep <= last_timestep or timestep <= first_timestep:
                for i in range(7+natoms):
                    ifile.readline()
                continue 
            # stop reading data once the targeted time step is exceeded
            snap_counter += 1
            if timestep == last_T:
                break_flag = True
            print 'current snapshot / all snapshots / timestep / file: ', snap_counter, nsnap, timestep, fn            
            # save the last time step in each file to compare it later on
            if i == last_snap_in_file[fn]-1:
                last_timestep = timestep
            print line1, line2
            hfile.write(line1)
            hfile.write(line2)
            for j in range(7):
                line = ifile.readline()
                hfile.write(line)
            # read in the coords
            for j in range(natoms):
                line = ifile.readline()
                line = line.split()
                aID = int(line[0]) - 1
                xs = float(line[2])
                ys = float(line[3])
                xs = xs - math.floor(xs)
                ys = ys - math.floor(ys)
                x[aID] = xs
                y[aID] = ys
            # transfrom the coords to integers and bytes to store
            for j in range(natoms):
                xi[j] = x[j]*256**2
                yi[j] = y[j]*256**2
                x1[j] = xi[j]/256
                y1[j] = yi[j]/256
                x2[j] = xi[j]%256
                y2[j] = yi[j]%256
            # write down the bytes as utf8 code to a file
            for j in range(natoms):
                u1 = unichr(x1[j])
                u2 = unichr(x2[j])
                u3 = unichr(y1[j])
                u4 = unichr(y2[j])
                ofile.write(u1)
                ofile.write(u2)
                ofile.write(u3)
                ofile.write(u4)
        ifile.close()
        hfile.close()
    ofile.close()
    os.system('sed -i \"s/nsnapshots, natoms = ' + str(nsnap) + ' ' + str(natoms) + '/nsnapshots, natoms = ' + str(snap_counter) + ' ' + str(natoms) + '/g\" out1.header')    
    return

##################################################################

def main():
    """ main function, called when script is started"""
    # compress the data
    first_T = determine_backup()
    compress_data(nsnap, first_T)
    return

##################################################################

if __name__ == '__main__':
    main()
