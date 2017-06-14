#Make a data input file for lammps

n=105
head=5
f = open('input_5_105.data', 'w')

#print 'scaled time is:',int((n/100.0)**3*50000000)

f.write("Single Polymer Chain with ")
f.write(str(n))
f.write(" beads \n \n")

f.write(str(n))
f.write(' atoms \n')
f.write(str(n-1))
f.write(' bonds \n')
f.write(str(n-2))
f.write(' angles \n \n')

f.write(str('1 atom types \n2 bond types \n1 angle types \n'))

f.write(str('\n0.0 '))
f.write(str(float(n+1)))
f.write(str(' xlo xhi \n0.0 10.0 ylo yhi \n-2.5 2.5 zlo zhi \n'))

f.write(str('\nMasses \n \n1 1.0 \n'))

f.write(str('\nAtoms \n'))

for i in range(n):
    type=1
    #if i+1>head:
    #    type=1
    f.write(str('\n'))
    f.write(str(i+1))
    f.write(str(' '))
    f.write(str(1))
    f.write(str(' '))
    f.write(str(type))
    f.write(str(' '))    
    f.write(str(float(i)))
    f.write(str(' '))
    f.write(str(0.0))
    f.write(str(' '))
    f.write(str(0.0))


f.write(str('\n \nBonds \n'))

for i in range(n-1):
    type=2
    if i+1>head:
        type=1
    f.write(str('\n'))
    f.write(str(i+1))
    f.write(str(' '))
    f.write(str(type))
    f.write(str(' '))
    f.write(str(i+1))
    f.write(str(' '))    
    f.write(str(i+2))

f.write(str('\n \nAngles \n'))

for i in range(n-2):
    type=1
    #if i+1>head:
    #    type=1
    f.write(str('\n'))
    f.write(str(i+1))
    f.write(str(' '))
    f.write(str(type))
    f.write(str(' '))
    f.write(str(i+1))
    f.write(str(' '))    
    f.write(str(i+2))
    f.write(str(' '))    
    f.write(str(i+3))

       
f.close()
