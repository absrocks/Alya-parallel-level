import sys
import math
import numpy as np


meshName=sys.argv[1]


coordfile       = '{}.coord'.format(meshName)
velofile        = 'VELOC.alya'
pressfile        = 'PRESS.alya'

fCoord          =open(coordfile,'r')
fVel            =open(velofile,'w')
fPress           =open(pressfile,'w')

pi = 3.14159265
L = 2*pi
V0 = 1.0

print('---| Start writing initial condition')
for line in fCoord:
    data=line.split()

    pid = int(data[0])
    dims = len(data)-1
    x   = float(data[1])
    y   = float(data[2])
    z   = float(data[3])

    vx =  V0*np.sin(x/(0.2*L))*np.cos(y/(0.2*L))*np.sin(z/(0.2*L))
    vy = -V0*np.cos(x/(0.2*L))*np.sin(y/(0.2*L))*np.sin(z/(0.2*L))
    vz = 0.0
    pr = 0.0

    fVel.write('{} {} {} {}\n'.format(pid,vx,vy,vz))
    fPress.write('{} {}\n'.format(pid,pr))
        
fCoord.close()
fVel.close()
fPress.close()

print('---| End writing initial condition')


