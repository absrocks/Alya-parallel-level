import sys
import math
import numpy as np


meshName=sys.argv[1]


coordfile       = '{}.coord'.format(meshName)
velofile        = 'VELOC.alya'
tempfile        = 'TEMPE.alya'
#enthfile        = 'ENTHA.alya'
#concePatt       = 'CON{:02d}.alya'

fCoord          =open(coordfile,'r')
fVel            =open(velofile,'w')
fTemp           =open(tempfile,'w')

c = 0.2000
pi = np.pi
tol = 0.001

print('---| Start writing initial condition')
for line in fCoord:
    data=line.split()

    pid = int(data[0])
    dims = len(data)-1
    x   = float(data[1])
    y   = float(data[2])
    #z   = float(data[3])

    vx = 1.0
    vy = 0.0
    #vz = 0.0
    #T = np.sin(2*pi*x/1)

    if float(data[1]) <= -(c+tol):
        T = 0.0
    elif float(data[1]) >= c+tol:
        T = 0.0
    else:
        T = 1.0

    fVel.write('{} {} {}\n'.format(pid,vx,vy))
    fTemp.write('{} {}\n'.format(pid,T))
        
fCoord.close()
fVel.close()
fTemp.close()

print('---| End writing initial condition')


