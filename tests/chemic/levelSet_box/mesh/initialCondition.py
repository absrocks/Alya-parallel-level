import sys
import math
import numpy as np


meshName=sys.argv[1]


coordfile       = '{}.coord'.format(meshName)
velofile        = 'VELOC.alya'
concePatt       = 'CON{:02d}.alya'

fCoord          =open(coordfile,'r')
fVel            =open(velofile,'w')
fCon01          =open(concePatt.format(1),'w')
fCon02          =open(concePatt.format(2),'w')
fCon03          =open(concePatt.format(3),'w')
fCon04          =open(concePatt.format(4),'w')


print('---| Start writing initial condition')
for line in fCoord:
    data=line.split()

    pid = int(data[0])
    dims = len(data)-1
    x   = float(data[1])
    y   = float(data[2])
    z   = float(data[3])

    wavelength = math.pi * 1.0 / 1000.0 / 2.0


    if x < wavelength:
        Z=0.0
    else:
        Z =  (0.5 - 0.5*math.cos(2.0 * math.pi * x/wavelength))



    vx = ( 0.1 * (y+0.000785398)*(y-0.000785398)*(z+0.000785398)*(z-0.000785398) * 2628091457000.0  ) 
    vy = 0.0
    vz = 0.0

    fVel.write('{} {} {} {}\n'.format(pid,vx,vy,vz))
    fCon01.write('{} {}\n'.format(pid,0.0))
    fCon02.write('{} {}\n'.format(pid,0.0))
    fCon03.write('{} {}\n'.format(pid,Z))
    fCon04.write('{} {}\n'.format(pid,0.0))
        




fCoord.close()
fVel.close()
fCon01.close()
fCon02.close()
fCon03.close()
fCon04.close()
print('---| End writing initial condition')


