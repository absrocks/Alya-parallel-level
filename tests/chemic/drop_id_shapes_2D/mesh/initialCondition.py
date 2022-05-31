import sys
import math
import numpy as np


meshName=sys.argv[1]


coordfile       = '{}.coord'.format(meshName)
velofile        = 'VELOC.alya'
tempfile        = 'TEMPE.alya'
enthfile        = 'ENTHA.alya'
concePatt       = 'CON{:02d}.alya'

fCoord          =open(coordfile,'r')
fVel            =open(velofile,'w')
fTemp           =open(tempfile,'w')
fEnth           =open(enthfile,'w')
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
#3D#    z   = float(data[3])

#
#  Dry air enthalpies:
#  298.0  K: -112.326666559
#  373.15 K: 76094.4256129
#  673.15 K: 389117.295797 
# 
    
    T = 373.15
    H = 76094.4256129
    
    vx = 0.0 
    vy = 0.0
#3D#    vz = 0.0


#3D#    fVel.write('{} {} {} {}\n'.format(pid,vx,vy,vz))
    fVel.write('{} {} {}\n'.format(pid,vx,vy))
    fTemp.write('{} {}\n'.format(pid,T))
    fEnth.write('{} {}\n'.format(pid,H))
    fCon01.write('{} {}\n'.format(pid,0.0))
    fCon02.write('{} {}\n'.format(pid,0.0))
    fCon03.write('{} {}\n'.format(pid,0.0))
    fCon04.write('{} {}\n'.format(pid,0.0))
        




fCoord.close()
fVel.close()
fTemp.close()
fCon01.close()
fCon02.close()
fCon03.close()
fCon04.close()
print('---| End writing initial condition')


