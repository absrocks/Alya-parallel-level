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


rho0 = 0.957753
rho1 = 0.121074


print('---| Start writing initial condition')
for line in fCoord:
    data=line.split()

    pid = int(data[0])
    dims = len(data)-1
    x   = float(data[1])
    y   = float(data[2])
    z   = float(data[3])

    wavelength = math.pi * 10.0 / 1000.0 

    if x < 0.5*wavelength:
        C = (0.5 - 0.5*math.cos(2.0 * math.pi * x/wavelength))
    else:
        C = 1.0

    rho = rho0 + C * (rho1-rho0)
    vleft = 0.0


    vx = vleft #* rho0 / rho
    vy = 0.0
    vz = 0.0

    CC   = 0.05 * (C*(1.0-C)) + C*C
    H    = -1.75179161e+05

    fVel.write('{} {} {} {}\n'.format(pid,vx,vy,vz))
    fTemp.write('{} {}\n'.format(pid,288.0))
    fEnth.write('{} {}\n'.format(pid,H))
    fCon01.write('{} {}\n'.format(pid,C))
    fCon02.write('{} {}\n'.format(pid,CC))
        


fCoord.close()
fVel.close()
fTemp.close()
fCon01.close()
fCon02.close()
print('---| End writing initial condition')


