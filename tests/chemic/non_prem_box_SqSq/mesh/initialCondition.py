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

Zst = 0.055
rholeft = 0.66112

flamname = 'cam_diff_flamelet_i56.dat'
flamelet = np.loadtxt(flamname,skiprows=1)
with open(flamname,'r') as f:
    line = f.readline()
    cols = line.split()
    zCol = cols.index('z')
    ycCol = cols.index('Yc')
    rhoCol = cols.index('rho')
    correctOrd = np.argsort(flamelet[:,zCol])



print('---| Start writing initial condition')
for line in fCoord:
    data=line.split()

    pid = int(data[0])
    dims = len(data)-1
    x   = float(data[1])
    y   = float(data[2])
    z   = float(data[3])

    wavelength = math.pi * 10.0 / 1000.0 

    
    mf1 = Zst * (1.0 - x / wavelength)
    mf2 = (3.0*Zst) * (0.5 - 0.5*math.cos(2.0 * math.pi * x/wavelength))
    Z = mf1 + mf2
    Yc  = 0.4 * np.interp([Z], flamelet[:,zCol][correctOrd], flamelet[:,ycCol][correctOrd] )[0]

    vleft = ( 0.1 * (y+0.00785398)*(y-0.00785398)*(z+0.00785398)*(z-0.00785398) * 262809145.7  )

    vx = vleft 
    vy = 0.0
    vz = 0.0

    ZZ   = 0.1 * Z * (1-Z) + Z*Z
    YcYc = 0.2 * Yc**2     + Yc * Yc

    

    

    #--| H_f:                      -4672509.67135
    #--| H_ox:                     -10205.8103954
    H = Z * (-4672509.67135) + (1-Z) * (-10205.8103954)

    fVel.write('{} {} {} {}\n'.format(pid,vx,vy,vz))
    fTemp.write('{} {}\n'.format(pid,288.0))
    fEnth.write('{} {}\n'.format(pid,H))
    fCon01.write('{} {}\n'.format(pid,Yc))
    fCon02.write('{} {}\n'.format(pid,YcYc))
    fCon03.write('{} {}\n'.format(pid,Z))
    fCon04.write('{} {}\n'.format(pid,ZZ))
        




fCoord.close()
fVel.close()
fTemp.close()
fCon01.close()
fCon02.close()
fCon03.close()
fCon04.close()
print('---| End writing initial condition')


