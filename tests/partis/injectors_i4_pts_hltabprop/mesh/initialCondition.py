import sys
import math
import numpy as np


meshName=sys.argv[1]


coordfile       = '{}.coord'.format(meshName)
concefile       = 'CONCE.alya'

fCoord          =open(coordfile,'r')
fCon            =open(concefile,'w')


print('---| Start writing initial condition')
for line in fCoord:
    data=line.split()

    pid = int(data[0])
    dims= len(data)-1
    x   = float(data[1])
    y   = float(data[2])
    z   = float(data[3])
    
    mf  = 0.0


    fCon.write('{} {} {} {} {}\n'.format(pid,0.0,0.0,mf,0.0))
        




fCoord.close()
fCon.close()
print('---| End writing initial condition')


