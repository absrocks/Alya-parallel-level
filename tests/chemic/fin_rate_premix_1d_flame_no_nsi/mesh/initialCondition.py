import sys
import os
import math
import numpy as np
import cantera as ct

meshName    = sys.argv[1]
mech        = sys.argv[2]
fuel        = sys.argv[3]
oxidizer    = sys.argv[4]
phi         = float(sys.argv[5])
Tleft       = float(sys.argv[6])
vleft       = float(sys.argv[7])
P           = float(sys.argv[8])

boundaries = []
for ii in range(9,len(sys.argv)):
    boundaries.append(sys.argv[ii])



coordfile       = '{}.coord'.format(meshName)

os.system('rm -rf Fields')
os.system('mkdir Fields')

tempfile        = 'Fields/TEMPE.alya'
enthfile        = 'Fields/ENTHA.alya'
concePatt       = 'Fields/CON{:02d}.alya'



gasReact                         = ct.Solution(mech)
gasReact.TP                      = (Tleft,P)
gasReact.set_equivalence_ratio(phi,fuel,oxidizer)
gasProdu                         = ct.Solution(mech)
gasProdu.TP                      = (Tleft,P)
gasProdu.set_equivalence_ratio(phi,fuel,oxidizer)
gasProdu.equilibrate('HP')

print(gasReact.report())
print(gasProdu.report())


fCoord          =open(coordfile,'r')
fTemp           =open(tempfile,'w')
fEnth           =open(enthfile,'w')
fCon            =[]

for ii in range(gasReact.n_species):
    fCon.append( open(concePatt.format(ii+1),'w') )


print('---| Start writing initial condition')
for line in fCoord:
    data=line.split()

    pid = int(data[0])
    dims = len(data)-1
    x   = float(data[1])
    y   = float(data[2])

    if x < 0.000001:
        weight = 0.0
    else:
        weight = 1.0
    

    H = gasReact.enthalpy_mass
    fEnth.write('{} {}\n'.format(pid,H))
    fTemp.write('{} {}\n'.format(pid,298.0))

    for ii in range(gasReact.n_species):
        YR = gasReact.Y[ii]
        YP = gasProdu.Y[ii]
        Y  = YR*(1.0-weight) + YP*weight
        fCon[ii].write('{} {}\n'.format(pid,Y))
        

fCoord.close()
fTemp.close()
for ii in range(gasReact.n_species):
    fCon[ii].close()
print('---| End writing initial condition')


#
# Write fileds file
#
with open('fields.dat','w') as f:
    f.write('  FIELD={}\n'.format(2))
    f.write('    INCLUDE ./Fields/ENTHA.alya\n')
    f.write('  END_FIELD\n\n')

    f.write('  FIELD={}\n'.format(3))
    f.write('    INCLUDE ./Fields/TEMPE.alya\n')
    f.write('  END_FIELD\n\n')
    for ii in range(gasReact.n_species):
        f.write('  FIELD={}\n'.format(ii+4))
        f.write('    INCLUDE ./Fields/CON{:02d}.alya\n'.format(ii+1))
        f.write('  END_FIELD\n\n')

#
# Write filed size file
#
with open('field_size.dat','w') as f:
    f.write('  FIELDS={}\n'.format( 3 + len(gasReact.Y)))
    f.write('    FIELD: 1, DIMEN= 2, NODES\n')
    f.write('    FIELD: 2, DIMEN= 1, NODES\n')
    f.write('    FIELD: 3, DIMEN= 1, NODES\n')
    for ii in range(gasReact.n_species):
        f.write('    FIELD: {}, DIMEN= 1, NODES\n'.format(ii+4))
    f.write('  END_FIELDS\n')

#
# Write species boundary conditions
#
with open('speciesBoundary.dat','w') as f:
    for ii in range(gasReact.n_species):
        f.write('  CODES, NODES, CLASS = {}\n'.format(ii+1))
        for bc in boundaries:
            f.write('{:>10}  1  {}\n'.format(bc, gasReact.Y[ii]))
        f.write('  END_CODES\n\n')


