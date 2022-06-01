import sys
import os
import math
import numpy as np
import cantera as ct

meshName    = sys.argv[1]
mech        = sys.argv[2]
fuel        = sys.argv[3]
oxidizer    = sys.argv[4]
Tin         = float(sys.argv[5])
P           = float(sys.argv[6])

fuelBoundaries = sys.argv[7].split(',')
oxiBoundaries = sys.argv[8].split(',')



coordfile       = '{}.coord'.format(meshName)

os.system('rm -rf Fields')
os.system('mkdir Fields')

tempfile        = 'Fields/TEMPE.alya'
enthfile        = 'Fields/ENTHA.alya'
concePatt       = 'Fields/CON{:02d}.alya'



gasFuel                          = ct.Solution(mech)
gasFuel.TPY                      = (Tin,P,fuel)

gasOxidizer                      = ct.Solution(mech)
gasOxidizer.TPY                  = (Tin,P,oxidizer)

gasProdu                         = ct.Solution(mech)
gasProdu.TP                      = (Tin,P)
gasProdu.set_equivalence_ratio(1.0,fuel,oxidizer)
gasProdu.equilibrate('HP')

print(gasFuel.report())
print(gasOxidizer.report())
print(gasProdu.report())


fCoord          =open(coordfile,'r')
fTemp           =open(tempfile,'w')
fEnth           =open(enthfile,'w')
fCon            =[]

for ii in range(gasProdu.n_species):
    fCon.append( open(concePatt.format(ii+1),'w') )


print('---| Start writing initial condition')
for line in fCoord:
    data=line.split()

    pid = int(data[0])
    dims = len(data)-1
    x   = float(data[1])
    y   = float(data[2])


    H = gasProdu.enthalpy_mass
    T = gasProdu.T
    fEnth.write('{} {}\n'.format(pid,H))
    fTemp.write('{} {}\n'.format(pid,T))

    for ii in range(gasProdu.n_species):
        Y = gasProdu.Y[ii]
        fCon[ii].write('{} {}\n'.format(pid,Y))
        

fCoord.close()
fTemp.close()
for ii in range(gasProdu.n_species):
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
    for ii in range(gasProdu.n_species):
        f.write('  FIELD={}\n'.format(ii+4))
        f.write('    INCLUDE ./Fields/CON{:02d}.alya\n'.format(ii+1))
        f.write('  END_FIELD\n\n')

#
# Write filed size file
#
with open('field_size.dat','w') as f:
    f.write('  FIELDS={}\n'.format( 3 + len(gasProdu.Y)))
    f.write('    FIELD: 1, DIMEN= 2, NODES\n')
    f.write('    FIELD: 2, DIMEN= 1, NODES\n')
    f.write('    FIELD: 3, DIMEN= 1, NODES\n')
    for ii in range(gasProdu.n_species):
        f.write('    FIELD: {}, DIMEN= 1, NODES\n'.format(ii+4))
    f.write('  END_FIELDS\n')

#
# Write species boundary conditions
#
with open('speciesBoundary.dat','w') as f:
    for ii in range(gasProdu.n_species):
        f.write('  CODES, NODES, CLASS = {}\n'.format(ii+1))
        for bc in fuelBoundaries:
            f.write('{:>10}  1  {}\n'.format(bc, gasFuel.Y[ii]))
        for bc in oxiBoundaries:
            f.write('{:>10}  1  {}\n'.format(bc, gasOxidizer.Y[ii]))
        f.write('  END_CODES\n\n')


