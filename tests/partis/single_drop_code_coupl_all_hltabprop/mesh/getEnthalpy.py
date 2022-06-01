#
# To use Cantera do:
# module load python/2.7.13
# source /gpfs/projects/bsc21/Cantera-Alya/cantera/INTEL/bin/setup_cantera
#

import cantera as ct

gas = ct.Solution('gri30.cti')

gas.TPX = 400.0, ct.one_atm, "N2:0.79, O2:0.21"

print(gas.report())

print('Enthalpy: {} J/(kgK)'.format(gas.enthalpy_mass))




