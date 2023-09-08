import sys
file_name = sys.argv[0]

sys.path.append("../../Geometries")
import standar_box_regular_mesh02 as mesh 

obj = mesh.Negociante


#------------------------------------------------------------| Salome2dat |---# 
sys.path.append("../../Salome2Alya")
from salome2dat import *

S = Salome2dat()
S.get_data(obj)
#S.save_data(file_name)


#-----------------------------------------------------| Initial Conditions |---# 
"""
Root: Re = 1.15e3 
"""
Raire          = 287.0 # N m/(kg K) (aire)
epsilon_aire   = 1.74e-5 # [Pa s]
epsilon_agua   = 1.0e-2
capacity_ratio = 7.0/5 

R0  = 1.2 
T0  = 300.0
P0  = 101000.0 
V0  = 100.0 

L0  = 1.0   
vis = epsilon_aire*5e3
print "vis: %f " % vis 
print "Re:  %fx1e3"% (V0*L0/vis*1e-3) 

Npts = S.get_npoints()
Pts  = S.get_points() 
G    = S.get_groups() 

Vel = np.zeros(Pts.shape)
Rho = np.zeros(Npts)
T   = np.zeros(Npts)

Vel[:,2] = V0
Rho[:]   = R0
T[:]     = T0 

#from vortex import Vortex
#Vrx = Vortex() 
#Vrx.a = 0.1  
#Vrx.vel(Vel, Pts)

from vortex2 import Vortex
Vrx = Vortex() 
Vrx.sound(T0, capacity_ratio)
print "Mach:", V0/Vrx.a 



Vrx.Ro = L0*0.1 
Vrx.C  = 1.0 # m**2/s
Vrx.vel(Vel, Pts)

Vrx.density(Rho, R0) 
Vrx.temp(T, T0, capacity_ratio) 

S.VTK.set_points_prop_vec(Vel, "Vel")
S.VTK.set_points_prop(Rho, "Rho")
S.VTK.set_points_prop(T, "T")
S.to_vtk(file_name, "faces")
S.to_vtk(file_name, "vols")


print np.average(Vel, axis=0)

#------------------------------------------------------------------| alya |---# 
S.Alya.NUMBER_OF_STEPS = 2000
S.Alya.set_prop_vals([0.0,0.0,0.0], "VELOCITY")
S.Alya.set_prop_vals(T0, "TEMPERATURE")

S.Alya.set_prop_vec(Vel, file_name, "VELOCITY") 
S.Alya.set_prop(Rho, file_name, "DENSITY") 
S.Alya.set_prop(T, file_name, "TEMPERATURE") 
## PRESS/DENS
##S.Alya.set_prop(P, file_name, "PRESSURE") # error!! 
## TEMPERATURE/ENERGY
##S.Alya.set_prop(T, file_name, "ENERGY") 

init_keys = {}
#init_keys["Inlet"]  = ["R", "T"] 
#init_keys["Outlet"] = ["R", "T"]
#init_keys["Walls"]  = ["V"]
init_keys["Walls"]  = ["n"]

init_vals = {}
#init_vals["Inlet"]  = [1.0*R0, T0] 
#init_vals["Outlet"] = [0.9*R0, T0] 
#init_vals["Walls"]  = [[0.0,0.0,0.0]] 
init_vals["Walls"]  = [0.0] 

S.periodic_boundaries_z("Outlet", "Inlet", file_name, 1e-5) 
S.to_alya(file_name)
S.Alya.set_initial_conditions_codes(S.groups_names, init_keys, init_vals)



#-----------------------------------------------------------------------------# 
print "OK! \n"
