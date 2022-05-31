#!/usr/bin/env python
import sys
import os 

#==================================================================| Parse |===#
from optparse import OptionParser
usage = "\n\t!! >> python ./%prog [options] arg1 arg2!!\n"
parser = OptionParser(usage=usage)
parser.add_option("-n", "--name", type="string", nargs=1,
                  dest="app_name", default=(''), help="application name")
parser.add_option("-m", "--mesh", type="string", nargs=1,
                  dest="mesh_path", default=('./'), help="mesh path")
parser.add_option("-f", "--files", type="string", nargs=3,
                  dest="files", default=('', '', ''), help="files")
(options, args) = parser.parse_args()


if(options.app_name=='' or 
   options.files[0]=='' or 
   options.files[1]=='' or 
   options.files[2]==''): 

  parser.print_help()
  print 
  sys.exit() 
#-------------------------------------------------------------------------||---#

#===================================================================| init |===#
#-------------------------------------------------------------------------||---#
PATH = os.getcwd()
ROOT = "/home/jmake/z2014/Alya/Repository/Alya2014Ago05/"
sys.path.append(ROOT+"/Thirdparties/libple/PLEPP/Wrappers/Python")

namei  = "XXX"
namej  = "YYY"
#-------------------------------------------------------------------------||---#


#=================================================================| Mpi4py |===#
from mpi4py import MPI 
world_comm = MPI.COMM_WORLD
world_rank = world_comm.Get_rank()
world_size = world_comm.Get_size()
#-------------------------------------------------------------------------||---#

#==================================================================| COMMs |===#
import Commdomm

app_type = "SYRTHES 4"
app_name = options.app_name

CD = Commdomm.CommDom() 
CD.init() 
CD.set_app_type(app_type);
CD.set_app_name(app_name);
CD.set_world_comm(world_comm)

local_comm = MPI.COMM_NULL
local_comm = CD.set_mpi_comms()
local_rank = local_comm.Get_rank()
local_size = local_comm.Get_size()
#-------------------------------------------------------------------------||---#

#===================================================================| MESH |===#

#-------------------------------------------------------------------------||---#
n_vertices_j    = 0
vertex_coords_j = Commdomm.darray( n_vertices_j )

n_vertices    = 0
vertex_coords = []

n_elements = 0
vertex_num = [] 
#-------------------------------------------------------------------------||---#

#-------------------------------------------------------------------------||---#
sys.path.append(ROOT+"/Thirdparties/libple/PLEPP/Tools/VTK")
import vtk_creator 
VTKout   = vtk_creator.Vtk_creator() 
Alya2vtk = vtk_creator.Vtk2Alya_cell() 
#-------------------------------------------------------------------------||---#


#-------------------------------------------------------------------------||---#
from read_alya_geo import Read_alya_geo, Get_cell_type
Alya2Cs = Get_cell_type()

ALYA_PTS   = options.mesh_path + options.files[0]
ALYA_ELEMs = options.mesh_path + options.files[1]
ALYA_TYPES = options.mesh_path + options.files[2]
ALYA_VTU   = ''

PTS     = Read_alya_geo(ALYA_PTS, "COORDINATES")
n_pts   = len(PTS)

ELEMs   = Read_alya_geo(ALYA_ELEMs, "ELEMENTS")
n_elems = len(ELEMs)

TYPES   = Read_alya_geo(ALYA_TYPES, "TYPES")
n_types = len(TYPES)
#-------------------------------------------------------------------------||---#

#-------------------------------------------------------------------------||---#
if(1==2):
  fout = open("xxx.dat", "w")
  print>> fout, n_elems 
  for i in range(n_elems):
    Elem   = ELEMs[i]
    n_Elem = len(Elem)
    for j in range((n_Elem)): print>> fout, Elem[j], 
    print>> fout, ""
  fout.close() 


import numpy as np
if(local_size>1 and 1==1):
  gIDs = [ [] for i in range(local_size) ] 
  parts = np.genfromtxt("xxx.dat.epart.%d" % local_size) 

  for i in range(n_elems): 
    if(parts[i]==local_rank): 
      gIDs[local_rank].append(i)

  lIDs   = gIDs[local_rank]
  n_lIDs = len(lIDs) 

  ELEMi  = [] 
  TYPESi = [] 
  for i in range(n_lIDs): 
    IDi = lIDs[i]
    ELEMi.append(  ELEMs[IDi] ) 
    TYPESi.append( TYPES[IDi] )
    
  ELEMs   = ELEMi 
  TYPES   = TYPESi
  n_elems = n_lIDs
#-------------------------------------------------------------------------||---#

#-------------------------------------------------------------------------||---#
vertex_coords = []
for i in range(n_pts):
  coord = PTS[i]
  VTKout.set_point(coord) 
  for pt in coord: vertex_coords.append(pt)


vertex_num  = []
vertex_type = []
for i in range(n_elems):
  element = ELEMs[i]
  for node in element: vertex_num.append(int(node)-0) 

  alya_type = TYPES[i][0]
  cs_type = Alya2Cs.get_saturne_cell(alya_type)
  vertex_type.append(cs_type)   

  vtk_cell  = [elemi-1 for elemi in element] 
  vtk_type  = Alya2vtk.get_vtk_cell(alya_type)  
  VTKout.set_cell(vtk_cell, vtk_type[2]) 

n_vertices = n_pts
n_elements = n_elems

print "|_n_vertices, n_elements: %d %d " %(n_vertices, n_elements)
print "                        |_types:", Alya2Cs.print_types() 
#-------------------------------------------------------------------------||---#

#================================================================| LOCATOR |===#
#-------------------------------------------------------------------------||---#
n_vertices_i    = n_vertices  
vertex_coords_i = Commdomm.darray( len(vertex_coords) )
for i in range( len(vertex_coords) ): vertex_coords_i[i] = vertex_coords[i]

n_elements_i    = n_elements
vertex_num_i    = Commdomm.iarray( len(vertex_num) )
for i in range( len(vertex_num) ): vertex_num_i[i] = vertex_num[i]

vertex_type_i   = Commdomm.iarray( n_elements )
for i in range( n_elements ): vertex_type_i[i] = vertex_type[i]
#-------------------------------------------------------------------------||---#

#-------------------------------------------------------------------------||---#
if(local_rank==0):
  n_vertices_j    = len(PTS)
  vertex_coords_j = Commdomm.darray( n_vertices_j*3 )

  k=0 
  for i in range(n_vertices_j):
    coord = PTS[i]
    for j in range(3): 
      vertex_coords_j[k] = coord[j]
      k += 1

local_comm.Barrier()
#-------------------------------------------------------------------------||---#


#-------------------------------------------------------------------------||---#
commij = MPI.COMM_NULL  
if( (CD.__get_app_name__() == namei) and (CD.__get_friends__(namej) == 1) ):
  commij = CD.get_mpi_commij(namej)
if( (CD.__get_app_name__() == namej) and (CD.__get_friends__(namei) == 1) ):
  commij = CD.get_mpi_commij(namei)

CD.locator_create2(local_comm, commij, 0.001)

CD.__locator_set_meshj__(n_vertices_j, vertex_coords_j) 
CD.__create_nodal_t__(vertex_coords_i, n_vertices_i, 
                      vertex_num_i, vertex_type_i,  n_elements_i)  
CD.save_dist_coords(local_rank)
#-------------------------------------------------------------------------||---#

#-------------------------------------------------------------------------||---#
import math 
Propi = []

if(app_name=="XXX"):
  for i in range(n_pts):
    Coord =  PTS[i]
    mod   =  1.0*math.sqrt( sum([coord*coord for coord in Coord]) ) 
    Propi.append( mod )
if(app_name=="YYY"):
  for i in range(n_pts):
    Coord =  PTS[i]
    mod   = -1.0*math.sqrt( sum([coord*coord for coord in Coord]) ) 
    Propi.append( mod )
#-------------------------------------------------------------------------||---#

#-------------------------------------------------------------------------||---#
local_comm.Barrier()

n_recv = CD.get_n_interior() 
n_send = CD.get_n_dist_points()

dist_locations_i = Commdomm.iarray( n_send   )
CD.__locator_get_dist_locations__( dist_locations_i ) 

dist_coords_j    = Commdomm.darray( n_send*3 )
CD.__locator_get_dist_coords__(    dist_coords_j    ) 

parent_element_i = Commdomm.iarray( n_elems )
CD.__get_parent__(parent_element_i, n_elems)

point_j  = Commdomm.darray(3)
shapef_j = Commdomm.darray(8)
prop_i   = Commdomm.darray(8)

CellTouched     = [-1]*n_elems
vols_coords = []
var_ij      = Commdomm.darray(n_send) 
for i in range(n_send):
  point_j[0] = dist_coords_j[i*3+0]
  point_j[1] = dist_coords_j[i*3+1]
  point_j[2] = dist_coords_j[i*3+2]

  ielem              = parent_element_i[ dist_locations_i[i]-1 ]-1  #<- C-style
  CellTouched[ielem] = 1 

  alya_type    = TYPES[ielem][0]
  cs_type      = Alya2Cs.get_saturne_cell(alya_type)

  element      = ELEMs[ielem]
  n_nodes      = len(element)

  vertices_i   = Commdomm.iarray(n_nodes) 
  for j in range(n_nodes): 
    vertices_i[j] = element[j]               #<- Fortran-style 
    prop_i[j]     = Propi[ element[j]-1 ]
  #"""
  CD.__interpolation__( cs_type,
                        1.0e-1, 
                        vertex_coords_i, 
                        vertices_i, 
                        point_j, 
                        shapef_j )
  #"""
  if(cs_type==4): 
    CD.__tetra_interpolation__(vertex_coords_i, vertices_i, point_j, shapef_j)

  vols_coords.append(shapef_j) 
  #"""
  propj = sum([ prop_i[k]*shapef_j[k] for k in range(n_nodes)  ])
  var_ij[i] = propj 
#"""
#-------------------------------------------------------------------------||---#


#-------------------------------------------------------------------------||---#
var_ji = Commdomm.darray(n_recv) 
CD.__locator_exchange_double_scalar__(var_ij, var_ji)

interior_list_j  = Commdomm.iarray( n_recv ) 
CD.__locator_get_interior_list__(  interior_list_j  ) 

Propj = [0 for i in range(n_pts)]
for i in range( n_recv ): Propj[ interior_list_j[i]-1 ] = var_ji[i]

PtsTouched = [-1]*n_pts
for i in range( n_recv ): PtsTouched[ interior_list_j[i]-1 ] = 1
#-------------------------------------------------------------------------||---#


#-------------------------------------------------------------------------||---#
import sys
ALYA_VTU2 = options.app_name
VTKout.set_points_prop_scalar(PtsTouched, "PtsTouched")
VTKout.set_cell_prop_scalar(CellTouched, "CellTouched")
VTKout.set_vtk_unstructured() 
VTKout.psave(ALYA_VTU2, local_size, local_rank) 
#-------------------------------------------------------------------------||---#



#-------------------------------------------------------------------------||---#
#-------------------------------------------------------------------------||---#


#=========================================================================||===# 
#=========================================================================||===# 


#====================================================================| End |===# 
local_comm.Barrier()
if(local_rank==0):print "OK!!\n"
#=========================================================================||===# 
