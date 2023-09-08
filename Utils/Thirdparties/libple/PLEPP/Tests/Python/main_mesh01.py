#!/usr/bin/env python
#
#=================================================================| Mpi4py |===#
from mpi4py import MPI 
world_comm = MPI.COMM_WORLD
world_rank = world_comm.Get_rank()
world_size = world_comm.Get_size()

#===================================================================| init |===#
import sys
import os 
PATH = os.getcwd()
ROOT = "/home/jmake/z2014/Cplng/Mpi4python/All_Beta/"
sys.path.append(ROOT+"/Wrappers/Python")

import Commdomm
import Loader_alya
import Read_file 

Types     = ["COORDINATES", "ELEMENTS", "CHARACTERISTICS"]

## 3d case, coinciding meshes, Temper 3-interaction
basenamei = ROOT+"/Tools/Temper/Mesh01/METIS02/xxx_%s%s.alya"
filenamej = ROOT+"/Tools/Temper/Mesh02/METIS01/xxx_COORDINATES000.alya" 


#=================================================================| COMMs |===#
app_type = "PYTH"
app_name = "MESH01"

CD = Commdomm.CommDom() 
CD.init() 
CD.set_app_type(app_type);
CD.set_app_name(app_name);
CD.set_world_comm(world_comm)

local_comm = MPI.COMM_NULL
local_comm = CD.set_mpi_comms()
local_rank = local_comm.Get_rank()
local_size = local_comm.Get_size()


rank = str( local_rank ).zfill(3) 

DATA = {}
for typei in Types: 
  filename = basenamei%(typei, rank) 

  P = Read_file.read_log_file()
  P.set_name(filename)
  P.run()
  n_cols = P.get_cols()
  n_rows = P.get_rows()

  Data = []
  for i in range(n_rows): 
    data = []
    for j in range(1,n_cols): data.append( P.get_data_ij(i, j) ) 
    Data.append( data ) 
  P.end()
  DATA[typei] = Data


Elements = DATA["ELEMENTS"]
n_elements = len(Elements) 

vertex_num = []
for element in Elements:
  for node in element: vertex_num.append(int(node)-0) 
  #vertex_num.append(-169)


Coords = DATA["COORDINATES"]
n_vertices = len(Coords) 

vertex_coords = []
for coord in Coords:
  for pt in coord: vertex_coords.append(pt)


Q = Read_file.read_log_file()
Q.set_name(filenamej)
Q.run()
n_cols = Q.get_cols()
n_rows = Q.get_rows()

Dataj = []
for i in range(n_rows): 
  data = []
  for j in range(1,n_cols): data.append( Q.get_data_ij(i, j) ) 
  Dataj.append( data ) 
Q.end()


n_vertices_j    = 0
vertex_coords_j = Commdomm.darray( 0 )


if(local_rank==0):
  n_vertices_j    = len( Dataj )
  vertex_coords_j = Commdomm.darray( n_vertices_j*3 )

  k=0
  for dataj in Dataj:
    for pt in dataj: 
      vertex_coords_j[k] = pt
      k += 1

local_comm.Barrier()


#===============================================================| LOCATOR |===#
val = -1
if(local_rank==0): val = 69
val = local_comm.bcast(val, root=0) 

comm_size = CD.get_mpi_commij_size();

commij = MPI.COMM_NULL  
if(comm_size == 1):
  namej     = "MESH02"
  commij    = CD.get_mpi_commij(namej)
  print

  # sendrecv test 
  n_recv = 0
  n_send = 1
  recv   = Commdomm.iarray(n_recv)
  send   = Commdomm.iarray(n_send)

  send[0] = val
  CD.__mpi_sendrecv__(send, n_send, recv, n_recv, local_comm, commij)

  # init locator
  rooti = Commdomm.iarray(1)
  rootj = Commdomm.iarray(1)
  CD.__mpi_get_roots__(local_comm, commij, rooti, rootj); 

  n_rank    = local_size
  root_rank = rooti[0]
  CD.locator_create(commij, n_rank, root_rank, 1e-3)

  # init meshes
  n_vertices_i    = n_vertices
  vertex_coords_i = Commdomm.darray( len(vertex_coords) )
  for i in range( len(vertex_coords) ): vertex_coords_i[i] = vertex_coords[i]

  n_elements_i    = n_elements
  vertex_num_i    = Commdomm.iarray( len(vertex_num) )
  for i in range( len(vertex_num) ): vertex_num_i[i] = vertex_num[i]

  CD.locator_set_mesh(n_vertices_i, n_elements_i, vertex_coords_i, vertex_num_i, 
                      n_vertices_j, vertex_coords_j)

  CD.exchange(world_rank); 
  local_comm.Barrier()

  CD.save_dist_coords(local_rank)

  #local_comm.Barrier()
  #if(local_rank==0):print "OK!!\n"
#====================================================================| End |===# 
