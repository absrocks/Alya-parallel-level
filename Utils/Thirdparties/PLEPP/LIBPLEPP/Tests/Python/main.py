#!/usr/bin/env python
#=========================================================================||===#

#===================================================================| init |===#
import sys 
app_type = "pwrapper"
app_name = sys.argv[1]


#=================================================================| Mpi4py |===#
import os 
PATH = os.getcwd()
ROOT = "/home/jmake/z2014/Cplng/Mpi4python/All_Beta/"
sys.path.append(ROOT+"/Wrappers/Python")
sys.path.append("/home/jmake/z2014/ParaView410/lib/paraview-4.1/site-packages/") 


#=========================================================================||===#
#=========================================================================||===#
import Commdomm


from mpi4py import MPI 
#import MPI 

print (MPI.get_vendor())
# Open MPI all handles are pointers
# MPICH2 they are plain 'int'.

world_comm = MPI.COMM_WORLD
world_rank = world_comm.Get_rank()
world_size = world_comm.Get_size()


CD = Commdomm.CommDom() 
CD.init() 
CD.set_app_type(app_type);
CD.set_app_name(app_name);
CD.set_world_comm(world_comm)

local_comm = MPI.COMM_NULL
local_comm = CD.set_mpi_comms()
local_comm.Barrier()


local_rank = local_comm.Get_rank()

val = -1
if(local_rank==0): val = 69


commij = MPI.COMM_NULL
if(app_name=="CCCC"):
  namej     = "AAAA"
  commij    = CD.get_mpi_commij(namej); 
  comm_size = CD.get_mpi_commij_size();
  print

  n_recv = 0
  n_send = 1
  recv   = Commdomm.iarray(n_recv)
  send   = Commdomm.iarray(n_send)

  send[0] = val
  CD.sendrecv_int(send, n_send, recv, n_recv, local_comm, commij)
  val = local_comm.bcast(val, root=0)


if(app_name=="AAAA"):
  namej     = "CCCC"
  commij    = CD.get_mpi_commij(namej); 
  comm_size = CD.get_mpi_commij_size();
  print

  n_recv = 1
  n_send = 0
  recv   = Commdomm.iarray(n_recv)
  send   = Commdomm.iarray(n_send)

  recv[0] = -1
  CD.sendrecv_int(send, n_send, recv, n_recv, local_comm, commij)
  val = recv[0]
  val = local_comm.bcast(recv[0], root=0) 


if(app_name=="PYTH"):
  namej = "4TRN"
  commij = CD.get_mpi_commij(namej); 
  comm_size = CD.get_mpi_commij_size();
  print

  n_recv = 1
  n_send = 0
  recv   = Commdomm.iarray(n_recv)
  send   = Commdomm.iarray(n_send)

  recv[0] = -1
  CD.sendrecv_int(send, n_send, recv, n_recv, local_comm, commij)
  val = recv[0]
  val = local_comm.bcast(recv[0], root=0) 
  

local_comm.Barrier()
print " wrank/val p:", world_rank, val, 
if(local_rank==0): print "<--"


#====================================================================| End |===# 
#world_comm.Barrier() 
#world_comm.Free()

if world_rank==0:
  print "OK!"
  print 
