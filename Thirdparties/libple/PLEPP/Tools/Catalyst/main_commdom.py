#!/usr/bin/env python
#
#===================================================================| init |===#
import sys
import os
PATH = os.getcwd()
ROOT = "/home/jmake/z2014/Cplng/Mpi4python/All_Beta/"
sys.path.append(ROOT+"/Wrappers/Python")
sys.path.append("/home/jmake/z2014/ParaView410/lib/paraview-4.1/site-packages/")


#=================================================================| Mpi4py |===#
from mpi4py import MPI 
print MPI.get_vendor()

world_comm = MPI.COMM_WORLD
world_rank = world_comm.Get_rank()
world_size = world_comm.Get_size()

Datatype =  MPI.Datatype


#=================================================================| COMMs |===#
import Commdomm

app_type = "PYTH"
app_name = "MESH01"

CD = Commdomm.CommDom() 
CD.init() 
CD.set_app_type(app_type);
CD.set_app_name(app_name);
CD.set_world_comm(world_comm)

local_comm = MPI.COMM_NULL
local_comm = CD.set_mpi_comms()
#local_rank = local_comm.Get_rank()
#local_size = local_comm.Get_size()





print "OK!! \n"
#========================================================================||===#
#========================================================================||===#
