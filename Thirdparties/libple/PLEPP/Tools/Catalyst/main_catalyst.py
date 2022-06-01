#===================================================================| init |===#
import sys
import os
PATH = os.getcwd()
sys.path.append("/home/jmake/z2014/Cplng/Mpi4python/All_Beta/Wrappers/Python")
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


#=================================================================| COMMs |===#
from paraview.simple import *
import vtkPVVTKExtensionsCorePython
import time

try: paraview.simple
except: from paraview.simple import *

from paraview import coprocessing
import vtkPVCatalystPython as catalyst #vtkCoProcessorPython


#pipeline = pythoncatalyst.vtkCPPythonScriptPipeline()


coProcessor = None


def CatalystInit():
  if(coProcessor==None):
    coProcessor = catalyst.vtkCPProcessor()

    #pipeline = vtkCPPythonScriptPipeline() 
    #coProcessor.AddPipeline(pipeline)


def CatalystFinalize():
  if(coProcessor!=None):
    coProcessor.Finalize()
  

def CatalystCoProcess(timeStep, time, grid): 
  datadescription = catalyst.vtkCPDataDescription()
  datadescription.AddInput("input")
  datadescription.SetTimeData(time, timeStep)

  #grid = vtk.vtkImageData()
  #array = vtk.vtkDoubleArray()
  #...
  #grid.GetPointData().AddArray(array)
  #...
  #
  #if(Processor->RequestDataDescription(dataDescription) != 0): 
    #inputdescription = datadescription.GetInputDescriptionByName("input").SetGrid(grid)
    ##inputdescription.SetGrid(pd)
  #
  #cpProcessor.CoProcess(dataDescription);
  #


print "OK!! \n"
