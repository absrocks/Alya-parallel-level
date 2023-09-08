#!/usr/bin/env python
import sys 
import numpy as np 

sys.path.append("/home/jmake/z2014/Alya/Repository/Tools/VTK_Develoment")
import vtk_multiblock 
import vtk_tools as tools

#Things = ["VelNorm", "Ryz"]
Things = ["VELOC", "VelNorm", "xxx"]

PATH = "/home/jmake/z2014/Alya/Runner/DaniCase/"
PATH+= "VELNORM01/"

for m in range(60, 100): 
  CASEin = PATH + "velNorm01_%d/velNorm01_%d_0_0.vtu" % (m, m) 
  
  VTU     = tools.__vtu_reader__(CASEin)
  n_Props = VTU.GetOutput().GetPointData().GetNumberOfArrays()
  
  Props = []
  for i in range(n_Props): 
    prop = VTU.GetOutput().GetPointData().GetArrayName(i)
    Props.append(prop) 

  for think in Things: 
    if(think in Props): Props.remove(think)

  for prop in Props:
    VTU.GetOutput().GetPointData().RemoveArray( prop )
    print "Deleted: \'%s\'" % prop 

  PATH02    = "/home/jmake/z2014/Alya/Runner/DaniCase/" 
  CASEout = PATH02 + "velNorm01_%s.vtu" % ( str(m).zfill(3) ) 
  tools.__save_vtu__(CASEout, VTU.GetOutput())

  #StreamTracer = tools.__vtkStreamTracer__( VTU.GetOutputPort() ) 
  #CASEout = PATH + "StreamTracer01_%s.vtu" % ( str(m).zfill(3) ) 
  #tools.__save_vtp__(CASEout, StreamTracer.GetOutput())  
  
  

print "OK! \n"
