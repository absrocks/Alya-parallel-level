#!/usr/bin/env python
import sys 
import numpy as np 

sys.path.append("/home/jmake/z2014/Alya/Repository/Tools/VTK_Develoment")
import vtk_multiblock 
import vtk_tools as tools


Ijk_max =  np.array([2,2,2])**9
Ijk_min =  np.array([0,0,0]) 
print Ijk_max


VTI = tools.__create_vti__(Ijk_max)
Ijk_max -= 1

import array

fileout = open("xxx.dat", "wb")
headerout = array.array('I')
headerout.append( Ijk_max[0] )
headerout.append( Ijk_max[1] )
headerout.append( Ijk_max[2] )
headerout.append( 1 )
headerout.tofile( fileout )

def linear(x, x0, y0, x1, y1):
    return y0 + (y1-y0)*(x-x0)/(x1-x0)


PATH = "/home/jmake/z2014/Alya/Runner/DaniCase/"
PATH+= "VELNORM01/"
for m in [8]: 
  CASE = PATH + "velNorm01_%d/velNorm01_%d_0_0.vtu" % (m, m) 
  print CASE
  
  VTU  = tools.__vtu_reader__(CASE)
  VTU  = VTU.GetOutput()
  
  DATA = tools.__get_pts_data__( VTU, "VelNorm")
  DATA = np.array(DATA)
  #DATA.tofile("velNorm01_%d.dat" % m)


  Pts = tools.__get_coords__(VTU) 
  Pts = np.array(Pts)
  Pts_min = Pts.min(axis=0)
  Pts_max = Pts.max(axis=0)
  #print Pts.min(axis=0), 
  #print Pts.max(axis=0)

  Pts[:,0] = linear(Pts[:,0], Pts_min[0], Ijk_min[0], Pts_max[0], Ijk_max[0])
  Pts[:,1] = linear(Pts[:,1], Pts_min[1], Ijk_min[1], Pts_max[1], Ijk_max[1])
  Pts[:,2] = linear(Pts[:,2], Pts_min[2], Ijk_min[2], Pts_max[2], Ijk_max[2])
  Pts = np.floor(Pts).astype(int)

  dataout = array.array('f')

  n_Pts = len(Pts)  
  for i in range(n_Pts):
    pt = Pts[i]
    VTI.SetScalarComponentFromFloat( pt[0],pt[1],pt[2],0, DATA[i]) 


  dims = VTI.GetDimensions()
  print dims 
  
  for i in range(dims[0]): 
    for j in range(dims[1]): 
      for k in range(dims[2]):
        val = VTI.GetScalarComponentAsDouble(k, j, i, 0)
        dataout.append( val )

  dataout.tofile(fileout)
  fileout.flush()

  tools.__save_vti__("velNorm01_%d" % m, VTI)


fileout.close()


print "OK! \n"
