#!/usr/bin/python
import os
import re
import glob
import sys
import math 


#===============================================================================================# 
def Lineal(x, p1=[0.0,0.0], p2=[0.0,0.0]):
  y  = ( p2[1] - p1[1] )/( p2[0]-p1[0] )
  y *= ( x - p1[0] )
  y += p1[1]
  return y


def Read_alya_geo(fname, INIT_KEY="COORDINATES"):
    fname = glob.glob(fname) 
    if(not len(fname)<1):
      fname = fname[-1]
      print "|_\'%s\'" % fname
    else:
      print "Error: \'%s\'! \n" % fname 
      exit(1)

    data = open(fname, "r")
    lines = data.readlines()
    data.close()

    nline = len(lines)

    #global INIT_KEY
    INIT_KEY = INIT_KEY.replace("_", "") 
    INIT_KEY = INIT_KEY.replace("-", "") 
    INIT_KEY = INIT_KEY.replace("&", "") 
    END_KEY  = "END" + INIT_KEY 

    ok  = False 
    IDs = []
    for i in range(nline):
      line = lines[i]
      if(not line.find(INIT_KEY)<0): IDs.append(i+1)
      if(not line.find(END_KEY)<0):  IDs.append(i+0) 

    IDX = []
    XYZ = [] 
    for i in range(IDs[0], IDs[1]-1):
      line = lines[i]
      line = line.strip()
      line = line.split()
      IDX.append( eval(line[0])-1 )
      XYZ.append( -666.66 )

    n_ids = len(IDX) 
    j = 0
    for i in range(IDs[0], IDs[1]-1):
      line = lines[i]
      line = line.strip()
      line = line.split()

      idx  = IDX[j]
      j   += 1 
      XYZ[idx] = [eval(val) for val in line[1:]]
      #print line

    #XYZ = []
    #for i in range(IDs[0], IDs[1]-1):
    #  line = lines[i]
    #  line = line.strip() 
    #  line = line.split()
    #  XYZ.append([eval(val) for val in line[1:]]) 
 
    print "  |_No elements:", len(XYZ)
    return XYZ 



#==============================================================================#
#==============================================================================#
#--------------------------------------------------------------------------||--#
ALYA_PTs = "../banc_2d_Test_1/../mesh_2d/fort.777"
PTs  = Read_alya_geo(ALYA_PTs, "COORDINATES") 
n_pts = len(PTs)
#print PTs 


ALYA_ELEMs = "../banc_2d_Test_1/../mesh_2d/fort.777"
ELEMs  = Read_alya_geo(ALYA_ELEMs, "ELEMENTS") 
n_elems = len(ELEMs)

ALYA_TYPES = "../banc_2d_Test_1/../mesh_2d/fort.777"
TYPES  = Read_alya_geo(ALYA_TYPES, "TYPES") 
n_types = len(TYPES)
#--------------------------------------------------------------------------||--#

#--------------------------------------------------------------------------||--#
#VTK_PATH = "/home/jmake/z2014/Cplng/Mpi4python/Old/PLEPP_2014Jun16/Tools/VTK"
#sys.path.insert(0, VTK_PATH)
import vtk_creator as Creator

VTU = Creator.Vtk_creator() 

for i in range(n_pts): 
  coord = PTs[i]
  x = coord[0]
  y = coord[1]
  z = 0.0
  VTU.set_point([x,y,z]) 

CelDic = Creator.Vtk2Alya_cell()
for i in range(n_elems):
  alya_cel = TYPES[i][0]
  vtk_cel  =  CelDic.get_vtk_cell(alya_cel) 

  elem = ELEMs[i]
  cel  = [vrtx-1  for vrtx in elem]
  VTU.set_cell(cel, vtk_cel[2])
#--------------------------------------------------------------------------||--#

#--------------------------------------------------------------------------||--#
#VTU.set_points_prop_scalar(W, "W")
VTU.set_vtk_unstructured() 
VTU.save_vtk_unstructured("geo_mesh")
#--------------------------------------------------------------------------||--#

#===============================================================================================#
print "OK!! \n\n"
