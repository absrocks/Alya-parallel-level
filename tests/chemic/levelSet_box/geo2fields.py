#!/usr/bin/python
#import numpy as np
import sys
import os
import glob
import math

#===============================================================================================# 
INIT_KEY  = "COORDINATES"
END_KEY   = "END_" + INIT_KEY

def Read_alya_geo(fname):
    data = open(fname, "r")
    lines = data.readlines()
    data.close()

    nline = len(lines)
    global INIT_KEY
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

    XYZ = []      
    for i in range(IDs[0], IDs[1]-1):
      line = lines[i]
      line = line.strip() 
      line = line.split()
      XYZ.append([eval(val) for val in line[1:]]) 
    
    print "  |_No elements:", len(XYZ)
    return XYZ 


def Write_file(fname, data, dime=1):
  ndata = len(data)
  fdata = open(fname, "w") 
  for i in range(ndata): 
    line = data[i] 
    print>> fdata, i+1, 
    if(dime>1):
      for j in range(dime): print>> fdata, line[j],  
    else: 
      print>> fdata, line, 
    print>> fdata
    
#===============================================================================================# 

FILE = sys.argv[1]
pts  = Read_alya_geo(FILE) 

VELOC = []
CON01 = []
CON02 = []
CON03 = []
CON04 = []

for pt in pts:
      x = pt[0]
      if ( x < 0.015 and x > 0.01 ):
        CON03.append( 1.0 )
      else:
        CON03.append( 0.0 )
      CON01.append( 0.0 )
      CON02.append( 0.0 )
      CON04.append( 0.0 )

Write_file("CON01.alya", CON01)
print "CON01"
Write_file("CON02.alya", CON02)
print "CON02"
Write_file("CON03.alya", CON03)
print "CON03"
Write_file("CON04.alya", CON04)
print "CON04"
#===============================================================================================# 
print "OK!! \n\n"
