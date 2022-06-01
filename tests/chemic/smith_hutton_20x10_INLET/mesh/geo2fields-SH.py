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

for pt in pts:
      u =  2.0 * pt[1] * (1 - pt[0]*pt[0] )
      v = -2.0 * pt[0] * (1 - pt[1]*pt[1] )
      VELOC.append( [u,v] )

Write_file("VELOC.alya", VELOC, 2)
print "VELOC" 
#===============================================================================================# 
print "OK!! \n\n"
