#!/usr/bin/python
import os
import re
import glob
import sys

INIT_KEY  = "COORDINATES"
END_KEY   = "END_" + INIT_KEY

#===============================================================================================# 
def Read_alya_geo(fname):
    fname += "/*.dom.geo"
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


def Read_alya_ensi_scalar(fname):
    fname = glob.glob(fname)
    
    if(len(fname)>=1):
      fname = fname[-1]
      print "|_File in: \'%s\'" % fname
    else:
      print "Error: \'%s\'!! \n" % fname 
      exit(1)

    data = open(fname, "r")
    lines = data.readlines()
    data.close()
    nline = len(lines)

    XYZ = []
    nsharps = 0 
    for i in range(nline):
      line = lines[i]
      line = line[:-1]
      if(line.find("#")<0): 
        XYZ.append(line) 

    print "  |_No elements:", len(XYZ)
    return XYZ


def Write_file(fname, data, dime=-1):
  ndata = len(data)
  fdata = open(fname, "w") 
  for i in range(ndata/dime): 
    line = data[i] 
    print>> fdata, i+1,  
    for j in range(dime): 
      xyz = data[i*dime+j]
      try:
        float(xyz)
      except ValueError:
        print "    ERROR. Line: \'%s\' \n" % (line)
        exit(1)
      
      print>> fdata, eval(xyz),  
    print>> fdata

def Write_file3D(fname, data, dime=3):
  ndata = len(data)
  fdata = open(fname, "w") 
  for i in range(ndata/3): 
    line = data[i]
    print>> fdata, i+1,  
    nvecs = ndata/3
    print>> fdata, eval(data[i]), eval(data[i + nvecs*1]), 
    if(dime==3): print>> fdata, eval(data[i + nvecs*2]), 
    print>> fdata


#===============================================================================================#

FILE_DIM  = [1, 1]
FILE_TYPE = ["TEMPE", "PRESS"]

ENSI_CASE = sys.argv[1]
ENSI_NUMS = sys.argv[2]

files_out = [] 
for file_dim, file_type in zip(FILE_DIM, FILE_TYPE):
  filein = ENSI_CASE +".ensi.*" + file_type +"*"+ ENSI_NUMS
  Field = Read_alya_ensi_scalar( filein ) 

  fileout = ENSI_CASE +"_field_" + file_type +"_"+ str(ENSI_NUMS).zfill(6) +".alya"
  Write_file(fileout, Field, file_dim) 
  files_out.append(fileout)

print "|_File out:\'%s\'" % files_out 


## 2D correction 
filein  = ENSI_CASE +".ensi.*VELOC*"+ ENSI_NUMS
fileout = ENSI_CASE +"_field_" + "VELOC" +"_"+ str(ENSI_NUMS).zfill(6) +".alya"
Field = Read_alya_ensi_scalar( filein ) 
Write_file3D(fileout, Field, dime=3)

## DENSITY correction
filein = ENSI_CASE +"_field_" + "PRESS" +"_"+ str(ENSI_NUMS).zfill(6) +".alya" 
PRESS = Read_alya_ensi_scalar( filein ) 

filein = ENSI_CASE +"_field_" + "TEMPE" +"_"+ str(ENSI_NUMS).zfill(6) +".alya"
TEMPE = Read_alya_ensi_scalar( filein )  

Rgas = 287.0  
RHO  = [] 
for f01, f02 in zip(PRESS, TEMPE):
  P = f01.split()[1]
  T = f02.split()[1]
  rho = eval(P)/Rgas/eval(T)
  RHO.append("%f" % rho)

fileout = ENSI_CASE +"_field_" + "DENSI" +"_"+ str(ENSI_NUMS).zfill(6) +".alya"
Write_file(fileout, RHO, 1)


#===============================================================================================#
print "OK!! \n\n"
#
#GEO_FILEIN = sys.argv[1]
#Rs = Read_alya_geo(GEO_FILEIN) 
#
