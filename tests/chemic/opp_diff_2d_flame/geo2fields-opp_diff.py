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

CON01 = []
CON02 = []
CON03 = []
CON04 = []
VELOC = []
DENSI = []
ENTHA = []
TEMPE = []

Zmax = 1.0
Zmin = 0.0
t    = 0.001
tvel = 0.001

uIn = 0.1

for pt in pts:
	y = pt[1]
	x = pt[0]
	yrel = (t/2.0 - y) / t
	yrelVel =min(1.0,max(-1.0,  (2.0*y) / tvel))
	xrelVel =min(1.0,max(-1.0,  (2.0*x) / tvel))

	yrelCent = max(0.0, 1.0-2.0*abs(yrel-0.5))

	Z = 0.18
	Yc = 1.0
	#Z = min(Zmax,max(Zmin, Zmin + yrel * (Zmax-Zmin) ))
	h = -112.32671356 + (-951485.125+112.32671356)*Z/Zmax
	#Yc = min(0.3,max(0.0, yrelCent * 0.3 ))

	U = uIn*xrelVel
	V = -1*uIn*yrelVel
	#U = 0.0 
	#V = 0.0 

	VELOC.append( [U,V] )
	CON01.append( Yc )
	CON02.append( 0.0 )
	CON03.append( Z )
	CON04.append( 0.0 )
	ENTHA.append( h )
	TEMPE.append( 298.0 )
#
## WRITE OUTPUT FILESnew
#                       
Write_file("VELOC.alya", VELOC, 2)
print "VELOC"        
Write_file("CON01.alya", CON01)
print "CON01"        
Write_file("CON02.alya", CON02)
print "CON02"        
Write_file("CON03.alya", CON03)
print "CON03"        
Write_file("CON04.alya", CON04)
print "CON04"        
Write_file("ENTHA.alya", ENTHA)
print "ENTHA"        
Write_file("TEMPE.alya", TEMPE)
print "TEMPE"        
                     
#===============================================================================================# 
print "OK!! \n\n"    
