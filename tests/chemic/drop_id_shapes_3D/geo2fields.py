## Usage:
## > modul load python/2.7.13 
## > python geo2fields.py <name>.geo.dat
##

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

#CON01 = []
#CON02 = []
CON03 = []
#CON04 = []
#VELOC = []
#DENSI = []
#ENTHA = []

### Circles
c_1_1 = -0.4
c_1_2 = 0.4
c_1_3 = 0.0
R_1 = 0.3

### Circles smooth boundary
cs_1_1 = 0.5
cs_1_2 = -0.75
cs_1_3 = -0.7
Rs_1 = 0.15
eps_1 = 0.5* (1.0 / 200.0)**0.9

### Rectangles
x_1_min = 0.1
x_1_max = 0.5
y_1_min = 0.1
y_1_max = 0.5
z_1_min = 0.1
z_1_max = 0.5

x_2_min = -0.9
x_2_max = -0.3
y_2_min = -0.3
y_2_max = -0.1
z_2_min = 0.1
z_2_max = 0.3

x_3_min = 0.1
x_3_max = 0.9
y_3_min = 0.7
y_3_max = 0.9
z_3_min = -0.9
z_3_max = -0.7

x_4_min = -1.0
x_4_max = -0.6
y_4_min = -0.9
y_4_max = -0.7
z_4_min = -1.0
z_4_max = -0.6

x_5_min = -0.6
x_5_max = 0.0
y_5_min = -0.7
y_5_max = -0.4
z_5_min = -0.6
z_5_max = -0.3

x_6_min = 0.1
x_6_max = 0.4
y_6_min = -0.4
y_6_max = -0.1
z_6_min = 0.1
z_6_max = 0.4

x_7_min = 0.35
x_7_max = 0.65
y_7_min = -0.45
y_7_max = -0.4
z_7_min = 0.4
z_7_max = 0.45

for pt in pts:
	x = pt[0]
	y = pt[1]
        z = pt[2]
        phi = 0.0
        
        ### Circles smooth boundary
        rs_1 = math.sqrt( (x - cs_1_1)**2 + (y - cs_1_2)**2 + (z - cs_1_3)**2 )
        if ( rs_1 > Rs_1 ):
                d_1 = -(rs_1-Rs_1)
        else:
                d_1 = Rs_1-rs_1
        phi = 0.5 * ( math.tanh( 0.5 * d_1 / eps_1 ) + 1.0 )
        
        ### Circles
	r_1 = ((x-c_1_1)**2+(y-c_1_2)**2+(z-c_1_3)**2)**0.5
        if r_1 <= R_1:
                phi = 1.0

        ### Rectangles
        if x>=x_1_min and x<=x_1_max and y>=y_1_min and y<=y_1_max and z>=z_1_min and z<=z_1_max :
                phi = 1.0

        if x>=x_2_min and x<=x_2_max and y>=y_2_min and y<=y_2_max and z>=z_2_min and z<=z_2_max :
                phi = 1.0

        if x>=x_3_min and x<=x_3_max and y>=y_3_min and y<=y_3_max and z>=z_3_min and z<=z_3_max :
                phi = 1.0

        if x>=x_4_min and x<=x_4_max and y>=y_4_min and y<=y_4_max and z>=z_4_min and z<=z_4_max :
                phi = 1.0

        if x>=x_5_min and x<=x_5_max and y>=y_5_min and y<=y_5_max and z>=z_5_min and z<=z_5_max :
                phi = 1.0

        if x>=x_6_min and x<=x_6_max and y>=y_6_min and y<=y_6_max and z>=z_6_min and z<=z_6_max :
                phi = 1.0

        if x>=x_7_min and x<=x_7_max and y>=y_7_min and y<=y_7_max and z>=z_7_min and z<=z_7_max :
                phi = 1.0



	#VELOC.append( [U,V] )
	#CON01.append( Yc )
	#CON02.append( 0.0 )
	CON03.append( phi )
	#CON04.append( 0.0 )
	#ENTHA.append( H )
#
## WRITE OUTPUT FILESnew
#                       
#Write_file("VELOC.alya", VELOC, 2)
#print "VELOC"        
#Write_file("CON01.alya", CON01)
#print "CON01"        
#Write_file("CON02.alya", CON02)
#print "CON02"        
Write_file("CON03.alya", CON03)
print "CON03"        
#Write_file("CON04.alya", CON04)
#print "CON04"        
#Write_file("ENTHA.alya", ENTHA)
#print "ENTHA"        
                     
#===============================================================================================# 
print "OK!! \n\n"    
