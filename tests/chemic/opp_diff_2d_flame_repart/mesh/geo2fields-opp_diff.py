#!/usr/bin/python
#import numpy as np
import sys
import os
import glob
import math
import numpy as np
from scipy.special import erf

#===============================================================================================# 
INIT_KEY  = "COORDINATES"
END_KEY   = "END_" + INIT_KEY

def strainFlow(x,y,a,width=np.inf):
	if abs(x) < width/2.0:
		U = a * x
		V = -1.0 * a * y
	else:
		U = a * width/2.0 * x/abs(x)
		V = 0.0

	return(U,V)



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


ymin = np.inf 
ymax = -np.inf 
for pt in pts:
	y = pt[1]
	ymin = min(ymin,y)
	ymax = max(ymax,y)


dy = ymax - ymin
a = 4.0 # 1/s

for pt in pts:
	y = pt[1]
	x = pt[0]

	zeta = y /  (ymax/7.0)
	if abs(x) < 0.5*dy:
		Z = 0.05
		Yc=	1.0
	else:
		Z= 0.5 * (1.0 - erf(zeta)) 
		Yc=	0.0
	H_f=   -4.67195581e+06
	H_ox=  -1.02005890e+04
	H = H_ox + (H_f-H_ox)*Z

		
	U,V = strainFlow(x,y,a,width=1.5*dy)

	#U = 0.0 
	#V = 0.0 

	VELOC.append( [U,V] )
	CON01.append( Yc )
	CON02.append( 0.0 )
	CON03.append( Z )
	CON04.append( 0.0 )
	ENTHA.append( H )

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


U,Vmax = strainFlow(0.0,ymin,a,width=dy)
U,Vmin = strainFlow(0.0,ymax,a,width=dy)

print('Vmax={}, Vmin={}'.format(Vmax,Vmin))
                     
#===============================================================================================# 
print "OK!! \n\n"    
