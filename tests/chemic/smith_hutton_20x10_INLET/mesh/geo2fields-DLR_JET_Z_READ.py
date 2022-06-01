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

p   = 101325.0
R_0 = 8.3144621 

#INLET1: WMEAN, LAMBDA, MU, CP1_1, CP1_2, CP1_3, CP1_4, CP1_5, CP1_6, CP2_1, CP2_2, CP2_3, CP2_4, CP2_5, CP2_6:
# 2.88503972e-02 2.63586968e-02 1.85386681e-05 9.79931797e+02 1.39253289e-01 -3.06351853e-04 6.98498829e-07 -3.60312904e-10 -2.96820501e+05 8.64974163e+02 4.28527545e-01 -1.75298449e-04 3.56653572e-08 -2.84913268e-12 -2.75968992e+05

#INLET2: WMEAN, LAMBDA, MU, CP1_1, CP1_2, CP1_3, CP1_4, CP1_5, CP1_6, CP2_1, CP2_2, CP2_3, CP2_4, CP2_5, CP2_6:
#1.67367477e-02 6.06863104e-02 1.57818608e-05 1.68455180e+03 1.28034827e-01 1.30677022e-03 -7.44606974e-07 7.06497222e-11 -1.50305697e+06 1.20852727e+03 1.79243927e+00 -6.73259900e-04 1.27068188e-07 -9.37492971e-12 -1.40113863e+06

cp1_1 = 9.79931797e+02
cp2_1 = 1.39253289e-01
cp3_1 = -3.06351853e-04
cp4_1 = 6.98498829e-07
cp5_1 = -3.60312904e-10
cp6_1 = -2.96820501e+05

MW_1  = 2.88503972e-02
T_1   = 298.0

h_1   = (((((cp5_1/5*T_1+cp4_1/4)*T_1+cp3_1/3)*T_1+cp2_1/2)*T_1+cp1_1)*T_1+cp6_1)
rho_1     = p * MW_1 / ( R_0 * T_1)
print"AIR: T, rho, h = ",T_1,rho_1,h_1

cp1_2 = 1.68455180e+03
cp2_2 = 1.28034827e-01
cp3_2 = 1.30677022e-03
cp4_2 = -7.44606974e-07
cp5_2 = 7.06497222e-11
cp6_2 = -1.50305697e+06

MW_2  = 1.67367477e-02
T_2   = 298.0

h_2   = (((((cp5_2/5*T_2+cp4_2/4)*T_2+cp3_2/3)*T_2+cp2_2/2)*T_2+cp1_2)*T_2+cp6_2)
rho_2 = p * MW_2 / ( R_0 * T_2)

print"FUEL: T, rho, h = ",T_2,rho_2,h_2

# REACTANTS: MIXTURE-AVERAGED QUANTITIES
Z   = 0.0                     # Stoichoimetric
h   = Z*h_2 + (1.0-Z)*h_1       # Enthalpy
rho = Z*rho_2 + (1.0-Z)*rho_1   # Density

u   = 1.1                     # Inlet velocity == flame Speed
T   = 298.0                      # Temperature (given manually) 

print"REACTANT PROPERTIES"
print"Z =",Z
print"T =",T
print"rho=",rho
print"u=",u
print"h=",h

# PRODUCTS: adiabatic conditions (given manually)
T_prod    = 1800  # guess
rho_prod  =  0.167 # guess 

# Mass conservation
u_prod    = rho*u/rho_prod

x1   = 0.01
x2   = 0.65
rmax = 0.08
zmin = 7.00000000e-02
zmax = 3.50000000e-01

inpcoldsol = False
readZfile  = False

if (not inpcoldsol ):
   for pt in pts:
         CON01.append( 0.0 )
         TEMPE.append( 298.0 )
         VELOC.append( [0.0, 0.0, 0.0] )
         CON03.append( 0.0 )
         CON02.append( 0.0 )
         CON04.append( 0.0 )
         ENTHA.append( h )
   Write_file("CON01.alya", CON01)
   print "CON01"
   Write_file("CON02.alya", CON02)
   print "CON02"
   Write_file("CON03.alya", CON03)
   print "CON03"
   Write_file("CON04.alya", CON04)
   print "CON04"
#   Write_file("VELOC.alya", VELOC, 3)
#   print "VELOC"
   Write_file("ENTHA.alya", ENTHA)
   print "ENTHA"
   Write_file("TEMPE.alya", TEMPE)
   print "TEMPE"
else:
  zloc = 1.0
  zfile = open('dummy.alya','w')
  if readZfile:
    if os.path.isfile('CON03.alya'):
      zfile.close()
      zfile = open('CON03.alya','r')
    else:
      raise Exception("CON03.alya does not exist")
  else:
    zfile.close()
    zfile = open('dummy.alya','r')
  for pt in pts:
    if readZfile:
      zstr = zfile.readline()
      ind = 0
      zstr = zstr.strip()
      zstr = zstr.split()
  #    for i in range(len(zstr)):
  #      if ' ' in zstr[i]:
  #        ind = i
  #    print ind,zstr,len(zstr)
  #    zloc = float(zstr[ind+1:len(zstr)])
      zloc = float(zstr[1])

    r = math.sqrt( pt[1]*pt[1] + pt[2]*pt[2])
    if ( pt[0] > x1 and pt[0] < x2 and r < rmax and zloc > zmin  and zloc < zmax):
      CON01.append( 1.0 )
      TEMPE.append( 2200.0 )
    else:
      CON01.append( 0.0 )
      TEMPE.append( 288.0 )
  zfile.close()
  os.remove('dummy.alya')
  Write_file("CON01.alya", CON01)
  print "CON01"
  Write_file("TEMPE.alya", TEMPE)
  print "TEMPE"
print "OK!! \n\n"
#===============================================================================================# 
