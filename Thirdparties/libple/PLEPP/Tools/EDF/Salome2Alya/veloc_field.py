import os 
import sys
import glob 
import math
import numpy as np

os.system('rm fluid_velocity_field')
veloc=open('coudes_nsi.ensi.VELOC-000048','r')
field=open('fluid_velocity_field', 'a')

lines = veloc.readlines()
nblines = (len(lines)) / 3
print nblines
for i in range(nblines):
  x = lines[i].strip()
  y = lines[nblines+i].strip()
  z = lines[2*nblines+i].strip()
  print >> field, i+1,x,y,z
