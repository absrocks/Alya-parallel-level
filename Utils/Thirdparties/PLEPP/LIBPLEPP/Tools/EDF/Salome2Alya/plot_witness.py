import numpy as np 
import sys
import os 
import glob 

def Gnuplot(files, col=5, ylabel=""):
  GNUPLOT = """ # created by J. MIGUEL ZAVALA-AKE
  set grid 
  set xlabel "t"
  set ylabel \""""+ylabel+"""\"
  set macro 
  COL = %d
  #set multiplot layout %d,1 rowsfirst
  set multiplot layout 1,1 rowsfirst
  %s 
  unset multiplot
  pause -1 
  """
  nfiles = len(files)
  data   ="plot "
  for i in range(nfiles-1): data += "\"%s\" u 1:COL w l,\\\n"%files[i]   
  data += "\"%s\" u 1:COL w l\n"%files[-1] 

  fname = open("plot_witness.gp", "w")
  print >> fname, GNUPLOT % (col, nfiles, data) 
  fname.close() 

  os.system("gnuplot %s"%fname.name)

if(len(sys.argv)>1): 
  fname = sys.argv[1] 
  print "|_got: \"%s\""%fname 
else:
  fname = glob.glob("*.nsa.wit")[0]
  print "|_found: \"%s\""%fname 



data = open(fname, "r")
lines = data.readlines() 
data.close() 

Times   = [] 
Witness = {} 
for line in lines:
  line = line.split() 
  if(line[0]=="#" and len(line)>1):
    if(line[1]=="Time"): Times.append(eval(line[3])) 
  elif(not line[0]=="#"):  
    Witness[line[0]] = []
Times = np.array(Times) 

for line in lines:
  line = line.split() 
  if(not line[0]=="#"): 
    Witness[line[0]].append([eval(val) for val in line[1:]]) 

fnames = [] 
for key,vals in Witness.items():
  vals = np.array(vals)
#  vals[:,:-1] /= 343.0 
#  vals[:,-1]  /= 101000.0 
  vals = np.column_stack((Times, vals)) 
  fname = "plot_witness%s.dat"%key 
  fnames.append(fname)
  np.savetxt(fname, vals)   
  print "  |_\'%s\'"%fname 


Gnuplot(fnames, 5, "P/P0") # veloz 
Gnuplot(fnames, 4, "Vz/a") # press 

print "|_OK!\n"

