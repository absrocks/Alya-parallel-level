#!/usr/bin/python

import numpy as np 
import sys
import os 
import glob 

from optparse import OptionParser
usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-P", "--Press", default=101325.0, dest="press") 
parser.add_option("-F", "--File", default="", dest="fname")
parser.add_option("-W", "--Wits", default="", dest="wits")


(options, args) = parser.parse_args()

Pref = options.press 

#--------------------------------------------------------------------------------------------# 
def KerDat(): 
  fname = glob.glob("*.ker.dat")
  if(len(fname)<1):
    print "WARNING: there no exist *.nsa.wit files! \n\n"
    pass 
  else: 
    fname = fname[0] 
    print "|_\'%s\'" % fname 

    data = open(fname, "r")
    lines = data.readlines()
    data.close()

    ok = False  
    for line in lines:  
      if(line.find("WITNESS_POINTS")>0):     ok = True 
      if(line.find("END_WITNESS_POINTS")>0): ok = False 
      if(ok): print line[:-1]


def Gnuplot(files, col=5, ylabel=""):
  GNUPLOT = """ # created by J. MIGUEL ZAVALA-AKE
  set grid 
  set xlabel "t"
  set ylabel \""""+ylabel+"""\"
  set macro 
  Pref = """+ str(Pref) +"""\n 
  COL = %d
  #set multiplot layout %d,1 rowsfirst
  set multiplot layout 1,1 rowsfirst
  %s 
  unset multiplot
  pause -1 
  """
  nfiles = len(files)
  data   ="plot \\\n"
  if(ylabel=="P"): data += " Pref t \"%s\" w l,\\\n" % str(Pref)  
  for i in range(nfiles-1): data += "\"%s\" u 1:COL w l,\\\n"% (files[i])
  data += "\"%s\" u 1:COL w l\n"%files[-1] 

  fname = open("plot_witness.gp", "w")
  print >> fname, GNUPLOT % (col, nfiles, data) 
  fname.close() 

  os.system("gnuplot %s"%fname.name)

#if(len(sys.argv)>1): 
#fname = sys.argv[1]
if( options.fname!="" ): 
  fname = options.fname
  print "|_got: \"%s\""%fname 
else:
  fname = glob.glob("*.nsa.wit")
  if(len(fname)<1): 
    parser.print_help()
    #print "ERROR: there no exist *.nsa.wit files! \n\n" 
    print "\n\n"
    sys.exit(1)

  if(len(fname)>1): 
    print "Choose one file *.nsa.wit:"  
    for filei in fname: print "   \'%s\'" % filei    
    print "\n"
    sys.exit(1)

  fname = fname[0]
  print "|_found: \"%s\""%fname 
#--------------------------------------------------------------------------------------------# 

fbase, fextension = os.path.splitext(fname) 
fbase = os.path.splitext(fbase)[0]

data = open(fname, "r")
lines = data.readlines() 
data.close() 

Times   = [] 
Witness = {}
Keys    = {} 
for line in lines:
  if(line.find("#")>=0): 
    if(line.find("Column")>=0): 
      tokens = line.split() 
      Keys[ tokens[1] ] = tokens[-1]

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
  #vals[:,:-1] /= 343.0 
  #vals[:,-1]  /= 101000.0 
  vals = np.column_stack((Times, vals)) 
  fname = fbase+"_witness%s.dat"%key 
  fnames.append(fname)
  np.savetxt(fname, vals)   
  print "  |_\'%s\'"%fname 

KerDat() 


for key,vals in Keys.items():
  print "  |_\'%s\': %s" % (key,vals)


if(options.wits!=""):
  iname = options.wits
  idx   = eval(Keys.get(iname))
  Gnuplot(fnames, idx, iname)  # press 

else: 
  for key,vals in Keys.items():
    if(key!="ISET"): 
      Gnuplot(fnames, eval(vals), key)  # press 


#Gnuplot(fnames, 5, "P")  # press 
#Gnuplot(fnames, 4, "Vz") # veloz 


#--------------------------------------------------------------------------------------------# 
#if(len(sys.argv)==3): 
#if(False):
#  f02 = sys.argv[2]   
#  f02, fextension = os.path.splitext(f02) 
#  f02 = os.path.splitext(f02)[0]
#  print "  |_Comparing: \'%s\'.nsa.wit \'%s\'.nsa.wit"% (fbase, f02)
#
#  fnames = [] 
#  for key,vals in Witness.items():
#    fname = f02+"_witness%s.dat"%key 
#    fnames.append(fname)
#
#    fname = fbase+"_witness%s.dat"%key 
#    fnames.append(fname)
#
#  Gnuplot(fnames, 4, "Vz/a") # press 
#  Gnuplot(fnames, 5, "P/P0") # veloz 

print "|_OK!\n"
#--------------------------------------------------------------------------------------------#
