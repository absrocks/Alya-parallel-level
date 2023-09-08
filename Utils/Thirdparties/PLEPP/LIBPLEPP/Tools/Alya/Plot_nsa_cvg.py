#!/usr/bin/python
import sys
import os 
import glob 

if(len(sys.argv)>1):
  fname = sys.argv[1]
  print "|_got: \"%s\""%fname
else:
  fname = glob.glob("*.nsa.cvg")
  if(len(fname)<1):
    print "ERROR: there no exist *.nsa.cvg files! \n\n"
    sys.exit(1)

  if(len(fname)>1):
    print "Choose one file *.nsa.cvg:"
    for filei in fname: print "   \'%s\'" % filei
    print "\n"
    sys.exit(1)

  fname = fname[0]
  print "|_found: \"%s\""%fname
  
#----------------------------------------------------------------------------# 
GNUPLOT = """
set log y 
set macro 
FILE_NAME = \"%s\"
set multiplot layout 2,2 rowsfirst
plot FILE_NAME u 4:5 t "Total Res" w lines ls 1 
plot FILE_NAME u 4:6 t "Mom. Res" w lines ls 1 
plot FILE_NAME u 4:7 t "Cont. Res" w lines ls 1 
plot FILE_NAME u 1:8 t "Ener. Res" w lines ls 1 
unset multiplot

pause -1 
"""

GNUPLOT = GNUPLOT % (fname)

FILE = open("nsa_props.gp", "w")
FILE.write(GNUPLOT)
FILE.close() 

os.system("gnuplot nsa_props.gp")

