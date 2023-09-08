#!/usr/bin/python
import sys
import os 
import glob 

if(len(sys.argv)>1):
  fname = sys.argv[1]
  print "|_got: \"%s\""%fname
else:
  fname = glob.glob("*.nsa.mxm")
  if(len(fname)<1):
    print "ERROR: there no exist *.nsa.mxm files! \n\n"
    sys.exit(1)

  if(len(fname)>1):
    print "Choose one file *.nsa.mxm:"
    for filei in fname: print "   \'%s\'" % filei
    print "\n"
    sys.exit(1)

  fname = fname[0]
  print "|_found: \"%s\""%fname
  
#----------------------------------------------------------------------------# 
GNUPLOT = """
#set log y
set grid  
set macro 
FILE_NAME = \"%s\"

Rgas = 8.314472e7  
set multiplot layout 2,2 rowsfirst
#plot FILE_NAME u 1:($6/($4*$8)/Rgas) t "Min/Rgas" w lines ls 1, "" u 1:($7/$5*$9/Rgas) t "Max/Rgas" w lines ls 2
plot FILE_NAME u 1:10 t "Min" w lines ls 1, "" u 1:11 t "Max" w lines ls 2 
plot FILE_NAME u 1:4 t "Dens Min" w lines ls 1, "" u 1:5 t "Dens Max" w lines ls 2 
plot FILE_NAME u 1:6 t "Pres Min" w lines ls 1, "" u 1:7 t "Pres Max" w lines ls 2 
plot FILE_NAME u 1:8 t "Temp min" w lines ls 1, "" u 1:9 t "Temp Max" w lines ls 2 
unset multiplot

pause -1 
"""

GNUPLOT = GNUPLOT % (fname)

FILE = open("nsa_props.gp", "w")
FILE.write(GNUPLOT)
FILE.close() 

os.system("gnuplot nsa_props.gp")

