#!/usr/bin/python
import csv
import sys

BASE_NAME = "csv2gp"

#------------------------------------------------------------------------------#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-Y', action='store', dest='Ys',
                    default=[], type=str, nargs='+',
                    help='DATA',
                    )

parser.add_argument('-X', action='store', dest='Xs',
                    default=[], type=str, nargs='+',
                    help='DATA',
                    )

parser.add_argument('-I', action='store', dest='Fin01',
                    default='*', type=str, nargs=1,
                    help='Fin',
                    )

parser.add_argument('-O', action='store', dest='Fout01',
                    default=[BASE_NAME], type=str, nargs=1,
                    help='Fout',
                    )

Results = parser.parse_args()
Xs     = Results.Xs[:]
Ys     = Results.Ys[:]
Fin01  = Results.Fin01[0]
Fout01 = Results.Fout01[0]

n_Xs  = len(Xs)
n_Ys  = len(Ys)

if((n_Xs==0)or(n_Ys==0)or(Fin01=='*')): 
  parser.print_help()
  print "EXIT!!\n"
  sys.exit(1)

BASE_NAME = Fout01 

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
import numpy as np

f           = open(Fin01, "r")
reader_dics = csv.DictReader(f)
headers     = reader_dics.fieldnames
n_headers   = len(headers)

PROPS = Xs+Ys 
for key in PROPS: 
  if not key in headers:
    print "ERROR: '%s' " % (key)
    print " KEYS: ", headers  
    sys.exit() 

COLS = { headers[i]:i for   i in range(n_headers) } 
IDx  = [ COLS[key]    for key in Xs ] 
IDy  = [ COLS[key]    for key in Ys ]
IDS  = IDx + IDy 
print "IDS:", IDS

ROWs   = [] 
for row_dic in reader_dics:
  rows = row_dic.values() 
  vals = [ eval(row_dic[idx]) for idx in PROPS ] 
  ROWs.append(vals)

n_ROWs = len(ROWs)
print "n_ROWs:", n_ROWs 

n_PROPS = len(PROPS)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
ROWs = np.array( ROWs )
np.savetxt("%s.dat"%BASE_NAME, ROWs, header=" ".join(PROPS) )

#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
Fout01 = open("%s.gp"%BASE_NAME, "w")
print>> Fout01, "#set terminal pngcairo size 350,262 enhanced font 'Verdana,10' "
print>> Fout01, "#set output '%s.png' " % BASE_NAME 
print>> Fout01, "#set terminal postscript eps enhanced color font 'Helvetica,10' "
print>> Fout01, "#set output '%s.eps' " % BASE_NAME
print>> Fout01, "#set xlabel 'Z [m]' "
print>> Fout01, "#set ylabel '{/Symbol s} [-]' "
print>> Fout01, "#set key left top"
print>> Fout01, "set grid"
print>> Fout01, "set macros"
print>> Fout01, "set ytics nomirror"
print>> Fout01, "#set log y"
print>> Fout01, "#set y2tics"
print>> Fout01 
for i in range(n_PROPS):
  key   = PROPS[i]
  key   = key.replace(":", '')
  txt01 = "%s = '%d' " % (key,i+1)
  print>> Fout01, txt01
print>> Fout01

YCOLS  = Ys 
XCOLS  = Xs * n_Ys
FNAMES = ["%s.dat"%BASE_NAME] * n_Ys
print>> Fout01, "plot \\"

for i in range(n_Ys):
  f = FNAMES[i]
  x = XCOLS[i]
  y = YCOLS[i]
  x = x.replace(":", '')
  y = y.replace(":", '')

  print>> Fout01, "'%s' " % ( f ),
  print>> Fout01, "ev 1 u @%s:@%s t '%d' w lp axes x1y1 lt 1 pt 1 " % (x, y, i+1),
  if(i<n_Ys-1): print>> Fout01, ",\\"
  else:         print>> Fout01, ""

print>> Fout01
print>> Fout01, "pause -1"
print>> Fout01, "# pt  6 circle"
print>> Fout01, "# pt 10 triangle down"
print>> Fout01, "# pt 12 diamond"
print>> Fout01, "\n\n# The big boss J. Miguel Zavala Ake!!"

Fout01.close()

print "OK!! -> \ngnuplot %s.gp\n" % (BASE_NAME)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
