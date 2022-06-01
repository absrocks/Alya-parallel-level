#-----------------------------------------------------------------------||---#
import glob

def read_file(fname, comment='', find=''):
  files = glob.glob(fname)
  if(len(files)==0):
    print "\n\n   ERROR: there is not \'%s\' \n\n" % (fname)
    sys.exit(1)

  fin = open(fname, "r")
  lines = fin.readlines()
  fin.close()

  nline = len(lines)
  Lines = []
  for i in range(nline):
    line = lines[i]
    line = ' '.join(line.split())
    if(comment != ''): line = line.split(comment)[0]
    #line = line.split()
    if(find != ''):
      if(line.find(find) != -1): Lines.append(line)
    else:
      Lines.append(line)

  return Lines
#-----------------------------------------------------------------------||---#
#-----------------------------------------------------------------------||---#
def Change_file(Fname="", Keys=[], Comment=''):
  if(len(glob.glob(Fname+"_orig")) == 0):
     os.system("cp %s %s" % (Fname, Fname+"_orig") )
  else:
     os.system("cp %s %s" % (Fname+"_orig", Fname) )

  Lines = read_file(Fname, comment=Comment)

  f01 = open(Fname+"_old", "w")
  f02 = open(Fname, "w")

  n_lines = len(Lines)
  for i in range(n_lines):
    line = Lines[i]
    print>> f01, line
    for key in Keys:
      if(line.find(key[0]) != -1):
        if(line.find(key[1]) != -1):
          line = key[-1]
    print>> f02, line

  f01.close()
  f02.close()
#-----------------------------------------------------------------------||---#


#-----------------------------------------------------------------------||---#
import sys
from optparse import OptionParser
#usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser() #(usage=usage)
parser.add_option("-C", "--case",     default="",  dest="ALYA_CASE")
parser.add_option("-N", "--domains",  default=-1,  dest="ALYA_NDOM", type=float)

(options, args) = parser.parse_args()

if(options.ALYA_CASE=="" or options.ALYA_NDOM==-1): 
  parser.print_help()
  print ""
  sys.exit()

ALYA_CASE = options.ALYA_CASE
ALYA_NDOM = options.ALYA_NDOM-1
#-----------------------------------------------------------------------||---#

#-----------------------------------------------------------------------||---#
import glob
import shutil 
import os
import math  

DIRS   = glob.glob("*PAR*")
n_DIRS = len(DIRS)

for dir in DIRS:
  shutil.rmtree(dir)

n_DIRS = 1+int(math.floor(ALYA_NDOM/100)) 
for i in range(n_DIRS+1):
  Path = "PAR%s" % (str(i*100).zfill(6)) 
  os.mkdir(Path)
#-----------------------------------------------------------------------||---#



#-----------------------------------------------------------------------||---#
DO_PARTITION = """$$ Create by J. MIGUEL ZAVALA-AKE 
$
RUN_DATA
  ALYA: %s 
END_RUN_DATA
$
PROBLEM_DATA
$
$$ 1) Run Alya.x xxx.dat 
$$ 2) Copy xxx.dat_orig > xxx.dat and 
$$ 3) Use this lines into xxx.dat in order to be able to use the partitions created  
$$    you must comment the N_DOMS-1 line and uncomment the TO_RUN 
$$ 4) mpirun -np SUBDOMAIN+1 Alya.x xxx.dat
$
  PARALL_SERVICE: On
\tPARTITION_TYPE:       FACES
\tFILE_HIERARCHY:       ON
\tFILE_OPEN_CLOSE:      Yes
\tVIRTUAL_FILE:         On, MAXIMUM_MEMORY=0.5
\t$TASK:                READ_PREPROCESS, BINARY                 $ TO RUN
\tTASK:                 ONLY_PREPROCESS, BINARY, SUBDOMAIN = %d $ N_DOMs-1
  END_PARALL_SERVICE
$
END_PROBLEM_DATA
"""
Fname = ALYA_CASE+".dat"
Change_file(Fname, [])

F01 = open(Fname, "w")
print>> F01, DO_PARTITION % (Fname, ALYA_NDOM) 
F01.close()
#-----------------------------------------------------------------------||---#


print "\n--| Loot at \'%s\' \n" % Fname  
print "OK!!"
#-----------------------------------------------------------------------||---#
#-----------------------------------------------------------------------||---#
