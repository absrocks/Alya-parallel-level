#=========================================================================||===#
#=========================================================================||===#
from optparse import OptionParser

import glob
PATH_FILE = "Alya_runner.txt"
flist = glob.glob(PATH_FILE) 
if( len(flist)!=1 ): 
  print "ERROR:", flist, "\n" 
  exit(1)
else: 
  print "FIND: \'%s\' \n" % ( flist[0] ) 

fin = open(PATH_FILE, "r")
lines = fin.readlines()
fin.close() 


TOKENS = ["ALYA_PATH", "ALYA_BIN"] 
DICT = {} 
for i in range( len(lines) ):
  line = lines[i].rstrip()
  line = line.replace("  ", "") 
  line = line.replace( " ", "")
  for token in TOKENS: 
    if( line.find(token) != -1 ): 
      ##print "\'%s\'" % line.split("=")[1] 
      DICT[ token ] = line.split("=")[1] 

 
ALYA_PATH = DICT["ALYA_PATH"] 
ALYA_BIN  = DICT["ALYA_BIN"] 

#=========================================================================||===#
#=========================================================================||===#
usage = "\t%prog [options] arg1 arg2!!\n"
parser = OptionParser(usage=usage)
parser.add_option("-n", "--nmpis", type="int", nargs=2, 
                  dest="nmpis", default=(0,0), help="mpi size") 
parser.add_option("-c", "--cases", type="string", nargs=2, 
                  dest="cases", default=('',''), help="mpi size") 

(options, args) = parser.parse_args()

#if len(args) != 2:
#  parser.print_help()
#  parser.error("\n  !!Incorrect number of arguments!!\n  exit\n\n")


RUNNERS = []
for nmpi, case in zip(options.nmpis, options.cases): 
  RUNNER = ""
  RUNNER += "-n %d " % nmpi
  RUNNER += "".join( [ALYA_PATH, ALYA_BIN] )
  RUNNER += " " + case
  RUNNER += ""
  RUNNERS.append(RUNNER) 

RUNNERS = "mpirun " + " : ".join( RUNNERS )
fout = open("runner.sh", "w")
print >> fout, RUNNERS 
fout.close() 
##print RUNNERS

#print args
#print options


import os 
os.system(RUNNERS) 

#=========================================================================||===#
#=========================================================================||===#
print "OK!!\n"

#
#if len(args) <= 0:
#  parser.error("Incorrect number of arguments!!\n")
#if options.verbose:
#  print "reading %s..." % options.filename
#
