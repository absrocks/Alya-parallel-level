#!/usr/bin/python
#------------------------------------------------------------------------------#
"""
make      metis4
make      libple
make      libplepp  SATURNE=1 ROOT_SATURNE=/home/bsc21/bsc21704/z2015/REPOSITORY/SATURNE333/
make -j 4 COMMDOM=1 SATURNE=1 ROOT_SATURNE=/home/bsc21/bsc21704/z2015/REPOSITORY/SATURNE333/

echo "#make libplepp  SATURNE=1 ROOT_SATURNE=/home/bsc21/bsc21704/z2015/REPOSITORY/SATURNE333/" >> COMMDOM_recompiler.sh
echo " make COMMDOM=1 SATURNE=1 ROOT_SATURNE=/home/bsc21/bsc21704/z2015/REPOSITORY/SATURNE333/" >> COMMDOM_recompiler.sh
"""

MAKE ="""
make      metis4
make      libple
make      libplepp   %s   

make -j 4 COMMDOM=%s %s  
"""
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
import argparse
parser = argparse.ArgumentParser()
#
#parser.add_argument('-X', action='append', dest='XXX',
#                    default=[], type=str, nargs="+",
#                    help='gcc|intel',
#                    )
#
parser.add_argument('-S', action='store', dest='ROOT_SATURNE',
                    default='-', type=str, nargs=1,
                    help='ROOT_SATURNE',
                    )

parser.add_argument('-C', action='store', dest='COMMDOM',
                    default='-', type=str, nargs=1,
                    help='COMMDOM=XXX',
                    )

Results      = parser.parse_args()
ROOT_SATURNE = Results.ROOT_SATURNE[0]
COMMDOM      = Results.COMMDOM[0]

#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
import os
import sys

if( os.path.isdir( ROOT_SATURNE ) ): 
  if( os.path.exists(ROOT_SATURNE+"Execs/include/ple_config.h") ): 
    pass  
  else:
    print "ERROR: dont find '%s' " % ( ROOT_SATURNE+"Execs/include/ple_config.h" )
    print "EXIT!!\n"
    sys.exit(1) 


SATURNE = ""
if( not ROOT_SATURNE=='-'): 
  SATURNE = "SATURNE=1 ROOT_SATURNE=%s" % ( ROOT_SATURNE )  

if( COMMDOM=='-'): 
   print "ERROR: no valid COMMDOM=%s' " % (  COMMDOM )
   print "EXIT!!\n"
   sys.exit(1)


MAKE = MAKE % ( SATURNE, COMMDOM, SATURNE  )  


from time import gmtime, strftime
fout = open("COMMDOM_recompiler.sh" , "w") 
print >> fout, "### ", strftime("%Y-%m-%d %H:%M:%S", gmtime())  
print >> fout, MAKE
fout.close() 
 

os.system( MAKE ) 
print "CREATED '%s' " % ("COMMDOM_recompiler.sh")
print "ok! \n"

#------------------------------------------------------------------------------#
if(True):
  pass 
else: 
  import sys
  parser.print_help()
  print "EXIT!!\n"
  sys.exit(1)
#------------------------------------------------------------------------------#
