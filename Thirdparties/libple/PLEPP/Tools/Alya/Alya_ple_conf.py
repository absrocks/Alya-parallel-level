#!/usr/bin/python
#-----------------------------------------------------------------------||---#

#-----------------------------------------------------------------------||---#
ALYA_BIN    = "/Executables/unix/"
METHIS_PATH = "/Thirdparties/metis-5.0.2_i8/"
PLEPP_PATH  = "/Thirdparties/libple/PLEPP/"
PLE_PATH    = "/lib/"

ALYA_MODULES = ["parall", "nastin"]

import os
ALYA_PATH = os.getcwd()
ALYA_BIN  = ALYA_PATH + ALYA_BIN
#-----------------------------------------------------------------------||---#


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
from optparse import OptionParser
import sys 
#usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser() #(usage=usage)
parser.add_option("-S", "--code_saturne", default="", dest="SATURNE_PATH")
parser.add_option("-P", "--ple", default="", dest="PLE_PATH")
parser.add_option("-C", "--configure", default="", dest="ALYA_CONFIGURE")
parser.add_option("-M", "--modules", default=ALYA_MODULES, dest="ALYA_MODULES", action="append")

(options, args) = parser.parse_args()
if(options.PLE_PATH=="" or options.ALYA_CONFIGURE==""):
  parser.print_help()
  print ""
  sys.exit()

#-----------------------------------------------------------------------||---#

#=======================================================================||===#
#================================================================| COMMs |===#

#-----------------------------------------------------------------------||---#
LIBS = {}
LIBS["PLEPP"]  = ALYA_PATH + PLEPP_PATH  + "libcommdom.a"
LIBS["METHIS"] = ALYA_PATH + METHIS_PATH + "/build/Linux-x86_64/libmetis/"+"libmetis.a"
LIBS["PLE"]    = options.PLE_PATH + PLE_PATH + "libple.a"

if(len(glob.glob(LIBS["PLE"])) == 0):
   print "\n\n   ERROR: there is not \'%s\'" % ("libple.a")
   print "         \'%s\'" % (LIBS["PLE"])
   print "\n\n"
   sys.exit(1)
#-----------------------------------------------------------------------||---#


#-----------------------------------------------------------------------||---#
fname = options.ALYA_CONFIGURE 
LINES = read_file(fname) 

COMPILERS = {}
COMPILERS["f90"]    = " -DV5METIS"
COMPILERS["fpp90"]  = " -DV5METIS"
COMPILERS["fomp90"] = " -DV5METIS"
COMPILERS["cpp"]    = ""
COMPILERS["link"]   = ""
COMPILERS["libs"]   = ""
COMPILERS["f77"]    = ""
COMPILERS["fpp77"]  = ""
COMPILERS["fa2p"]   = ""

for line in LINES:
  line = line.split("#")[0]
  if(line.find("==") != -1): 
    split   = line.split("==")
    n_split = len(split)
    if(n_split>0): 
      key = split[0] 
      val = ""
      if(n_split==2): val = split[1] 
      COMPILERS[key] = val + COMPILERS[key] 

COMPILERS["libs"] = ""
for key, val in LIBS.items(): 
  COMPILERS["libs"] += " "+ val 

f01 = open(ALYA_BIN + "x_configure.txt", "w")
for key, val in COMPILERS.items():
  print>> f01, key +"== "+ val   
f01.close()
#-----------------------------------------------------------------------||---#

#-----------------------------------------------------------------------||---#
fname = ALYA_PATH + METHIS_PATH + "/include/metis.h"
Lines = read_file(fname)
key = ["#define", "IDXTYPEWIDTH"] 

f01 = open(fname+"_old", "w")
f02 = open(fname, "w")

n_lines = len(Lines)
for i in range(n_lines):
  line = Lines[i]
  print>> f01, line 
  if(line.find(key[0]) != -1): 
    if(line.find(key[1]) != -1): 
      line = line.replace("64", "32") 
  print>> f02, line

f01.close()
f02.close()
#-----------------------------------------------------------------------||---#

#-----------------------------------------------------------------------||---#
fname = ALYA_PATH + PLEPP_PATH + "/Makefile.in"
Lines = read_file(fname)
Keys  = [ ["ROOT_PLE",  "=", options.PLE_PATH], 
          ["HOME",      "=", ALYA_PATH + PLEPP_PATH], 
          ["FCOMPILER", "=", COMPILERS["link"]]
        ]

f01 = open(fname+"_old", "w")
f02 = open(fname, "w")

n_lines = len(Lines)
for i in range(n_lines):
  line = Lines[i]
  print>> f01, line
  line = line.split("#")[0]
  for key in Keys:
    if(line.find(key[0]) != -1):
      if(line.find(key[1]) != -1):
        if(line.find("$") == -1): 
          line = ' '.join(key) 
  print>> f02, line

f01.close()
f02.close()
#-----------------------------------------------------------------------||---#


#-----------------------------------------------------------------------||---#
MAKEFILE ="""#J. MIGUEL ZAVALA AKE  
ALYA_ROOT   = %s
EXECS       = $(ALYA_ROOT)/Execs/unix
METIS       = $(ALYA_ROOT)/%s
PLEPP       = $(ALYA_ROOT)/%s

all: 
\t@echo
\t@echo '\t\tmake [alya|metis5|plepp|clean]' 
\t@echo


alya:
\t@$(MAKE) -C $(EXECS)
\t#@mv $(EXECS)/Alya.x $(EXECS)/Execs 


metis5: 
\t@$(MAKE) -C $(METIS) distclean
\t@$(MAKE) -C $(METIS) config 
\t@$(MAKE) -C $(METIS) 


plepp:
\t@$(MAKE) -C $(PLEPP) -f Makefile.all fortran 


clean:
\t@$(MAKE) -C $(EXECS) clean
"""
f01 = open("Makefile", "w")
f01.write(MAKEFILE % (ALYA_PATH, METHIS_PATH, PLEPP_PATH) )
f01.close()
#-----------------------------------------------------------------------||---#


#-----------------------------------------------------------------------||---#
COMMAND  = ""
COMMAND += "./configure "
COMMAND += "-lib -x "
COMMAND += "-f=%s " % "x_configure.txt"
COMMAND += ' '.join( options.ALYA_MODULES )
COMMAND += ""

os.chdir(ALYA_BIN)
os.system(COMMAND) 
#-----------------------------------------------------------------------||---#

#-----------------------------------------------------------------------||---#
print "|_[Alya]  Directory: \'%s\' " % ALYA_PATH, 
print "ok!!"

os.chdir(ALYA_PATH)
os.system("make")
#-----------------------------------------------------------------------||---#
