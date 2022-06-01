#!/usr/bin/python
#-----------------------------------------------------------------------||---#

#-----------------------------------------------------------------------||---#
ALYA_BIN    = "/Executables/unix/"
METHIS_PATH = "/Thirdparties/metis-5.0.2_i8/"
PLEPP_PATH  = "/Thirdparties/libple/PLEPP/"
PLE_PATH    = "/lib/"

ALYA_MODULES    = ["parall", "nastin", "temper"]
ALYA_DIRECTIVES = ["V5METIS"] 


import os
ALYA_PATH = os.getcwd()
ALYA_BIN  = ALYA_PATH + ALYA_BIN
#-----------------------------------------------------------------------||---#

#-----------------------------------------------------------------------||---#
import sys
My = sys.argv[0]
My = os.path.splitext( os.path.basename(My) )[0]
F01 = open(My+".log", "w")
#-----------------------------------------------------------------------||---#

#-----------------------------------------------------------------------||---#
MAKEFILE ="""#J. MIGUEL ZAVALA AKE  
ALYA_ROOT   = %s
EXECS       = $(ALYA_ROOT)/Executables/unix/
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
\t@$(MAKE) -C $(PLEPP) -f Makefile.all libplepp 


clean:
\t@$(MAKE) -C $(EXECS) clean
"""
f01 = open("Makefile", "w")
f01.write(MAKEFILE % (ALYA_PATH, METHIS_PATH, PLEPP_PATH) )
f01.close()
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
n_ALYA_MODULES    = len(ALYA_MODULES)
n_ALYA_DIRECTIVES = len(ALYA_DIRECTIVES)

from optparse import OptionParser
import sys 
#usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser() #(usage=usage)
#parser.add_option("-S", "--code_saturne", default="", dest="SATURNE_PATH")
parser.add_option("-P", "--ple",        default="",              dest="PLE_PATH")
parser.add_option("-C", "--configure",  default="",              dest="ALYA_CONFIGURE")
parser.add_option("-M", "--modules",    default=ALYA_MODULES,    dest="ALYA_MODULES",    action="append")
parser.add_option("-D", "--directives", default=ALYA_DIRECTIVES, dest="ALYA_DIRECTIVES", action="append")
(options, args) = parser.parse_args()

message = ""
message+="%s " % (sys.argv[0])
message+="-P%s " % (options.PLE_PATH)
message+="-C%s " % (options.ALYA_CONFIGURE)
for i in range(n_ALYA_MODULES, len(ALYA_MODULES)): message+="-M%s " %(options.ALYA_MODULES[i])
for i in range(n_ALYA_DIRECTIVES, len(ALYA_DIRECTIVES)): message+="-D%s " %(options.ALYA_DIRECTIVES[i])
message+="\n"
print>> F01, message 

if(options.PLE_PATH=="" or options.ALYA_CONFIGURE==""):
  parser.print_help()
  print ""
  sys.exit()
#-----------------------------------------------------------------------||---#


#-----------------------------------------------------------------------||---#
LIBS = {}
LIBS["PLEPP"]  = ALYA_PATH + PLEPP_PATH      + "/Wrappers/Fortran/" + "libcommdom.a"
LIBS["METHIS"] = ALYA_PATH + METHIS_PATH     + "/build/Linux-x86_64/libmetis/"+"libmetis.a"
LIBS["PLE"]    = options.PLE_PATH + PLE_PATH + "libple.a"

if(len(glob.glob(LIBS["PLE"])) == 0):
   message = "" 
   message+= "\n\n   ERROR: there is not \'%s\' \n" % ("libple.a")
   message+= "         \'%s\'" % (LIBS["PLE"]) 
   message+= "\n\n"
   print message 
   print>> F01, message 
   sys.exit(1)
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
Fname = ALYA_PATH + "/Sources/kernel/parall/" + "par_code_split_universe.f90" 
Keys  = [ 
          ["implicit none",         "", "\nuse mod_commdom_edf\nimplicit none\n"], 
          ["PAR_UNIVERSE_SIZE > 1", "", "\ncall par_commdom_init()\nif(PAR_UNIVERSE_SIZE > 1) then\n"]
        ]
Change_file(Fname, Keys) 


Fname = ALYA_PATH + "/Sources/kernel/domain/" + "domain.f90"
Keys  = [
          ["implicit none",  "", "\nuse mod_commdom_edf\nimplicit none\n"],
          ["end subroutine", "", "\ncall par_commdom_set_mesh02()\nend subroutine domain\n"]
        ]
Change_file(Fname, Keys)


Fname = ALYA_PATH + "/Sources/modules/temper/" + "Temper.f90"
Keys  = [
          ["implicit none",  "", "\nuse mod_commdom_edf\nuse def_domain\nuse def_temper\nimplicit none\n"],
          ["ITASK_DOITER",  "", "\nif(kfl_regim_tem==0) call par_commdom_cfd_sync_apps(dtinv)\ncase(ITASK_DOITER)\n"],
          ["tem_doiter()",    "", "\nif(kfl_regim_tem==0) call par_commdom_cfd_coupling_vars( TEMPE(1:npoin,1) )\ncall tem_doiter()\n"]
        ]
Change_file(Fname, Keys)




f01 = ALYA_PATH + PLEPP_PATH + "/Tools/Alya/" + "mod_commdom_edf.f90"
f02 = ALYA_PATH + "/Sources/kernel/parall/"   + "mod_commdom_edf.f90"
os.system("ln -s %s %s" %(f01, f02) )


f01 = ALYA_PATH + PLEPP_PATH + "/Tools/Alya/" + "tem_commdom.f90"
f02 = ALYA_PATH + "/Sources/modules/temper/"  + "tem_commdom.f90"
os.system("ln -s %s %s" %(f01, f02) )

#-----------------------------------------------------------------------||---#

#-----------------------------------------------------------------------||---#
fname = options.ALYA_CONFIGURE
LINES = read_file(fname)

#DIRECTIVES = " -DV5METIS"
DIRECTIVES = " ".join(["-D"+directive for directive in options.ALYA_DIRECTIVES])


COMPILERS = {}
COMPILERS["f90"]    = " " + DIRECTIVES
COMPILERS["fpp90"]  = " " + DIRECTIVES
COMPILERS["fomp90"] = " " + DIRECTIVES
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
      #COMPILERS[key] = val + COMPILERS[key]
      COMPILERS[key] = val + COMPILERS.get(key,"")

COMPILERS["libs"]  = ""
COMPILERS["libs"] += " " + LIBS["METHIS"]
#COMPILERS["libs"] += " -lstdc++"
COMPILERS["libs"] += " " + LIBS["PLEPP"]
COMPILERS["libs"] += " " + LIBS["PLE"]
COMPILERS["libs"] += " -lstdc++ -lmpi -lmpi_cxx"


f01 = open(ALYA_BIN + "x_configure.txt", "w")
for key, val in COMPILERS.items():
  print>> f01, key +"== "+ val
f01.close()
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
COMMAND  = ""
COMMAND += "./configure "
COMMAND += "-x "
COMMAND += "-f=%s " % "x_configure.txt"
COMMAND += ' '.join( options.ALYA_MODULES )
COMMAND += ""

os.chdir(ALYA_BIN)
os.system(COMMAND) 
print>> F01, COMMAND
#-----------------------------------------------------------------------||---#

#-----------------------------------------------------------------------||---#
print "|_[Alya]  Directory: \'%s\' " % ALYA_PATH, 
print "ok!!"

os.chdir(ALYA_PATH)
os.system("make")
#-----------------------------------------------------------------------||---#

#-----------------------------------------------------------------------||---#
import time 
print>> F01, "\nDirectory: \'%s\', " % ALYA_PATH,  
print>> F01, time.strftime("%d/%m/%Y"), ",",
print>> F01, "OK!"
F01.close()
#-----------------------------------------------------------------------||---#
