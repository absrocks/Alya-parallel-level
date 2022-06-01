import shutil
import glob
import os
import numpy as np 
PWD = os.getcwd()

def Replacer(fout=None, cmd=None, Dic=None):
   dummy = cmd
   for key, val in Dic.iteritems():
     dummy = dummy.replace(key,val)

   fout = open(fout, "w")
   print>> fout, dummy
   fout.close()


CMD01_PLEPP="""
#BSUB -n XXX_RANK  
#BSUB -R"span[ptile=16]"
#BSUB -o out.txt  
#BSUB -e err.txt 
#BSUB -J TAU_FF_RANKxSS_RANK   
#BSUB -W 00:XXX_TIME  
# #BSUB -q bsc_case  


export TAU_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/TAU/TAU2251/EXEC2/x86_64
export TAU_MAKEFILE=$TAU_PATH/lib/Makefile.tau-icpc-mpi-pdt
export PATH=$PATH:$TAU_PATH/bin
#
#export TAU_TRACE=1 
#
#export TAU_CALLPATH=1
#export TAU_CALLPATH_DEPTH=10 
#
export TAU_SUMMARY=1  
export TAU_SYNCHRONIZE_CLOCKS=1 
export TAU_COMPENSATE=1 
#export TAU_TRACK_MESSAGE=1 # To enable MPI message statistics
 

ALYA_PATH=XXX_PATH

ALYAi=F_RANK  
ALYAj=S_RANK 

time mpirun \\
-np $ALYAi  \\
 $ALYA_PATH/Alya.x dlrJet --name DIRIC  \\
 : \\
-np $ALYAj \\
 $ALYA_PATH/Alya.x solidDLR --name NEUMA 

# ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016/Executables/plepp03

"""


CMD02="""#  
export TAU_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/TAU/TAU2251/EXEC2/x86_64
export TAU_MAKEFILE=$TAU_PATH/lib/Makefile.tau-icpc-mpi-pdt
export PATH=$PATH:$TAU_PATH/bin

Exec01()
{
  rm F$1xS$2.ppk
  paraprof --pack F$1xS$2.ppk
  rm profile.*
  paraprof -v F$1xS$2.ppk

  pprof -t -f profile.Mean            > F$1xS$2.tau   
  pprof -t -f profile.Mean,AllThreads > F$1xS$2.tau2
 #paraprof F$1xS$2.ppk  
}

ALYA_PATH=XXX_PATH

ALYAi=F_RANK  
ALYAj=S_RANK  

Exec01 $ALYAi $ALYAj  

"""

CMD03="""# 
#Plot_tau.py -F TAUXX*/*.tau

fun01()
{
  cd $1  
  bash TAU.sh
  cd ..
}

""" 

MESH = "/gpfs/vesta-fs0/projects/ATPESC2016/JMAKE/RUNNER/DLR/PARTITIONS/SOLID/"
#PATH = "/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016SEP06/Executables/gLI8M5_hpctk/" 
#PATH = "/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016/Executables/plepp03"
#PATH = "/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016/Executables/tau/"
#PATH = "/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016/Executables/tau_PLEPP03"
#PATH="/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016DEC23/Executables/plepp03" # 2016DIC26 
#PATH="/home/bsc21/bsc21704/z2017/REPOSITORY/ALYA_2017JAN07/Executables/plepp03b" #2017JAN07
PATH="/home/bsc21/bsc21704/z2017/REPOSITORY/ALYA_2017JAN10/Executables/plepp03" #2017JAN14 

EXEC = "Alya.x"


RUN     = True 
RUN     = False 
TIME    = 19  
F_CASE  = "dlrJet"
S_CASE  = "solidDLR"
#
S_RANKS = [9, 8, 7, 6, 5, 4, 3, 2]  
F_RANKS = [9]*len(S_RANKS)
NAME    = "TAU12_07"
#
CMD01   = CMD01_PLEPP  

Dirs = [] 
for i,j in zip(F_RANKS,S_RANKS):
   # numero de procesadores 
   F_RANK = np.power(2,i) 
   S_RANK = np.power(2,j)

   # directorio  
   DIR = "%s_F%04dxS%04d" % (NAME, F_RANK,S_RANK)
   if os.path.exists(DIR): shutil.rmtree(DIR)
   os.makedirs(DIR)
   Dirs.append(DIR)

   # links 
   for case in ["%s.*"%F_CASE, "%s.*"%S_CASE, "*.alya"]: 
     for f in glob.glob(case):
       src = "%s/%s"%(PWD,f)
       dst = "%s/%s"%(DIR,f)
       os.symlink(src, dst)

   # cambio de directorio 
   os.chdir(DIR)
   print os.getcwd()

   # comandos...
   #---------------------------------------------------------# 
   F_ALYA   = "%s %s" % (EXEC, F_CASE)
   S_ALYA   = "%s %s" % (EXEC, S_CASE)

   TO_REPLACE  = {"XXX_PATH":PATH,   "XXX_TIME":str(TIME), "XXX_RANK":str(F_RANK+S_RANK),  
                    "F_CASE":F_CASE,   "F_RANK":str(F_RANK), 
                    "S_CASE":S_CASE,   "S_RANK":str(S_RANK) }
   Replacer("RUNNER.tau", CMD01, TO_REPLACE)


   TO_REPLACE  = {"XXX_PATH":PATH,  "F_RANK":str(F_RANK),  "S_RANK":str(S_RANK) }  
   Replacer("TAU.sh", CMD02, TO_REPLACE)

   if(RUN): os.system("bsub < RUNNER.tau")
   #---------------------------------------------------------# 

   os.chdir(PWD)


#
fout = open("%s.sh"%NAME, "w")
print>> fout, CMD03  
for dir in Dirs: 
  print>> fout, "fun01 %s " % dir  
print>> fout, "#tar -cvf %s.tar %s*/*.ppk" % (NAME, NAME)
print>> fout, "#tar -rvf %s.tar %s*/*.tau" % (NAME, NAME)
print>> fout, "#tar -rvf %s.tar %s*/*.log" % (NAME, NAME)
fout.close()
#

#from time import sleep
#sleep(3.0)
#os.system("bjobs")

#
# 2**np.array(range(4,11))  
#  array([  16,   32,   64,  128,  256,  512, 1024])
# 
# 403311 *8*1/ 2**np.array(range(4,11))
#   array([201655, 100827,  50413,  25206,  12603,   6301,   3150])
#
# 806201 *8*1/ 2**np.array(range(4,11))
#  array([403100, 201550, 100775,  50387,  25193,  12596,   6298])
#

##
#
# ~/z2016/RUNNER/PLEPPs/PLEPP03/plepp03/CHT_ALYA01/TAU.py 
# 
#paraprof --pack S%d.ppk
#pprof
#paraprof profile.* &  
#tauprof -> mean(right click), METRIC=TIME, Value=Inclusive percent, Unit:seconds!! (options-> selec metric -> inclusive, Show Values as Percent) 
#                                    
## RING 
#export TAU_COMM_MATRIX=1
#
### Tracing with Jumpshot (export TAU_TRACE=1)  
# 1)
# tau_treemerge.pl
# 2)
# tau2slog2 tau.trc tau.edf -o tau.slog2
# 3)
# tar -czvf tau.tar tau.slog2
#
## Tracing with VAMPIR:
# tau2otf tau.trc tau.edf app.otf -n 4 -z 
# vampir app.otf 
 
# Scalability
#   https://www.cs.uoregon.edu/research/tau/docs/newguide/bk01pt01ch06s07.html
#   PerfExplorer (https://www.cs.uoregon.edu/research/tau/docs/newguide/bk01ch04s03.html) 
# 1)
# perfdmf_configure --create-default
# 2)
# perfexplorer_configure
# 3)
# perfdmf_loadtrial -a "app" -x "experiment" -n "S032" TAU000032/S32.ppk
# perfdmf_loadtrial -a "app" -x "experiment" -n "S064" TAU000064/S64.ppk  
# ... 
# 4)
# PerfExplorer & 
# 
