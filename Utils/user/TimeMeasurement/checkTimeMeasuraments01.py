#!/usr/bin/python
import os
import glob
import fileinput
import sys
import time
import numpy as np

PWD = os.getcwd()
"""
+2016DIC30. CREATED
+2017JAN10. FROM: 
            /Users/poderozo/Google Drive/z2016/Drafts01/ART03_2016OCT/DATA/CHT_PLEPP04
"""
#-------------------------------------------------------------------------------------#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-F', action='store', dest='DIRECTORIES',
                    default=[], type=str, nargs='+',
                    help='USAGE: -F "DIR*" ',
                    )

#(options, args) = parser.parse_args()
options = parser.parse_args()

if(
    options.DIRECTORIES != []
  ):
  pass
else:
  parser.print_help()
  exit(1)

#-------------------------------------------------------------------------------------#
NOSE_MODUL  = 'MODULX'
NOSE_TASKS  = 'TASKSX'

TASKs = {-1:'------',
          1:'REAPRO',  2:'TURNON',  3:'INIUNK',
          4:'TIMSTE',  5:'BEGSTE',  6:'DOITER',
          7:'CONCOU',  8:'CONBLK',  9:'NEWMSH',
         10:'ENDSTE', 11:'FILTER', 12:'OUTPUT',
         13:'TURNOF', 14:'BEGITE', 15:'ENDITE',
         16:'MATRIX', 17:'DOOPTI', 18:'ENDOPT',
         19:'BEGZON', 20:'ENDZON',
         21:' SETDT', 22:'SETMES', 23:' ALLYA', 24:'DYNAMI'
        }

MODULEs = {  0:'KERNEL',  1:'NASTIN',  2:'TEMPER',
             4:'TURBUL',  7:'ALEFOR', 10:'SOLIDZ',
            19:'CHEMIC',  5:'PARALL', 30:'KERMOD',
            31:' SETDT', 32:'SETMES', 33:' ALLYA', 34:'DYNAMI'
          }

def longest_substr(lst):
    longest = None
    for word in lst:
        for i in range(len(word)):
            for j in range(i+1, len(word)+1):
                if ((longest is None or (j - i > len(longest))) and
                    sum(word[i:j] in w for w in lst) > 1):
                    longest = word[i:j]
    return longest

#--------------------------------------------------------------------------||--#
class ReadSizesFile():
  def __init__(self, _filename):
    fname       = os.path.basename(_filename)
    self.fname, file_extension = os.path.splitext(fname)

    self.TimesxRanks = None 
    self.DATA01      = np.loadtxt( _filename )
    print "|_'%s'" % _filename,  self.DATA01.shape, ",", 

    self.MPI_SIZE  = self.DATA01.shape[0]

    self.coupled   = self.DATA01[0:,1]
    self.Ncoupled  = np.count_nonzero( self.coupled ) 
    self.Idcoupled = np.nonzero( self.coupled )[0] 
    self.coupled   = self.coupled[self.Idcoupled]

    print "Ncoupled/MPI_SIZE:", self.Ncoupled,"/",self.MPI_SIZE, "=", 1.0*self.Ncoupled/self.MPI_SIZE  


  def getData(self):
    return self.DATA01[:,1].copy() 


  def getCoupled(self):
    return self.Idcoupled  


  def __XXX__(self, TimesxRanks): 
    self.TimesxRanks = TimesxRanks 


  def plotDicAsHistGp(self, TimesxRanks, name=''):

    fname  = self.fname
    fname  = fname.lower()

   #TimesxRanks = self.TimesxRanks.copy() 
   #print self.Ncoupled, TimesxRanks.shape, TimesxRanks.min(), TimesxRanks.max() 

    toSave = np.vstack( (self.coupled, TimesxRanks) )
    toSave = toSave.T 

    fout = "%s%s.%s" % (fname, "_hist01", "dat")
    with open(fout,'w') as f: 
      np.savetxt(f, toSave, fmt='%e')
      print " |_>'%s/%s'" % (".", f.name)

    X = toSave[:,0]
    Y = toSave[:,1]
    H, xedges, yedges = np.histogram2d( X,Y )
    H = H.T  # Let each row list bins with common y range.
    X, Y = np.meshgrid(xedges, yedges)


    fout = "%s%s.%s" % (fname, "_hist02", "dat")
    with open(fout,'w') as f:
      for x1,y1,h1 in zip(X, Y, H):                             
        for x2,y2,h2 in zip(x1, y1, h1):
          print>> f, " %e %e %e " % (x2,y2,h2)  
        print>> f

      print " |_>'%s/%s'" % (".", f.name)

    return 




#--------------------------------------------------------------------------||--#
class prettyfloat(float):
    def __repr__(self):
      return "%0.4e" % self

#--------------------------------------------------------------------------||--#
class ReadTimesFile():
  def __init__(self, _filename, _path, nonRoot=False, Idx=None, _debug=False):
    self.NonSum    = [' ALLYA', 'SETMES', ' SETDT']
    self.NonPlot   = [' ALLYA', 'SETMES', ' SETDT', ' TOTAL', 'NRANKi']
    self.percetile = [0.0, 25.0,75.0,50.0, 100.0]


    self.nFiles = 0
    fname       = os.path.basename(_filename)
    self.fname, file_extension = os.path.splitext(fname)
    #self.path   = os.path.dirname(_filename)
    self.path   = _path  

    self.IdxByTimes = {}
    self.IdxByModul = {}
    self.IdxByTasks = {}

    self.DATA01      = np.loadtxt( _filename )             #  1  2  3  4  5          1 + 5          2 + 5, ...,        n + 5
    self.TimesxRanks = self.DATA01[:,6-1:].copy()          # [b, i, m, t, w, TimesxRanks_1, TimesxRanks_2, ..., TimesxRanks_n]  
    self.nRanks      = self.TimesxRanks.shape[1]
    self.TimesxRanks = None 

    if nonRoot and Idx is None: 
      self.TimesxRanks = self.DATA01[:,2+5-1:].copy()      # [TimesxRanks_2, TimesxRanks_3, ..., TimesxRanks_n]   
      print "|_ROOT REMOVED!!:", self.TimesxRanks.shape 

    if not nonRoot and not Idx is None:
      self.TimesxRanks = self.DATA01[:,1+5-1:].copy()      # [TimesxRanks_1, TimesxRanks_2, ..., TimesxRanks_n] 
      self.TimesxRanks = self.TimesxRanks[:,Idx]
      print "|_ IDX REMOVED!!:", self.DATA01.shape, "->", self.TimesxRanks.shape

    if nonRoot and not Idx is None:
      print "|_ ERROR: nonRoot and Idx!!"
      exit(1)

    print "|_'%s'" % _filename,  self.DATA01.shape
    self.__checkData__( nott=1, debug=_debug) 


  def getData(self):
    return self.DATA01[:,5:].copy()


  def __checkData__(self, nott=-1, debug=False):
    if(debug):    print " # ittim, modul, current_task, current_when, MPI_SIZEx(F) " 
    if(nott!=-1): print "|_STEP:%d REMOVED!!" % (nott) 

    Dic = {}
    for idx,Data01 in enumerate(self.DATA01):
      b,i,m,t,w  = Data01[[0,1,2,3,4]].astype(int)
      if not i == nott:
        modul      = MODULEs.get(m,NOSE_MODUL)
        tasks      =   TASKs.get(t,NOSE_TASKS)

        self.__fillDic__(self.IdxByTimes,     i, idx) 
        self.__fillDic__(self.IdxByModul, modul, idx)
        self.__fillDic__(self.IdxByTasks, tasks, idx)

        #if(debug):         print "  ittim, modul, task:", str(i).zfill(3), modul, tasks
        print "  ittim, modul, task:", str(i).zfill(3), modul, tasks

    ## Row selection / [Time/Modul/Task] 
    self.IdxByTimes   = { k:np.array(v) for k,v in self.IdxByTimes.iteritems() }
    self.IdxByModul   = { k:np.array(v) for k,v in self.IdxByModul.iteritems() }
    self.IdxByTasks   = { k:np.array(v) for k,v in self.IdxByTasks.iteritems() }
    #
    ## Filtrando datos...  
    ## eliminando time==0 en Times (es KERMOD) 
    if self.IdxByTimes.has_key(0): self.IdxByTimes.pop(0)
    ## eliminando ' ALLYA', 'SETMES' de Times pues estos NO son parte de los pasos de tiempo
    for k,v in self.IdxByTimes.iteritems(): 
      for toDelete in [' ALLYA', 'SETMES']: 
        v = np.extract(v!=self.IdxByModul[toDelete], v) 
      self.IdxByTimes[k] = v


    ## Gather data per processor / [Time/Modul/Task] 
    self.TimesByRanks = self.__verticalReduction__( self.IdxByTimes, self.TimesxRanks[:,:]                     )
    self.ModulByRanks = self.__verticalReduction__( self.IdxByModul, self.TimesxRanks[:,:], nonSum=self.NonSum )
    self.TasksByRanks = self.__verticalReduction__( self.IdxByTasks, self.TimesxRanks[:,:], nonSum=self.NonSum )

    ## Processors statistic / [Time/Modul/Task]. [percentiles, MEAN, TOTAL] 
    self.TimesStatist = self.__horizontalReduction__( self.TimesByRanks, debug=True )
    self.ModulStatist = self.__horizontalReduction__( self.ModulByRanks, debug=True ) 
    self.TasksStatist = self.__horizontalReduction__( self.TasksByRanks, debug=True ) 

    #self.plotDicAsHistGp( self.TasksStatist ) 
 
    return 


  def __fillDic__(self, Dic, key, val):  
    if key in Dic:
        Dic[key].append( val )
    else:
        Dic[key] = [ val ]
    return 


  def __verticalReduction__(self, Dic, Data, nonSum=None, debug=False): 
    Aux = { k:Data[v,:].sum(axis=0) for k,v in Dic.iteritems() } 

    if nonSum is None:  
      pass 
    else: 
       yesSum        = list( set(Aux.keys()).difference(nonSum) )
       total         = [ Aux[k] for k in yesSum ] 
       total         = sum(total)
       Aux[" TOTAL"] = total  
       Aux["NRANKi"] = self.nRanks  

    if(debug): self.__printDic__(Aux)

    return Aux  


  def __horizontalReduction__(self, Dic, debug=False):
    Aux = { k:self.__horizontStatistic__(v) for k,v in Dic.iteritems() }
    if(debug): self.__printDic__(Aux)

    return Aux


  def __horizontStatistic__(self, Data, percetile=None, debug=False ):
    if percetile is None: percetile = self.percetile
    Aux         = { p:np.percentile(Data,p) for p in percetile }
    Aux["MEAN"] = np.average(Data)  

    if(debug): self.__printDic__(Aux) 

    return Aux


  def __printDic__(self, Dic):
#    for i,(k,v) in enumerate(Dic.iteritems(),0):
    for i,k in enumerate( sorted(Dic.keys()) ,0):
      v = Dic[k] 
      print "%s" % str(i+1).zfill(2), k, 
      if   isinstance(v,dict): 
        print {k2:"% 0.4e"%v2 if isinstance(v2,float) else v2 for k2,v2 in v.iteritems()}
      elif isinstance(v,list):
        print ["% 0.4e"%v2 if isinstance(v2,float) else v2 for v2 in v]
      else:
        print v 

    print 


  def plotDicAsHistGp(self, Dic ):
    AuxDic = self.__Dic2Lists__( Dic )
    #self.__printDicTrans__(AuxDic)

    header = ["PERCEN"] + Dic.keys()  
    fname  = self.__saveDic__(AuxDic, header)
    self.__gpHistFromDic__( header, fname, nonPlot=self.NonPlot )

    return AuxDic 


  def __printDicTrans__(self, Dic):
    keys  = Dic.keys() 
    vals  = Dic.values()                  # values 
    valsT = [list(x) for x in zip(*vals)] # Trasnposed values 
    print keys 
    for i,v in enumerate(valsT): 
      print map(prettyfloat,v)  


  def __Dic2Lists__(self, Dic):
    AuxDic = {}
    for i1,(k1,v1) in enumerate(Dic.iteritems(),0):
      if isinstance(v1,dict):
        for i2,(k2,v2) in enumerate(v1.iteritems(),0): 
          self.__fillDic__(AuxDic, k2, v2)
      if isinstance(v1,list):
        pass

    return AuxDic 


  def __saveDic__(self, Dic, colNames=None):
    """
    Input:  Dic    = {A:[a1,a2,...,an], B:[b1,b2,...,bn],...,Z:[c1,c2,...,cn]} 
            Header = [name1, name2, ..., namen] 
    Output: File  
            # name1 name2, ..., namen  
            A    a1    a2  ...     an  
            B    b1    b2  ...     bn  
            ... 
            Z    z1    z2  ...     zn 
    """
    toSave = [] 
    for i,(k,v) in enumerate(Dic.iteritems(),0):
      toSave.append( [k] + v ) 

    fname = "%s.%s" % (self.fname+"_"+str(self.nFiles).zfill(3), "dat")
    fname = fname.lower() 
#    fname = "%s/%s" % (self.path, fname)
    fout  = open(fname, "w")
    if not colNames is None: 
      header = [ "'%d:%s'" %(i+1,k) for i,k in enumerate(colNames) ]
      header = " ".join(header); #print header 
      print>> fout, "#", header  

    for row in toSave: 
       for v in row: 
         print>> fout, v,
       print>> fout  

    print " |_>'%s/%s'" % (self.path, fout.name) 
    fout.close()  

    return fout.name


  def __gpHistFromDic__(self, colNames, gpFile, nonPlot=None):
    yesKeys  = list( set(colNames).difference(nonPlot) )
    colKeys  = { i+1:k   for i,k in enumerate(colNames,0) }
    colNums  = {   k:i+1 for i,k in enumerate(colNames,0) } 
    toPlot   = [ colNums[k] for i,k in enumerate(yesKeys,0) ] 
    toPlot   = sorted(toPlot) 
    #print toPlot 

    msn      = ""
    msn     += "set style data histograms\n"
    msn     += "set style histogram rowstacked; set boxwidth 1 relative \n"
    msn     += "#set style histogram cluster gap 1\n"
    msn     += "set style fill solid 1.0 border -1\n"
    msn     += "set datafile separator ' ' \n"  
    msn     += "set key outside \n"
    msn     += "set grid \n"
    msn     += "set multiplot layout 1,1 \n"
    msn     += "set yrange [:] \n"
    msn     += " \n"

    msn     += "plot '%s/%s' \\\n" % (self.path, gpFile)
    for i,v in enumerate(toPlot,0): 
      k = colKeys[v] 
      #print i, v, k  
      if(i==0): 
        continue  
      if(i==1):
        msn += "   using %d:xtic(1) t '%s' ,\\\n" % ( v,k ) 
      elif(i==len(toPlot)-1): 
        msn += "'' using %d         t '%s'    \n" % ( v,k ) 
      else: 
        msn += "'' using %d         t '%s' ,\\\n" % ( v,k ) 

    msn     += " \n"
    msn     += "unset multiplot \n"
    msn     += "pause -1 \n"
    #print msn 

    fname, fextension = os.path.splitext(gpFile)
    fname = fname + ".gp"
    fout  = open(fname, "w")
    print>> fout, msn 
    fout.close()
    print " |_>'%s/%s'" % (self.path, fout.name)
 
    return 


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
np.set_printoptions(precision=2)
#
DIRs  = glob.glob(options.DIRECTORIES[0])
DIRs  = sorted(DIRs)
DIRs  = [ dir for dir in DIRs if os.path.isdir(dir) ]
print DIRs

Q1s = []  
for dir in DIRs:
  PWD = os.getcwd()
  os.chdir(dir)
  path, localDir = os.path.split(os.getcwd())

  quartile50th = []   
  for Type in ["DIRIC", "NEUMA"]: 
    FILEs = sorted( glob.glob(Type+"*") )
    print "+'%s'" % localDir, FILEs, "<-"

    if( len(FILEs)>0 ): 
      FileSta = [f for f in FILEs if "sta" in f]
      FileTms = [f for f in FILEs if "tms" in f]

      if( len(FileSta)==0 ):
        print "ERROR: '%s' !!" % ("01, EXIT")
#        exit()
      if( len(FileTms)==0 ):
        print "ERROR: '%s' !!" % ("02, EXIT")
#        exit()

      S = None
      if( len(FileSta)>0 and len(FileTms)>0 ):
        S = ReadSizesFile(FileSta[0])
        T = ReadTimesFile(FileTms[0], dir, nonRoot=False, Idx=S.getCoupled() )
        tasksStatist = T.plotDicAsHistGp( T.TasksStatist )
        quartile50th.append( (T.nRanks,tasksStatist[50.0]) )
        S.plotDicAsHistGp( T.TasksByRanks[' TOTAL'], FileTms[0] )


      T = None 
      if( len(FileTms)>0 ): 
        T = ReadTimesFile(FileTms[0], dir, nonRoot=True,  Idx=None           ) 
        tasksStatist = T.plotDicAsHistGp( T.TasksStatist ) 
        quartile50th.append( (T.nRanks,tasksStatist[50.0]) )
       #S.plotDicAsHistGp( T.TasksByRanks[' TOTAL'], FileTms[0] )

      print 

  if(len(quartile50th)>0): Q1s.append( quartile50th )  

  os.chdir(PWD)

ToSave = [[],[]] 
for i,q1s in enumerate(Q1s): 
    Dranks, Ddata = 0, []
    if len(q1s)>0:  
      Dranks, Ddata = q1s[0] 

    Nranks, Ndata = 0, []  
    if len(q1s)>1: 
      Nranks, Ndata = q1s[1]

    Tranks =  Dranks  + Nranks 
    Ddata  = [Tranks] + Ddata + [Nranks]
    Ndata  = [Tranks] + Ndata + [Dranks]

    Ddata  = [ float(val) for val in Ddata]
    Ndata  = [ float(val) for val in Ndata]

    ToSave[0].append( np.array(Ddata) )
    ToSave[1].append( np.array(Ndata) )

with open("quartile50th.direc",'w') as f: 
 for array in ToSave[0]:
   for val in array: print>> f, val,
   print>> f  
 f.close() 

with open("quartile50th.neuma",'w') as f: 
 for array in ToSave[1]:
   for val in array: print>> f, val,
   print>> f
 f.close()
 
