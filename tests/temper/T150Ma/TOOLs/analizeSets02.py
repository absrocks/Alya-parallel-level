#!/usr/bin/env python
import sys
import os 
import os.path
import numpy as np 
import json 
import time
import multiprocessing
#-------------------------------------------------------------------------||---#
full_path = os.path.realpath(__file__)
full_path = os.getcwd()
path, filename = os.path.split(full_path)

date = time.strftime("%Y%b%d")
date = date.upper()

n_cpus = multiprocessing.cpu_count()

#-------------------------------------------------------------------------||---#
def LookinForFile( fname, extension ):
  path, filename = os.path.split( fname )
  filename, file_extension = os.path.splitext( filename )
  filename = "%s/%s%s" % (full_path, filename, extension)
 
  result = False  
  import glob
  fname = glob.glob( filename )
  if(len(fname)<1):
    print "WARNING: there no exist '%s' files! \n\n" % ( extension ) 
  elif(len(fname)>1):
    print "WARNING: there exist %d '%s' '%s' files, choose one!! \n\n" % (len(fname), extension, " ".join(fname) )
  else:
    fname = fname[0]
   #print "|_\'%s\'" % fname
    result = fname 

  if(not result): exit(1)
  return result 

#-------------------------------------------------------------------------||---#
def AnalizeSetFile( _Fin ): 
 #print _Fin 
  key  = "Time"
  Time = []
  with open(_Fin) as f:
    for i,line in enumerate(f):
      line = line.strip()
      line = line.split()
      if key in line:
        Time.append( eval(line[-1]) )
  Time = np.array(Time)

  Data01 = np.loadtxt(_Fin)
  Wits   = {}
  for i,data in enumerate(Data01):
    idx  =  data[0]
    if idx in Wits:
      Wits[idx].append(i)
    else:
      Wits[idx]  = [i]

  Data02 = [] 
  for key in sorted(Wits.keys()):
    value = Wits[key]
    Data01[value,0] = Time[:]
    Data02.append( Data01[value] )

  Data02 = np.array(Data02)
# print Data02.shape 

  return Data02  


#-------------------------------------------------------------------------||---#
def KerDat( fname, Keys ):
    data = open(fname, "r")
    lines = data.readlines()
    data.close()

#    Dict = { key:None for key in Keys} 
    Dict = {} 
    for key in Keys:
      Dict[key] = None
      for line in lines:
        line = line.split("#")[0]
        line = line.replace("$", "")
        ok   = False
        if(line.find(key)>0): 
          line = line[:-1]
          line = line.split(key)
          line = line[1]
          line = line.replace(' ', '')
          Dict[key] = line 

    return Dict 

#-------------------------------------------------------------------------||---#
def GetStrouhal( _signal, dt, d0, v0 ):  
    n_sample = _signal.shape[0] 
    dB    = np.abs(np.fft.fft( _signal ))
    dB    = 20 * np.log10(dB)
    idx   = np.argmax(dB)
    freqs = np.fft.fftfreq(n_sample+0, d=dt) # sample frequencies
    freq  = np.abs( freqs[idx] )
    St    = Strouhal(freq,d0,v0)
   #print "  |_Fr:%f Hz, St:%f, dt:%e" % (freq,St,dt)    
    return St  

#-------------------------------------------------------------------------||---#
def GetRangesIds( _signal ): 
  roots   = [ i for i,dummy in enumerate(_signal[:-1]) if(_signal[i]*_signal[i+1]<=0) ]
  n_roots = len(roots)
  Idx     = [ i for i in range(0,n_roots,2)]
  Idx     = [ (roots[top],roots[bottom]) for top,bottom in zip(Idx[:-1], Idx[1:]) ]

  return Idx  

def Analize( _signal, _Ids ): 
  Avr=[]; Rms=[]; Max=[]; Min=[]; Tra=[];
  for idx in _Ids:
    section = _signal[idx[0]:idx[1]]
    Avr.append( np.mean(section) )
    Rms.append( np.sqrt(np.dot(section,section)/section.shape[0]) )
    Max.append( np.amax(section) )
    Min.append( np.amin(section) )
    Tra.append( np.trapz(section) )

  Min = np.array(Min)   
  Avr = np.array(Avr) 
  Max = np.array(Max)   
  Rms = np.array(Rms) 
  Tra = np.array(Tra)   

 #print "  |_Min:%f, Avr:%f, Max:%f, Rms:%f" % ( Min.mean(), Avr.mean(), Max.mean(), Rms.mean() ) 
  print "  |_Min:%f, Avr:%f, Max:%f, Rms:%f, Tra:%f " % ( Min[-1], Avr[-1], Max[-1], Rms[-1], Tra[-1] )

  result = np.zeros( (len(_Ids),6) )
  result[:,0] = range( len(_Ids) )
  result[:,1] = Min[:] 
  result[:,2] = Avr[:]  
  result[:,3] = Max[:] 
  result[:,4] = Rms[:]  
  result[:,5] = Tra[:] 

  return result 

#-------------------------------------------------------------------------||---#
DragLift = lambda _force, _rho, _d, _v: _force/(0.5*_rho*_d*_v**2)
Forced   = lambda  _time,   _A, _f    : [ 0.0, _A * np.sin( 2 * 3.14159 * _f * _time ) ]
Strouhal = lambda _f, _d , _v : _f * _d / _v
Nusselt  = lambda _dT, _L, _k, _T2, _T1: np.where( (_T2-_T1)!=0, _dT*_L/_k /(_T2-_T1), 0.0 )

#-------------------------------------------------------------------------||---#

#==================================================================| Parse |===#
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-F", "--File", default="", dest="File")
parser.add_argument('-Y', action='store', dest='Ys',
                    default=[], type=str, nargs='+',
                    help='Columns',
                    )
parser.add_argument('-O', action='store', dest='Operations',
                    default=[], type=str, nargs='+',
                    help='Operations',
                    )
options = parser.parse_args()

#-------------------------------------------------------------------------||---#
if(
options.File == ''     
  ):  
  parser.print_help()
  print 
  sys.exit() 

Fin01 = options.File

#-------------------------------------------------------------------------||---#
main_dat = LookinForFile( Fin01, ".dat")

if(main_dat):
  Dic01 = KerDat(main_dat, ['TIME_STEP_SIZE'])

  dt = Dic01['TIME_STEP_SIZE']
  dt = dt.replace(":","")
  dt = eval(dt)
 #print "  |_dt:%e" % dt

#-------------------------------------------------------------------------||---#
ker_dat = LookinForFile( Fin01, ".ker.dat")

if(ker_dat): 
  Dic01 = KerDat(ker_dat, ['Re =', 'Ma =', 'Pr =', 'Tr =']) 

  Ma, Re, Pr, Tr = None,None,None,None 
  U, Ti, rho, mu, L, mu, cp, k, Tw = None, None,None,None,None,None,None,None,None  
  for k,v in Dic01.items():
    if("Ma" in k):  Ma = eval(v)
    if("Re" in k):  Re = eval(v)
    if("Pr" in k):  Pr = eval(v)
    if("Tr" in k):  Tr = eval(v)

    v = v.replace("/", " ")
    v = v.replace("*", " ")
    v = v.replace("np.sqrt", " ")
    v = v.replace("(", " ")
    v = v.replace(")", " ")  
    v = v.strip()
    v = v.split()
    if("Tr" in k): [Tw, Ti] = [eval(x) for x in v];  
    if("Ma" in k): [U, dummy, dummy, dummy, dummy, Ti] = [eval(x) for x in v]; 
    if("Re" in k): [rho, U, mu, L] = [eval(x) for x in v];  
    if("Pr" in k): [mu, cp, kappa] = [eval(x) for x in v]; 

  print "  |_Ma:%f, Re:%f, Pr:%f, Tr:%f, dt/(L/U):%f" % ( Ma, Re, Pr, Tr, (L/U)/dt ) 


#-------------------------------------------------------------------------||---#
nsi_set  = LookinForFile( Fin01, "-boundary.nsi.set")
if(nsi_set):  
  nsi_data = AnalizeSetFile(nsi_set)

  Set = 3 
  Cl  = DragLift( nsi_data[Set,:,3]+nsi_data[Set,:,6], rho, L, U ) 
  Ids = GetRangesIds( Cl ) 

  print "    |_[%s]" % ("Stro"), 
  St = [ GetStrouhal(Cl[r[0]:r[1]],dt,L,U) for i,r in enumerate(Ids)]
  St = np.array(St)
  print "   |_Hz:%f" % St[-1]

  print "    |_[%s] " % ("Lift"), 
  Cl  = Analize(Cl, Ids) 

  print "    |_[%s] " % ("Drag"),  
  Cd  = Analize( DragLift( -nsi_data[Set,:,2]-nsi_data[Set,:,5], rho, L, U), Ids ) 

  np.savetxt("Cd.dat", Cd)   

#-------------------------------------------------------------------------||---#
tem_set  = LookinForFile( Fin01, "-boundary.tem.set")
if(tem_set): 
  tem_data = AnalizeSetFile(tem_set)

  print "    |_[%s] " % ("Nslt"), 
  Analize( Nusselt(tem_data[Set,:,3]/(np.pi*L),L,kappa,Tw,Ti), Ids )  


#-------------------------------------------------------------------------||---#
#for i,range in enumerate(Ids):
#  init,end = range 
#  print Cd[init:end].shape 


#=========================================================================||===#
