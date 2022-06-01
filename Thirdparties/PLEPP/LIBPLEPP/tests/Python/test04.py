#!/usr/bin/env python
import sys
import os 
import os.path
import numpy as np 
import json 
import time 
import shutil
#-------------------------------------------------------------------------||---#
HEADER  = '\033[95m'
OKBLUE  = '\033[94m'
OKGREEN = '\033[92m'
ENDC    = '\033[0m'
pcolor  = lambda _str: HEADER + '%s'%(_str) + ENDC 

#-------------------------------------------------------------------------||---#
FullPath = os.path.realpath(__file__)
FullPath = os.getcwd()
Path, Filename = os.path.split(FullPath)

date = time.strftime("%Y%b%d")
date = date.upper() 

DIR_OU = FullPath + "/XXX_" + date  

if os.path.isdir(DIR_OU) and os.access(DIR_OU, os.R_OK):
  shutil.rmtree(DIR_OU)
  print "+'%s' removed! " % ( DIR_OU )

os.makedirs( DIR_OU )
#Filename = DIR_OU +"/"+ Filename 

#==================================================================| Parse |===#
from optparse import OptionParser
usage = "\n\t >> python ./%prog [options] arg1 arg2!!\n\t >> PATH_TO_FILES=MESH_PATH/APP_NAME+FILES \n"
parser = OptionParser(usage=usage)
#
parser.add_option("-n", "--name", type="string", nargs=1,
                  dest="app_name", default=(''), help="application name")
#
parser.add_option("-m", "--mesh", type="string", nargs=1,
                  dest="mesh_path", default=('./'), help="mesh path")
#
parser.add_option("-r", "--Reynolds", type="float", nargs=1,
                  dest="Reynolds", default=(-1.0), help="Reynolds") 
#
parser.add_option("-i", "--Ti", type="float", nargs=1,
                  dest="Ti", default=(-1.0), help="Input Temperature")
#
parser.add_option("-w", "--Tw", type="float", nargs=1,
                  dest="Tw", default=(-1.0), help="Wall Temperature")
#
parser.add_option("-d", "--Density", type="float", nargs=1,
                  dest="Density", default=(-1.0), help="Input Density")
#
parser.add_option("-k", "--Conductivity", type="float", nargs=1,
                  dest="Conductivity", default=(-1.0), help="Input Conductivity")
#
parser.add_option(     "--Viscosity", type="float", nargs=1,
                  dest="Viscosity", default=(-1.0), help="Input Viscosity")
#
parser.add_option("-s", "--Strouhal", type="float", nargs=1,
                  dest="Strouhal", default=(-1.0), help="Strouhal Frequency")
#
parser.add_option("-f", "--Forced", type="float", nargs=1,
                  dest="Forced", default=(-1.0), help="Forced Frequency")
#
parser.add_option("-a", "--Amplitude", type="float", nargs=1,
                  dest="Amplitude", default=(-1.0), help="Adimensional Amplitude ")
#
parser.add_option('-x', '--Files', type=str,
                  dest="my_dict", default='{}', help="Files")
(options, args) = parser.parse_args()

my_dictionary = json.loads(options.my_dict)

if(options.app_name=='' or len(my_dictionary)==0 or 
  options.Reynolds==-1.0 or 
  options.Ti==-1.0 or 
  options.Tw==-1.0 or  
  options.Strouhal==-1.0 or 
  options.Forced==-1.0 or
  options.Amplitude==-1.0 or  
  options.Density==-1.0 or 
  options.Conductivity==-1.0 or 
  options.Viscosity==-1.0 
  ):   
  parser.print_help()
  print 
  sys.exit() 

#-------------------------------------------------------------------------||---#
BASE_NAME = options.mesh_path +"/"+ options.app_name  
DATAs     = {}  

for k, v in my_dictionary.iteritems():
  file = BASE_NAME+v 
  if os.path.isfile(file) and os.access(file, os.R_OK):
    Data = np.loadtxt(file)
    print "\t size:%d '%s'" % ( Data.shape[0], file ) 
    DATAs[k] = Data
  else:
    print "Either file is missing or is not readable '%s' " % file 
    sys.exit() 

#=================================================================| Mpi4py |===#
from mpi4py import MPI
print (MPI.get_vendor())
world_comm = MPI.COMM_WORLD

#-------------------------------------------------------------------------||---#
import Commdomm

class PyPlepp: 
  def __init__(self, _type, _namei, _namej, _comm, _dim=3):  
    self.dim = _dim  
    self.CD  = Commdomm.CommDom()
    self.CD.init()
    self.CD.set_app_type(   _type  ) 
    self.CD.set_app_name(   _namei ) 
    self.CD.set_world_comm( _comm  ) 

    self.lcomm = self.CD.set_mpi_comms()
    self.lrank = self.lcomm.Get_rank()
    self.lsize = self.lcomm.Get_size()

    self.wcomm = _comm  
    self.wrank = _comm.Get_rank()
    self.wsize = _comm.Get_size()
 
    self.commij = MPI.COMM_NULL
    if( (self.CD.__get_app_name__() == _namei) and (self.CD.__get_friends__(_namej) == 1) ):
      self.commij = self.CD.get_mpi_commij(_namej)
    if( (self.CD.__get_app_name__() == _namej) and (self.CD.__get_friends__(_namei) == 1) ):
      self.commij = self.CD.get_mpi_commij(_namei)
    self.lcomm.Barrier()


  def SetGeometryI(self, n_vertices, vertex_coords, n_elements, n_vertex_type, vertex_num, vertex_type):
    self.n_elements_i    = n_elements
    self.vertex_num_i    = Commdomm.iarray( n_elements*n_vertex_type )
    self.vertex_type_i   = Commdomm.iarray( n_elements               )

    k = 0
    for i in range(n_elements):
      self.vertex_type_i[i]  = int(vertex_type[i])    
      for j in range(n_vertex_type): 
        self.vertex_num_i[k] = int(vertex_num[i][j])
        k += 1

    self.n_vertices_i    = n_vertices
    self.vertex_coords_i = Commdomm.darray( n_vertices*self.dim )

    k = 0 
    for i in range(n_vertices): 
      for j in range(self.dim): 
        self.vertex_coords_i[k] = vertex_coords[i][j]
        k += 1  


  def SetGeometryJ(self, n_vertices, vertex_coords):
    self.n_vertices_j    = n_vertices
    self.vertex_coords_j = Commdomm.darray( n_vertices*self.dim )

    k = 0
    for i in range(n_vertices):
      for j in range(self.dim):
        self.vertex_coords_j[k] = vertex_coords[i][j]
        k += 1


  def Runner(self):  
    if(self.commij != MPI.COMM_NULL): 
      self.CD.locator_create2(self.lcomm, self.commij, 1e-3)

      self.CD.locator_set_cs_mesh(self.n_vertices_i,
                                  self.n_elements_i,
                                  self.vertex_coords_i,
                                  self.vertex_num_i,
                                  self.vertex_type_i,
                                  self.n_vertices_j,
                                  self.vertex_coords_j, 
                                  self.dim);
      self.CD.save_dist_coords(self.lrank)

    self.n_send = self.CD.get_n_dist_points()
    self.n_recv = self.CD.get_n_interior()

    self.lcomm.Barrier()


  def GetCoordsJ(self):  
    dist_coords_j    = Commdomm.darray( self.n_send*self.dim )
    self.CD.__locator_get_dist_coords__(    dist_coords_j    )
    return self.F2C( self.n_send, self.dim, dist_coords_j ) 


  def Bcast(self, n_send, varij, n_recv): 
    send = Commdomm.darray(n_send) 
    recv = Commdomm.darray(n_recv)

    for i in range(n_send): send[i] = varij[i]
    self.CD.__mpi_sendrecv_real__( send, n_send, recv, n_recv, self.lcomm, self.commij);
    self.CD.__mpi_bcast_real__(                  recv, n_recv, self.lcomm, self.commij);

    varji = np.zeros(n_recv)
    for i in range(n_recv): varji[i] = recv[i] 

    return varji  


  def ExchangeDouble(self, send, n_dof):
    n_send = self.n_send * n_dof
    n_recv = self.n_recv * n_dof 
    var_ij = Commdomm.darray(n_send)
    var_ji = Commdomm.darray(n_recv)

    for i in range(n_send): var_ij[i] = send[i]
    self.CD.__locator_exchange_double_scalar__(var_ij, var_ji, n_dof)

    recv = np.zeros(n_recv)
    for i in range(n_recv): recv[i] = var_ji[i] 

    return recv 


  def F2C(self, n_var, n_dof, var_ji):
    """
    [x1 y1] [x2 y2] ... [xn yn] 
    """
    aux = np.zeros( (n_var,n_dof) )  
    for i in range(n_var): 
      for k in range(n_dof): 
        aux[i,k] = var_ji[i*n_dof + k] 

    return aux.copy()  


  def C2F(self, n_var, n_dof, var_ij):
    aux = np.zeros( n_var*n_dof )
    k = 0
    for i in range(n_var):
      for j in range(n_dof):
        aux[k] = var_ij[i][j]
        k += 1 

    return aux.copy()

#=========================================================================||===#
DragLift = lambda _force, _rho, _d, _v: _force/(0.5*_rho*_d*_v**2)  
Forced   = lambda  _time,   _A, _f    : [ 0.0, _A * np.sin( 2 * 3.14159 * _f * _time ) ]  
Strouhal = lambda _f, _d , _v : _f * _d / _v
Nusselt  = lambda _dT, _L, _k, _T2, _T1: np.where( (_T2-_T1)!=0, _dT*_L/_k /(_T2-_T1), 0.0 )

t_size = int(1e6) 
t_max  = 2e-5 * 20

Re0    = options.Reynolds     #  
T_i    = options.Ti           # 303.0 
T_w    = options.Tw           # 373.0 
f_0    = options.Strouhal     # 560.0  
Frate  = options.Forced 
Arate  = options.Amplitude
#
r_0    = options.Density      # kg/m3, 1.177  
k_0    = options.Conductivity # W/M/K **  
mu_0   = options.Viscosity    # Ns/m2 **  == mu_air x 100  
d_0    = 7e-3                 # m     ** 
# 
v_0    = Re0 / (d_0 * r_0 / mu_0)  
f_1    = Frate * f_0      # Forced frequency  
d_1    = Arate * d_0      # Max amplitude   
#
# ** M-P Errera, F. Duchaine JCP 2016  
#
print "PROPS-> Re:%f Vel:%f Ti:%f Tw:%f Rho:%f fS:%f f0:%f ymax:%f " % ( Re0, v_0, T_i, T_w, r_0, f_0, f_1, d_1)
#=========================================================================||===#
PTS      = DATAs["COORDINATES"][:,1:] 
ELEMENTS = DATAs["ELEMENTS"][:,1:]
TYPES    = DATAs["TYPES"][:,1:]

DIM   = 2 
DOF   = DIM + 1 
IAM   = "PLEPP";
NAMEi = "DIRIC";
NAMEj = "NEUMA";

PP = PyPlepp( IAM, NAMEi, NAMEj, MPI.COMM_WORLD, DIM )
PP.SetGeometryI( PTS.shape[0], PTS, ELEMENTS.shape[0], ELEMENTS.shape[1], ELEMENTS, TYPES )  
PP.SetGeometryJ( PTS.shape[0], PTS ) 
PP.Runner()

fname01 = Filename + "_%s.dat" % date  
if os.path.exists(fname01): os.remove(fname01) 

if(1):
  time    = 0.0;
  displ   = np.zeros( (     2,DIM) )  
  signal  = np.zeros( (t_size,DIM) )  
#  roots   = np.zeros(  t_size      )
  Roots = [ 0 ] 
  for itime in range(t_size):
    #---------------------------------------------------------------------||---#
    dt    = 1.0/PP.Bcast(1, [1.0], 1) 
    time += dt[0];
    print "time:", time, dt[0] 

    #--------------------------------------------------------------| <-FSI |---#
    displace   = Forced( time, d_1, f_1 )
    dummy      = np.zeros( (PP.n_send,DOF) )
    dummy[:,0] = 0.0 #displace[0] - displ[0][0]
    dummy[:,1] = 0.0 #displace[1] - displ[0][1]

    #-----------------------------------------------------------| EXCHANGE |---#
    varij = PP.C2F( PP.n_send, DOF, dummy )
    varji = PP.ExchangeDouble( varij, DOF ) 
    varji = PP.F2C( PP.n_recv, DOF, varji )
    ptsji = PP.GetCoordsJ()

    #--------------------------------------------------------------| FSI<- |---#
    draglift = DragLift( varji[:,:DIM], r_0, d_0, v_0 ) 
    signal[itime,:] = draglift.sum(axis=0) 

    dB    = np.abs(np.fft.fft( signal[:itime+1,1] ))
    dB    = 20 * np.log10( dB )
    idx   = np.argmax(dB) 
    freqs = np.fft.fftfreq(itime+1, d=dt[0]) # sample frequencies
    f_2   = np.abs( freqs[idx] )  
    St    = Strouhal(f_2,d_0,v_0) 
    print "Fr1: %f, St1: %f" % (f_2,St), 

    meanCd, minCl, rmsCl, maxCl = 0.0, 0.0, 0.0, 0.0  
    if(itime>0):
      rmsCl  = np.sqrt( np.dot(signal[:itime,1],signal[:itime,1])/itime )   
      maxCl  = np.amax( signal[:itime,1] )
      minCl  = np.amin( signal[:itime,1] ) 

      meanCd = np.mean( signal[:itime,0] )
      print "meanCd1:%f, rmsCl1:%f , maxCl1:%f " % (meanCd, rmsCl, maxCl),  

      rootCl = signal[itime,1] * signal[itime-1,1] <= 0.0
     #print "rootCl", pcolor(rootCl)
      if(rootCl): Roots.append( itime )

      if( len(Roots)==3 ):  
        t2 = Roots[-1] #   End time   
        t1 = Roots[-3] # Start time 
        # 
        dummy = np.zeros( (t2-t1,3) )  
        dummy[:,0] = np.arange(t1,t2) 
        dummy[:,1] = signal[t1:t2,0]  
        dummy[:,2] = signal[t1:t2,1] 
        np.savetxt(DIR_OU+"/yyy%s.xxx" % ( str(itime).zfill(6) ), dummy )  

        # Lift  
        signal2 = signal[t1:t2, 1 ] 
        Nsignal = signal2.shape[0]
       #print "t2,t1,N:", t2,t1, Nsignal,

        maxCl2  = np.amax( signal2 )
        minCl2  = np.amin( signal2 ) 
        rmsCl2  = np.sqrt( np.dot(signal2,signal2)/Nsignal )

        dB2     = np.abs(np.fft.fft(signal2))
        dB2     = 20 * np.log10( dB2 )
        idx2    = np.argmax( dB2 )
        freqs2  = np.fft.fftfreq(Nsignal, d=dt[0]) # sample frequencies
        f_22    = np.abs( freqs2[idx2] )
        St2     = Strouhal(f_22,d_0,v_0)
        print " Fr2:%f, St2:%f" % (f_22,St2),  

        intgrl2 = np.trapz( signal2 )

        # Drag  
        signal2 = signal[t1:t2, 0 ]
        meanCd2 = np.mean( signal2 ) 
        cosa01 = np.abs(maxCl2 - np.abs(minCl2)) <= 1e-3 
        cosa02 = np.abs(intgrl2) <= 1e-3   
        print " meanCd2:%f, rmsCl2:%f, maxCl2:%f, intgrl2:%f cosa:%d cosa2:%d" % (meanCd2, rmsCl2, maxCl2, intgrl2, cosa01, cosa02) 

        Roots = [t2] 

    displ[0] = displace 
    varji[:,:DIM] = draglift  

    #--------------------------------------------------------------| CHT<- |---#
    nusselt = varji[:,DIM] #np.abs( varji[:,-1] )   
    nusselt = Nusselt( nusselt, d_0, k_0, T_w, T_i )  
    minNu, meanNu, rmsNu, maxNu = 0.0, 0.0, 0.0, 0.0  
    if(itime>0):
      minNu  = np.amin( np.abs(nusselt) )
      meanNu = np.mean( np.abs(nusselt) )
      rmsNu  = np.sqrt( np.dot( np.abs(nusselt),np.abs(nusselt) )/itime )
      maxNu  = np.amax( np.abs(nusselt) )
      print "minNu:%f, meanNu:%f, rmsNu:%f, maxNu: %f" % (minNu, meanNu, rmsNu, maxNu)

    varji[:,DIM] = nusselt 
 
    #---------------------------------------------------------------| SAVE |---#
    if os.path.exists(fname01):
      fout01  = open(fname01, "a")
    else:
      fout01  = open(fname01, "w")
      print>> fout01, "#    1,   2,  3,      4,     5,     6,    7,      8,      9,    10,   11,   12,    13,    14,    15 "  
      print>> fout01, "# time, f_2, St, meanCd, minCl, rmsCl, maxCl, minNu, meanNu, rmsNu, maxNu, drag, lift, dispX, dispY "  
    print>>   fout01,    time, f_2, St, meanCd, minCl, rmsCl, maxCl, minNu, meanNu, rmsNu, maxNu,  

    for dummy in signal[itime,:]: print>> fout01, dummy, 
    for dummy in displace: print>> fout01, dummy,
    print>> fout01
    fout01.close()

    #phi = np.arctan2(ptsji[:,0], ptsji[:,1]) # paraview order!! 
    #Idx = np.argsort(phi)
    fname02 = DIR_OU+"/coords_%s.xxx" % ( str(itime).zfill(6) ) 
    fout02  = open(fname02, "w") 
    print>> fout02, "#    1, 2, 3,   4"
    print>> fout02, "# time, x, y, phi"
    for idx in range(PP.n_send):
      print>> fout02, time,   
      pt  = ptsji[idx]
      for dummy in  pt: print>> fout02, dummy, 
      print>> fout02  
    fout02.close() 

    Ftot = varji[:,:DIM].sum(axis=1) 
    Idx  = np.argsort(Ftot) 
    fname02 = DIR_OU+"/props_%s.xxx" % ( str(itime).zfill(6) )
    fout02  = open(fname02, "w")
    print>> fout02, "#   1,    2,  3,  4,  5 "
    print>> fout02, "# idx, time, fx, fy, nu "
    for i in range(PP.n_recv):
      idx = i #Idx[i]
      var = varji[idx]
      print>> fout02, i, time,
      for dummy in var: print>> fout02, dummy,
      print>> fout02
    fout02.close()
    #---------------------------------------------------------------------||---#
    #---------------------------------------------------------------------||---#

#-------------------------------------------------------------------------||---#

#-------------------------------------------------------------------------||---#


#=========================================================================||===#
#=========================================================================||===#
