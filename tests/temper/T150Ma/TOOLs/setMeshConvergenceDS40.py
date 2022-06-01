import os
ROOT = os.getcwd()

#-------------------------------------------------------------------------------------#
import numpy as np

#-------------------------------------------------------------------------------------#
bashCommand01 = """
Replacer01.py \
-I ../../BASE01/ \
-O %s \
-d'{
    "XXX_NAME":"%s",
    "XXX_DT":"%e", 
    "XXX_REY":"%f", 
    "XXX_VEL":"%f", 
    "XXX_DEN":"%f", 
    "XXX_VIS":"%e", 
    "XXX_SPH":"%e", 
    "XXX_CON":"%e", 
    "XXX_TIN":"%f", 
    "XXX_TWALL":"%f", 
    "XXX_MESH":"%s",
    "XXX_ST":"%f", 
    "XXX_FR":"%f", 
    "XXX_AR":"%f",   
    "YYY_VIS":"%e",  
    "YYY_CON":"%e",  
    "YYY_SPH":"%e", 
    "ZZZ_CONDU":"%e" 
   }'
"""


Rair     = 8.314621
Mair     = 0.0289531
gamma    = 7.0/5 

Vel      = lambda _Re, _rho, _L, _mu :  _Re / ( _rho * _L / _mu )
Vis      = lambda _Re, _rho, _L, _U  : _rho * _U * _L / _Re  
Vel      = lambda _Ma, _T            : _Ma  * np.sqrt( gamma * Rair / Mair * _T ) 
Rho      = lambda _P,  _T            : _P  / ( _T * Rair / Mair  )   
Kappa    = lambda _mu, _cp, _Pr      : _mu * _cp / _Pr  

# Phys. Fluids, Vol. 15, No. 7, July 2003  
SUTHE01  = lambda _T, _C1, _C2    : _C1 * _T**(3.0/2) / (_T+_C2)
CONDU01  = lambda _T, _B1, _S1    : _B1 * _T**(3.0/2) / ( _T + _S1 * (1e-12)**(1.0/_T) )
#
MU01     = lambda _T: SUTHE01(_T, 1.4500e-6, 110.0)
K01      = lambda _T: CONDU01(_T, 0.6325e-5, 245.4)

Ti       =     288.15  
P0       =  101325.00  
cp_in    =    1016.56  
#  
D        =  7e-3        
Ma       =  0.010 # 0.025  
Pr       =  0.71  
times    =  1.0/60 
Strouhal =  0.20  
#
if(1): 
  PATH   = "/home/bsc21/bsc21704/z2016/RUNNER/BENCHMARKS02/SUNDEN02/MESHES03/"
  MESHES = {
           "COARSE01/H400L325M05T/h400l325m05t":1 #, "CERFACS01/Cerfacs02":2  
           }  # PLEPP02.sh -> B  
#
Ksf    =  1.0 
TSs    = [ 1.5 ] #[ 298.15/288.15 ] 
REs    = [ 029.6, 040.0, 080.0, 100.0, 120.0, 148.3, 200.0 ] 
Arates = [  0.0  ]
Frates = [  0.0  ]

for m,Fr in enumerate(Frates):
  for l,Ar in enumerate(Arates): 
    for Mk,k in MESHES.iteritems():
      for j,Ts in enumerate(TSs):
        for i,Re in enumerate(REs): 
          v_in      = Vel(Ma,Ti)  
          rho_in    = Rho(P0,Ti)  
          mu_in     = Vis(Re,rho_in,D,v_in)  
          k_in      = Kappa(mu_in,cp_in,Pr)
          dt        = times * D / v_in 
          Tw        = Ts * Ti 
          #
          DirName   = "R%04d" % Re 
          mesh      = PATH + Mk 
          print DirName,  
          print " Re:%f, Vel:%f, Rho:%f, mu:%e, dt:%e, Ti:%f, Tw:%f, Fr:%f, Ar:%f" % (Re, v_in, rho_in, mu_in, dt, Ti, Tw, Fr, Ar)  
          # 
          rho0      = rho_in    
          mu0       = MU01(Ti) 
          k0        = K01( Ti)  
          cp0       = Pr * k0 / mu0  
          COMMAND01 = bashCommand01 % ( DirName, DirName, dt, Re, v_in, rho_in, mu_in, cp_in, k_in, Ti, Tw, mesh, Strouhal, Fr, Ar, mu0, k0, cp0, k_in*Ksf)
          os.system( COMMAND01 )
          os.chdir( DirName )
          print os.getcwd() 
          os.system("bsub < RUNNER.sh") 
          os.chdir(ROOT)

print "OK! "


"""
#cp /home/bsc21/bsc21704/z2016/RUNNER/BENCHMARKS02/SUNDEN01/CPLNG17_02/JMSHI01/setMeshConvergenceDS40.py . 

Exec.py -C "AlyaClean.x" -F "R0148T12M01G0*"
Exec.py -C "Alya2pos.x vortex2D" -F "R0148T12M01G0*"
Exec.py -C "Alya2pos.x Interior01" -F "R0148T12M01G0*"
Exec.py -C "bsub < RUNNER.sh" -F "R0148T12M01G0*"
Exec.py -C "Plot_sets.py -F 'vortex2D-boundary.tem.set'" -F "R0148T12M01G0*" 
Exec.py -C "Plot_sets.py -F 'Interior01-boundary.tem.set'" -F "R0148T12M01G0*"
Plot_multignuplot.py -F 'R0148T12M01G0*/vortex2D_tem_set004.dat' -X '(log($1))' -Y '($3)' '($4)'
Plot_multignuplot.py -F 'R0148T12M01G0*/Interior01_tem_set002.dat' -X '(log($1))' -Y '($3)' '($4)'

Plot_multignuplot.py -F 'XXX*/R0148T15M01_01/vortex2D-boundary.tem.set' -X '(($1==4)?($0):(1/0))' -Y '($4/(3.14159*7e-3))' -Ry ':4e4' -E "bash bash.sh"
Plot_multignuplot.py -F 'XXX*/R0148T15M01_01/vortex2D_tem_set004.dat' -X  '($1)' -Y  '($4/(3.14159*7e-3))' '($3)' -Rx ':0.3' -Ry ':4e4' ':'

gnuplot> plot [][0:10] "R0148T15M01_02/vortex2D-boundary.tem.set" u (($1==4)?($0):(1/0)):($4/(3.14159*7e-3)*7e-3/(409.5-273.0)/0.2860763) w lp t "Nu", 770.937/(3.14159*7e-3)*7e-3/(409.5-273.0)/0.2860763

"""
