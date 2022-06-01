import numpy as np

Data00 = np.loadtxt("XXX_COORDINATES.alya")
Data01 = np.loadtxt("YYY_2014_11_7.dat") 
Data02 = np.loadtxt("HEATF_2014_11_7_7.dat")
Data03 = np.loadtxt("TEMPE_2014_11_7_7.dat")

Tb = Data01[:, 0]/Data01[:, 3]
qw = Data02[:,-1]
Tw = Data03[:,-1]

T0 = 273.0         # [K] 
Ts0= 273.0*3       # [K] 
D  = 2.0           # [m]
k  = 1.2e-6        # [W/m/K]
hx = qw/(Tw - Tb)  # <- q_w = k dt/dr|_{r=R} = h_x ( T_w - Tb )
Nu = hx*D/k    # 

rho = 1.2          # [kg/m3]
Um  = 3.5e-4       # [m/s]
mu  = 1.8e-5       # [Pa.s]
cp  = 100.5        # [J/Kg/K]
Pr  = mu*cp/k             # Prandl 
Re  = rho*D*Um/mu         # Reynolds
Br  = mu*Um**2/k/(T0-Ts0) # Brinkman  
x   = Data00[:,-1]/(D*Re*Pr) 

print "Re:", Re, ",", 
print "Pe:", Re*Pr, ",", 
print "Br:", Br


Data = []
Data.append( x  ) 
Data.append( Tb ) 
Data.append( Tw ) 
Data.append( qw ) 
Data.append( Nu ) 
#Data.append( ) 

np.savetxt("nusselt_number01.dat", np.transpose(Data))

