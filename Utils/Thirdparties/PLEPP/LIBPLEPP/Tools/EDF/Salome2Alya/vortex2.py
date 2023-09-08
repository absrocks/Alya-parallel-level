import numpy as np 


class Vortex:
  def __init__(self): 
    """
    Isentropic flow: dS = 0 
                     dS = dQ/T    (reversible processes)
                     dS=0 -> dQ=0 (adiabatic processes)
                     
    dU = -pdV + dQ ; Cv = dU/dT             -> Cv dT = -pdV \__ gamma = Cp/Cv = -dp/p / dV/V 
    dH =  dU + d(pV) = Vdp + dQ; Cp = dH/dT -> Cp dT =  Vdp /   with pV = RT and Cp = Cv + R (ideal gas)
                                                                1/R = gamma M**2 = rho T /p
    """
    self.a     =-1.0           # sound speed
    self.C     = 1.0           # cte 
    self.Rc    = [0.0,0.0,0.0] # vortex center 
    self.Ro    = 1.0           # vortex radio  


    print "\n+Vortex:"
    
    
  def vel(self, V, R):
    """
    psi(x,y) = gamma exp[ -1/2 (|R-Rc|/a)**2 ]
    Ux = d phi/dy
    Uy =-d phi/dx
    Uz = 0
    Umax = gamma/a/sqrt(e)
    """
    rx = 1
    ry = 2
    rz = 0 
    Rx = R[:,rx] - self.Rc[rx]
    Ry = R[:,ry] - self.Rc[ry]
    Ro2   = self.Ro*self.Ro
    modR2 = Rx*Rx + Ry*Ry  
    modR  = np.sqrt(modR2)
    self.PSI  = modR2 
    self.PSI /= 2.0*Ro2
    self.PSI  = np.exp( -self.PSI ) * self.C 
    #self.PSI *= modR / self.a*self.a 

    V[:,rz] +=  0.0  
    V[:,rx] += -Ry * self.PSI / Ro2
    V[:,ry] +=  Rx * self.PSI / Ro2 

    a2 = self.a*self.a 
    self.PSI *=  self.PSI     # Psi**2
    self.PSI *= -0.5/(a2*Ro2) # -0.5*1.0/(a*Ro)**2 * Psi**2 

    """
    Rho/Rho0 = exp[ -0.5*(a*Ro)**-2 * Psi**2 ] 
    """
    self.Rho = np.exp( self.PSI ) 
    #self.Rho = self.PSI 

    print "|_Vortex.vel %s m/s," % str(np.average(V, axis=0)), 
    #print np.max(V, axis=0)
    #print np.min(V, axis=0)
    print "std %s " % str(np.std(V, axis=0))


  def density(self, Rho, Rho0=1.2): 
    Rho[:] = Rho0 * self.Rho
    print "|_Vortex.density %f Kg/m3" % np.average(Rho) 
    

  def temp(self, T, T0=300.0, gamma=1.4):
    """
    T/T0 = (rho/rho0)**(gamma-1)
    """
    T[:] = T0*(self.Rho)**(gamma-1.0)
    print "|_Vortex.temp %f K" % np.average(T) 
    

  def sound(self, T0, gamma):
    Rgas = 287.0 
    self.a = np.sqrt(gamma*Rgas*T0) 
    print "|_Vortex.sound_speed %d [m/s]" % self.a
    
    return self.a 
  
      
    
