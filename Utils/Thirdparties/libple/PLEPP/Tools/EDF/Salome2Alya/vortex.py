import numpy as np 


class Vortex:
  def __init__(self): 
    self.a     = 1.0 # [m]
    self.gamma = 1.0 # [m2/s]
    self.Rc    = [0.0,0.0,0.0] 
    print "\n+Vortex:"
    
    
  def vel(self, V, R):
    """
    psi(x,y) = gamma exp[ -1/2 (|R-Rc|/a)**2 ]
    Ux = d phi/dy
    Uy =-d phi/dx
    Uz = 0
    Umax = gamma/a/sqrt(e)
    """
    self.PSI  = (R[:,1] - self.Rc[1])**2 + (R[:,2] - self.Rc[2])**2 
    self.PSI /= 2.0*self.a*self.a
    self.PSI  = self.gamma*np.exp( -self.PSI )        

    V[:,0] +=  0.0  
    V[:,1] += -(R[:,2] - self.Rc[2]) * self.PSI / (self.a*self.a)
    V[:,2] +=  (R[:,1] - self.Rc[1]) * self.PSI / (self.a*self.a)

    
    self.Vel = V
    print "|_Vortex.vel"


  def press(self, P, rho=1.0, P0=101325.0):
    """
    rho Vj dVi/drj + dP/dri = 0 
    P0  = 101325 N/m2 
    rho(aire) = 1.0 Kg m-3
    rho(water) = 1000.0 Kg m-3
    """
    P[:] = P0 - 0.5*rho*(self.gamma/self.a)**2 * self.PSI
    print "|_Vortex.press"


  def total_energy(self, E, T=300.0, Cv=1210.0): 
    """
    e = Cp T + p/rho = Cv T
    Cv(Water300K) = 4.17000 J cm-3 K-1 
    Cv(Air300K)   = 0.00121 J cm-3 K-1 = 1210 J m-3 K 
    """
    E[:]  = Cv * T
    E[:] += (self.Vel**2).sum(axis=1) 
    
    print "|_Vortex.energy"
    
    
    
    
