import numpy as np 
"""
 ma = F - kx - c v  ->  x'' + (c/m) x'   + (k/m) x = F/m  
                        x'' + (2 eta w0) + w0**2 x = F/m 
"""


k   = 5.790;    # N/m2 
c   = 0.325e-3; # g/s2 -> 1kg/1e3g 
m   = 2.979e-3; # g -> kg 
D   = 0.160e-2; # cm -> m 

k = 10.0 
c =  2.0
m =  1.0 

w0  = np.sqrt(k/m)       + 0.0j 
eta = c/2.0/np.sqrt(m*k) + 0.0j 
print "w0, eta:", w0, eta

x0  = -0.5 #D
v0  =  4.0 #0.0 


# GENERAL SOLUTION  
gamma_p = -w0*eta + w0 * np.sqrt( eta**2 - 1.0+0.0j ) 
gamma_m = -w0*eta - w0 * np.sqrt( eta**2 - 1.0+0.0j) 

x = [x0]
v = [v0]

A =  (gamma_p * (x[0]+0.0j) - (v[0]+0.0j))/(gamma_m - gamma_p) + x[0]+0.0j 
B = -(gamma_p * (x[0]+0.0j) - (v[0]+0.0j))/(gamma_m - gamma_p)
print "A,B:", A, B

x01 = lambda t: A * np.exp( gamma_p*t ) + B * np.exp( gamma_m*t )   

wd  = w0 * np.sqrt(1-eta**2)
x02 = lambda t: ( x[0] * np.cos(wd*t) + (eta*w0*x[0] + v[0])/wd * np.sin(wd*t) ) * np.exp(-eta*w0*t) 

t   = np.linspace(0.0, 150.0, num=1e5)

#np.savetxt("analytic.dat", np.array([t, x01(t).real, x02(t).real]).T )  


# TRANSFER FUNCTION 
# http://isites.harvard.edu/fs/docs/icb.topic251677.files/notes23.pdf 
#
# x'' + 2b x' + w0**2 x = A cos wt  
#
# H(iw) = G(w) exp( -i phi(w) )  
#
# X = Xp + Xh. general + homogeneous  
#
 
w0 = np.sqrt(10.0) #np.sqrt(k/m)  
b  = 2.0/2 # c/2/m  
A  = 1.0 #D*1      
w  = 2.0 #w0*20  
 
G  = lambda w: 1.0/np.sqrt( (w0**2-w**2)**2 + (2*b*w)**2 ) 
P  = lambda w: 1.0/np.arctan( (w0**2 - w**2)/(2*b*w) )  
H  = lambda w: G(w) * np.exp( -P(w)*1.0j )
Xp = lambda w,t: np.real( H(w) * A * np.exp( w*t*1.0j ) )  


def f(U,t, _m, _c, _k, _F=0.0): 
    dUdt[0] =          U[1];                       
    dUdt[1] =  -_k/m * U[0] - _c/_m * U[1] + _F/_m; 

np.savetxt("analytic.dat", np.array([t, x02(t).real, Xp(w,t)]).T )  
 
