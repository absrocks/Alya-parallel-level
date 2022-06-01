import sys


def mu(t):
	return ( 1.4283e-6 * t**1.5 / (t + 110.5)  )

def rho(t,p=101325):
	return p/(t*287.0)


Pr 	= 0.71
L 	= 1
dT 	= 720
T0 	= 600
Ra 	= 1e6

mu0 = mu(T0)
rho0=rho(T0)
g = Ra/(Pr * dT * L**3 * rho0**2 ) * T0 * mu0**2

lambda_ref = mu(T0)*1005/0.71


verbose=False
if verbose:
	print('{:<30} {:>19.12f}'.format('rho0',rho0))  
	print('{:<30} {:>19.12f}'.format('g',g))
	print('{:<30} {:>19.12f}'.format('lambda(240)',mu(240)*1005/0.71))
	print('{:<30} {:>19.12f}'.format('lambda(960)',mu(960)*1005/0.71))
	print('{:<30} {:>19.12f}'.format('lambda_ref',lambda_ref))



#######################
# Nusselt
#######################

if len(sys.argv)>1:
	q = abs(float(sys.argv[1]))
	h  = q/(L*dT)
	nu = h*L/lambda_ref
	print('{:<6} {:>8.4f}'.format('Nu:',nu))

	if len(sys.argv)>2:
		q2 = abs(float(sys.argv[2]))
		h2  = q2/(L*dT)
		nu2 = h2*L/lambda_ref
		print('{:<6} {:>8.4f}'.format('Nu:',nu2))

		print('{:<6} {:>8.4f}'.format('Err:',abs(abs(q)-abs(q2))*2/(abs(q)+abs(q2))))












