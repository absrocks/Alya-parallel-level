import numpy as np


def ringLine(center,normal,r,dr,dn):
	#
	line = '    RING   '
	line += 'X={}, '.format(center[0])
	line += 'Y={}, '.format(center[1])
	line += 'Z={}, '.format(center[2])
	line += 'NX={}, '.format(normal[0])
	line += 'NY={}, '.format(normal[1])
	line += 'NZ={}, '.format(normal[2])
	line += 'R={}, '.format(r)
	line += 'DR={}, '.format(dr)
	line += 'DN={}'.format(dn)
	line += '\n'

	return line


r25cm 	= np.linspace( 22e-2, 28e-2,  7)
r50cm 	= np.linspace( 45e-2, 55e-2, 11)

print(r25cm)
print(r50cm)

nWit 	= len(r25cm)+len(r50cm)
shift 	= 0.0
dr	= 1e-2
dn	= 1e-1


with open('geomWitness.dat','w') as f:
	
	f.write('WITNESS, GEOMETRY NUMBER = {}\n'.format(nWit))
	#
	# 5 mm
	#
	for r in r25cm:
		f.write( ringLine([25e-2+shift,0,0],[1,0,0],r,dr,dn) )
	for r in r50cm:
		f.write( ringLine([50e-2+shift,0,0],[1,0,0],r,dr,dn) )

	f.write('END_WITNESS\n')


