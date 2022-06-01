import os

def rewriteField(name,val):
	oldname = '{}.alya'.format(name)
	newname = '{}.new.alya'.format(name)
	
	fOld = open(oldname,'r')
	fNew = open(newname,'w')
	
	line = True
	
	while line:
		line = fOld.readline()
		data = line.split()
		if len(data) > 1:
			ip  = int(data[0])
			
			output = '{}'.format(ip)
			for v in val:
				output += ' {}'.format(v)
			output += '\n'

			fNew.write(output)
	
	fOld.close()
	fNew.close()
	
	command = 'mv {0} {1}'.format(newname, oldname)
	os.system(command)


def rewriteH(HZ0, HZ1):
	oldname = 'CON03.alya'
	newname = 'ENTHA.alya'
	
	fOld = open(oldname,'r')
	fNew = open(newname,'w')
	
	line = True
	
	while line:
		line = fOld.readline()
		data = line.split()
		if len(data) > 1:
			ip  = int(data[0])
			Z   = float(data[1])

			output = '{}'.format(ip)
			output += ' {}'.format( HZ0 + (HZ1-HZ0) * Z )
			output += '\n'

			fNew.write(output)
	
	fOld.close()
	fNew.close()
	

def rewriteYc(Zst, Ycst):
	zname = 'CON03.alya'
	ycname = 'CON01.alya'
	
	fz = open(zname,'r')
	fyc = open(ycname,'w')
	
	line = True
	
	while line:
		line = fz.readline()
		data = line.split()
		if len(data) > 1:
			ip  = int(data[0])
			Z   = float(data[1])

			if Z <= Zst:
			    eta = Z/Zst
			    yc = eta * Ycst
			else:
			    eta = (1.0-Z)/(1.0-Zst)
			    yc = eta * Ycst

			output = '{}'.format(ip)
			output += ' {}'.format( yc )
			output += '\n'

			fyc.write(output)
	
	fz.close()
	fyc.close()


rewriteYc(0.055,0.025)


