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
	


#rewriteField('VELOCnew', [0.0,0.0])
rewriteField('CON01', [1.0])
rewriteField('CON02', [0.0])
rewriteField('TEMPE', [300.0])
rewriteField('CON03', [0.0456364948709])
rewriteField('CON04', [0.0])
rewriteH(642313.976105, -2031363.72727)

#--| H_f:                      -2031363.72727
#--| H_ox:                     642313.976105

