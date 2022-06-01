import sys

perName = sys.argv[1]

fPer = open(perName,'r')
binding = {}
for line in fPer:
	data = line.split()
	if len(data) >1:
		mast = int(data[0])
		targ = int(data[1])
		if targ in binding.keys():
			binding[mast] = targ
			#print('----| m:{} t:{}'.format(targ,mast))
		else:
			binding[targ] = mast
			#print('----| m:{} t:{}'.format(mast,targ))
fPer.close()



for k in sorted(binding.keys()):
	mast = binding[k]
	targ = k
	if mast in binding.keys():
		shihan = binding[mast]
		for k2 in binding.keys():
			if binding[k2] == mast:
				binding[k2] = shihan
			


fPer = open(perName,'w')
for k in sorted(binding.keys()):
	fPer.write('{} {}\n'.format(binding[k],k))
fPer.close()

