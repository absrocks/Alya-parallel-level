import sys
import os
meshName=sys.argv[1]
boufile='{}.fix.bou'.format(meshName)
realbou='{}.fix.bou.real'.format(meshName)


real_boundaries=[]
for i,arg in enumerate(sys.argv):
    if i > 1:
        real_boundaries.append(int(arg))

print('Real boundaries to leave in {}: {}'.format(realbou, real_boundaries))


f=open(boufile,'r')
g=open(realbou,'w')

for line in f:
    data = line.split()
    if len(data)>1:
        if int(data[1]) in real_boundaries:
            g.write(line)

f.close()
g.close()



