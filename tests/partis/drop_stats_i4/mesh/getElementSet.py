import sys

name=sys.argv[1]
with open('{}.dims.dat'.format(name),'r') as f:
    line = f.readline()
    line = f.readline()
    line = line.split()
    nelem = int(line[1])

with open('{}.set.elm'.format(name),'w') as f:
    for i in range(nelem):
        f.write('{} 1\n'.format(i+1))

nslave = 2
nelemslave = (nelem // nslave) + 1
with open('{}.partfie.dat'.format(name),'w') as f:
	j = 1
	for i in range(nelem):
		f.write('{} {}\n'.format(i+1,j))
		if (i+1) % nelemslave == 0:
			j += 1

