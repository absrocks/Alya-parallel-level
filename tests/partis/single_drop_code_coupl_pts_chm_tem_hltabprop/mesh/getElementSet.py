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

