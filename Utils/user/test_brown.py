#!/usr/bin/python

'''
======================================
This test allows to check if Brownian diffusivity
is correctly working in partis module in serie or
in parallel, which is specially critical due to 
diffculties generating good pseudorandom parallel
numbers.

Author: Edgar Olivares
Date: 15th November 2017

To validate results compare output with theoretical
result as shown next in gnuplot's example.

--------------
USING GNUPLOT FOR 2D CASE
--------------

D=5.35E-06
t=0.01
set xrange[0:0.0015]
plot 'test_brown.txt' u 1:($2/2500.0) w p,  1.0 - exp( -x**2/(4.0*D*t) )
======================================
'''

from math import sqrt
import sys

def usage():
    print '|| Name of the input is required without extension ||'
    print '|| ./command input                                 ||'
    sys.exit()

'''
This function opens the input file
and generates a dictionary with
  key:   # particle
  value: position
for future stats calculation.
'''
def readfile():
    fi = open(name+'.pts.res','r')
    for l in fi:
        if not l.startswith('#'):
            fields = l.split()            
            time   = float(fields[0]) # particle time
            pid    = int(fields[1])   # particle number
            x      = float(fields[2]) # x-coordinate 
            y      = float(fields[3]) # y-coordinate
            # Only take into account last time step
            if time == 0.01: 
                if not pid in d:
                    d[pid] = [x,y]
    fi.close()

'''
This fuction calculates
the dispersion of the particles
from its initial position
   (x,y) = (0,0)
to be able to compare with theory.
'''
def stats():
    o = open('test_brown.txt','w')
    for jj in range(1000):
        toler = jj/999.0*0.005
        ii = 0
        for key,coord in d.iteritems():
            xerro = 0
            for idime in range(2):
                xerro = xerro + (coord[idime]-0.005)**2
            xerro = sqrt(xerro)
            if xerro < toler:
                ii = ii + 1
        o.write(str(toler)+'\t'+str(ii)+'\n')
    o.close()

if __name__ == "__main__":
    if sys.argv[-1].split('/')[-1] == 'test_brown.py':
        usage()
    else:
        name = sys.argv[-1]
    d={}
    print '## --> Opening '+name+'.pts.res file'
    readfile()
    print '## --> Calculating statistics'
    stats()
    print '## --> Output test_brown.txt has been succesfuly created'
    print '## --> See you soon.'
