#!/usr/bin/python

import sys
import random
import math

def usage():
    print '+---------------------------------------------------------------+'
    print '|            Negative Jacobian element rotator                  |'
    print '|         for P1 and Q2 element (2D Tri and Quad)               |'
    print '|  Developed by Alfonso Santiago (alfonso.santiago@bsc.es)      |'
    print '|             Extended to P1 and Q2 elements                    |'
    print '|         by Matias Rivero (matias.rivero@bsc.es)               |'
    print '|                         | BSC |                               |'
    print '| Usage:                                                        |'
    print '|       ./command input output                                  |'
    print '|                                                               |'
    print '|To check that the software has worked OK just execute this     |'
    print '|script again to the modified file, and the message \"Some       |'
    print '|lines were modified!\" shouldn\'t appear.                        |'
    print '|                                                               |'
    print '|Ultra beta version.                                            |'
    print '+---------------------------------------------------------------+'
    sys.exit()


def calc_det(CoordNod,EleNod,elem):

    xjacm_f=[[0.0 for x in xrange(2)]for x in xrange(2)]
    deriv_f=[[0.0 for x in xrange(4)]for x in xrange(2)]

    deriv_f[0][0] = -0.25  
    deriv_f[0][1] = 0.25
    deriv_f[0][2] = 0.25
    deriv_f[0][3] = -0.25
    deriv_f[1][0] = -0.25
    deriv_f[1][1] = -0.25
    deriv_f[1][2] = 0.25
    deriv_f[1][3] = 0.25

    xjacm_f[0][0] = 0.0
    xjacm_f[0][1] = 0.0
    xjacm_f[1][0] = 0.0
    xjacm_f[1][1] = 0.0
    for m in range(4):
         xjacm_f[0][0] = xjacm_f[0][0] + CoordNod[EleNod[elem-1][m]-1][0]*deriv_f[0][m]
         xjacm_f[0][1] = xjacm_f[0][1] + CoordNod[EleNod[elem-1][m]-1][0]*deriv_f[1][m]
         xjacm_f[1][0] = xjacm_f[1][0] + CoordNod[EleNod[elem-1][m]-1][1]*deriv_f[0][m]
         xjacm_f[1][1] = xjacm_f[1][1] + CoordNod[EleNod[elem-1][m]-1][1]*deriv_f[1][m]

    gpdet_f = (xjacm_f[0][0]*xjacm_f[1][1]) - (xjacm_f[1][0]*xjacm_f[0][1])

    xc = (CoordNod[EleNod[elem-1][0]-1][0]+CoordNod[EleNod[elem-1][1]-1][0]+CoordNod[EleNod[elem-1][2]-1][0]+CoordNod[EleNod[elem-1][3]-1][0])/4.0
    yc = (CoordNod[EleNod[elem-1][0]-1][1]+CoordNod[EleNod[elem-1][1]-1][1]+CoordNod[EleNod[elem-1][2]-1][1]+CoordNod[EleNod[elem-1][3]-1][1])/4.0

    x1 = CoordNod[EleNod[elem-1][0]-1][0] - xc 
    y1 = CoordNod[EleNod[elem-1][0]-1][1] - yc 
    x2 = CoordNod[EleNod[elem-1][1]-1][0] - xc 
    y2 = CoordNod[EleNod[elem-1][1]-1][1] - yc 
    x3 = CoordNod[EleNod[elem-1][2]-1][0] - xc 
    y3 = CoordNod[EleNod[elem-1][2]-1][1] - yc 
    x4 = CoordNod[EleNod[elem-1][3]-1][0] - xc 
    y4 = CoordNod[EleNod[elem-1][3]-1][1] - yc 
  
    atan2_1 = math.atan2(y1,x1) 
    atan2_2 = math.atan2(y2,x2) 
    atan2_3 = math.atan2(y3,x3) 
    atan2_4 = math.atan2(y4,x4) 
   
    if(atan2_1 < 0.0):
         atan2_1 = atan2_1 + 2*math.pi
    if(atan2_2 < 0.0):
         atan2_2 = atan2_2 + 2*math.pi
    if(atan2_3 < 0.0):
         atan2_3 = atan2_3 + 2*math.pi
    if(atan2_4 < 0.0):
         atan2_4 = atan2_4 + 2*math.pi

    logical = ((atan2_2 > atan2_1) & (atan2_1 > atan2_4) & (atan2_4 > atan2_3))

    if(not(logical)):
         gpdet_f = 0.0 

    return gpdet_f


if(len(sys.argv)< 2):
    usage()
elif(sys.argv[1]=='-v'):
    usage()


f  = open(sys.argv[1], 'r')
fw = open(sys.argv[2], 'w')

FlagNodesPerElement=False
FlagElements=False
FlagCoordinates=False
FlagBoundaries=False
FlagSkewSystems=False
FlagModif=False

NumNodes=0
NodesElem=[]
CoordNod=[]
EleNod=[]

count_tri = 0
count_quad = 0

gpcar=[[0.0 for x in xrange(3)]for x in xrange(3)]

xjacm=[[0.0 for x in xrange(2)]for x in xrange(2)]
deriv=[[0.0 for x in xrange(4)]for x in xrange(2)]

deriv[0][0] = -0.25  
deriv[0][1] = 0.25
deriv[0][2] = 0.25
deriv[0][3] = -0.25
deriv[1][0] = -0.25
deriv[1][1] = -0.25
deriv[1][2] = 0.25
deriv[1][3] = 0.25


print 'Reading the original file...'
#--------------------------INICIO ESCRITURA DEL ARCHIVO-------------------------------
for line in f:
    line = line.strip()
    if (line.startswith('END_')):
        FlagNodesPerElement=False
        FlagElements=False
        FlagCoordinates=False
        FlagBoundaries=False
        FlagSkewSystems=False

    elif (line.startswith('NODES_PER') or FlagNodesPerElement):
        FlagNodesPerElement=True
        line=line.split()
        new=[]
        if line[0]!='NODES_PER_ELEMENT':
            new.append(int(line[1]))
            NodesElem.append(new)

    elif (line.startswith('ELEMENTS') or FlagElements):
        FlagElements=True
        line=line.split()
        new=[]
        if line[0]!='ELEMENTS':
            if NodesElem[int(line[0])-1][0]==4:
                new.append(int(line[1]))
                new.append(int(line[2]))
                new.append(int(line[3]))
                new.append(int(line[4]))
                count_quad = count_quad + 1
            if NodesElem[int(line[0])-1][0]==3:
                new.append(int(line[1]))
                new.append(int(line[2]))
                new.append(int(line[3]))
                count_tri = count_tri + 1
            EleNod.append(new)

    elif (line.startswith('COORDINATES') or FlagCoordinates):
        FlagCoordinates=True
        line=line.split()
        new=[]
        if line[0]!='COORDINATES':
            new.append(float(line[1]))
            new.append(float(line[2]))
            CoordNod.append(new)

    elif (line.startswith('BOUNDARIES') or FlagBoundaries):
        FlagBoundaries=True

    elif (line.startswith('SKEW_') or FlagSkewSystems):
        FlagSkewSystems=True

    else:
       print 'La linea leida no esta entre las opciones'



#----------------------------FIN LECTURA DEL ARCHIVO----------------------------------
f.close()
f  = open(sys.argv[1], 'r')

print 'Writing the new file...'
print 'In this mesh there are', count_tri, 'Triangles and', count_quad, 'Quadrangles'
#--------------------------INICIO ESCRITURA DEL ARCHIVO-------------------------------

FlagNodesPerElement=False
FlagElements=False
FlagCoordinates=False
FlagBoundaries=False
FlagSkewSystems=False

for line in f:
    line = line.strip()
    if (line.startswith('END_')):
        FlagNodesPerElement=False
        FlagElements=False
        FlagCoordinates=False
        FlagBoundaries=False
        FlagSkewSystems=False
        fw.write(line)
        fw.write('\n')

    elif (line.startswith('NODES_PER') or FlagNodesPerElement):
        FlagNodesPerElement=True
        fw.write(line)
        fw.write('\n')

    elif (line.startswith('ELEMENTS') or FlagElements):
        FlagElements=True
        line=line.split()

        if line[0]!='ELEMENTS':
            elem=int(line[0])
            if NodesElem[elem-1][0]==4:
                items = [0,1,2,3]
                while (calc_det(CoordNod,EleNod,elem) <= 0.0):
                    FlagModif=True
                    random.shuffle(items)   
                    tmp0 = EleNod[elem-1][items[0]]
                    tmp1 = EleNod[elem-1][items[1]]
                    tmp2 = EleNod[elem-1][items[2]]
                    tmp3 = EleNod[elem-1][items[3]]
                    EleNod[elem-1][0] = tmp0
                    EleNod[elem-1][1] = tmp1
                    EleNod[elem-1][2] = tmp2
                    EleNod[elem-1][3] = tmp3
                fw.write(line[0])
                fw.write(" ")
                fw.write(str(EleNod[elem-1][0]))
                fw.write(" ")
                fw.write(str(EleNod[elem-1][1]))
                fw.write(" ")
                fw.write(str(EleNod[elem-1][2]))
                fw.write(" ")
                fw.write(str(EleNod[elem-1][3]))
                fw.write('\n')
            if NodesElem[elem-1][0]==3:
                gpdet = (-CoordNod[EleNod[elem-1][0]-1][0]+CoordNod[EleNod[elem-1][1]-1][0])*(-CoordNod[EleNod[elem-1][0]-1][1]+CoordNod[EleNod[elem-1][2]-1][1]) \
                       -(-CoordNod[EleNod[elem-1][0]-1][1]+CoordNod[EleNod[elem-1][1]-1][1])*(-CoordNod[EleNod[elem-1][0]-1][0]+CoordNod[EleNod[elem-1][2]-1][0])
                if (gpdet < 0.0):
                    FlagModif=True                    
                    tmp1 = EleNod[elem-1][1]
                    tmp2 = EleNod[elem-1][2]
                    EleNod[elem-1][1] = tmp2 
                    EleNod[elem-1][2] = tmp1
                fw.write(line[0])
                fw.write(" ")
                fw.write(str(EleNod[elem-1][0]))
                fw.write(" ")
                fw.write(str(EleNod[elem-1][1]))
                fw.write(" ")
                fw.write(str(EleNod[elem-1][2]))
                fw.write('\n')

        else:
            fw.write(line[0])
            fw.write('\n')

    elif (line.startswith('COORDINATES') or FlagCoordinates):
        FlagCoordinates=True
        fw.write(line)
        fw.write('\n')
   
    elif (line.startswith('BOUNDARIES') or FlagBoundaries):
        FlagBoundaries=True
        fw.write(line)
        fw.write('\n')

    elif (line.startswith('SKEW_') or FlagSkewSystems):
        FlagSkewSystems=True
        fw.write(line)
        fw.write('\n')
    
    else:
       print 'La linea leida no esta entre las opciones'

if (FlagModif):
    print 'Some lines were modified!' 
else:
    print 'No line was modified, the mesh was OK.'

print 'Bye! :)'

f.close()
fw.close()


