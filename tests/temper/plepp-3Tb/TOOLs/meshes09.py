# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v7.6.0 with dump python functionality
###
import numpy as np
F = lambda _n, _l, _L: int( np.floor( np.abs(1.0*_L/_l-0.0)*_n ) )

r06 =  200  
 
n01 =  64  
l01 =  0.5 
l02 =  1.0 
n02 = int( F(n01,l01,l02)*1.0 ) 
l03 =  0.45 
n03 = int( F(n01,l01,l03)*3.0 )  
#
l06 =  0.15
n06 = int( F(n01,l01,l06)*2.0 )
#
print "nX01, nX02, nYs:", n01, n02, n03
print "nY06", n06

r = lambda _S, _N: _S**(1.0/(_N-1.0)) #Scale Factor S 
L = lambda _L, _S, _N: _L*(1-r(_S,_N))/(1-r(_S,_N)**_N) 

r01    = 1.0/10
l01_01 = L(l01,r01,n01)  
r01_01 = r(    r01,n01) 
Segs01 = [ l01_01*r01_01**i for i in range(n01) ]   
print "Segs01", sum(Segs01) #, Segs01[-1]

r03    = 1.0/r01  
l03_01 = L(l02,r03,n02)
r03_01 = r(    r03,n02)
Segs03 = [ l03_01*r03_01**i for i in range(n02) ]
print "Segs03", sum(Segs03) #, Segs03[0]

r04    = 1.0/r03 
l04_01 = L(l02,r04,n02)
r04_01 = r(    r04,n02)
Segs04 = [ l04_01*r04_01**i for i in range(n02) ]
print "Segs04", sum(Segs04) #, Segs04 

r05    = 1.0/r04 
l05_01 = L(l02,r05,n02)
r05_01 = r(    r05,n02)
Segs05 = [ l05_01*r05_01**i for i in range(n02) ]
print "Segs05", sum(Segs05) #, Segs05 


import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, r'/home/jmake/z2016/BNDRYLAYER01/PYs')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
Vertex_1      = geompy.MakeVertex(0, 0, 0)
Vertex_2      = geompy.MakeVertex(0, 0.45, 0)
Line_1        = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
Extrusion_1   = geompy.MakePrismDXDYDZ(Line_1, 4, 0, 0)
Translation_1 = geompy.MakeTranslation(Line_1, 0.5, 0, 0)
Translation_2 = geompy.MakeTranslation(Line_1, 1.0, 0, 0)
Translation_3 = geompy.MakeTranslation(Line_1, 2.0, 0, 0)
Translation_4 = geompy.MakeTranslation(Line_1, 3.0, 0, 0)
Partition_1 = geompy.MakePartition([Extrusion_1], [Translation_1, Translation_2, Translation_3, Translation_4], [], [], geompy.ShapeType["FACE"], 0, [], 0)

All = geompy.CreateGroup(Partition_1, geompy.ShapeType["FACE"])
geompy.UnionIDs(All, [2, 19, 26, 12, 33])

Ys = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Ys, [4, 23, 30, 9, 16, 37])
X01 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(X01, [11, 7])
X02 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(X02, [14, 18])
X03 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(X03, [25, 21])
X04 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(X04, [28, 32])
X05 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(X05, [35, 39])
#
Inlet = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Inlet, [4])
Outlet = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Outlet, [37])
Top01 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Top01, [14, 21, 28, 35, 7])
Bottom01 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Bottom01, [11])
Bottom02 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Bottom02, [18])
Bottom03 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Bottom03, [25, 32])
Bottom04 = geompy.CreateGroup(Partition_1, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Bottom04, [39])
# 
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( Partition_1, 'Partition_1' )
#
geompy.addToStudyInFather( Partition_1, All, 'All' )
geompy.addToStudyInFather( Partition_1, Ys, 'Ys' )
geompy.addToStudyInFather( Partition_1, X01, 'X01' )
geompy.addToStudyInFather( Partition_1, X02, 'X02' )
geompy.addToStudyInFather( Partition_1, X03, 'X03' )
geompy.addToStudyInFather( Partition_1, X04, 'X04' )
geompy.addToStudyInFather( Partition_1, X05, 'X05' )
#
geompy.addToStudyInFather( Partition_1, Inlet, 'Inlet' )
geompy.addToStudyInFather( Partition_1, Outlet, 'Outlet' )
geompy.addToStudyInFather( Partition_1, Top01, 'Top01' )
geompy.addToStudyInFather( Partition_1, Bottom01, 'Bottom01' )
geompy.addToStudyInFather( Partition_1, Bottom02, 'Bottom02' )
geompy.addToStudyInFather( Partition_1, Bottom03, 'Bottom03' )
geompy.addToStudyInFather( Partition_1, Bottom04, 'Bottom04' )
#
# New mesh 
Face_1 = geompy.MakeFaceHW(2, l06, 1)
geompy.TranslateDXDYDZ(Face_1, 2, -l06/2, 0)
Vertex_2 = geompy.MakeVertex(2, 0, 0)
Vertex_1 = geompy.MakeVertex(2, -l06, 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_2)
#
Partition_2 = geompy.MakePartition([Face_1], [Line_1], [], [], geompy.ShapeType["FACE"], 0, [], 0)
X03_02 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(X03_02, [7, 11])
X04_02 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(X04_02, [18, 14])
Inlet02 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Inlet02, [4])
Outlet02 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Outlet02, [16])
All02 = geompy.CreateGroup(Partition_2, geompy.ShapeType["FACE"])
geompy.UnionIDs(All02, [12, 2])
Top02 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Top02, [14, 7])
Bottom05 = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Bottom05, [11, 18])
Yes = geompy.CreateGroup(Partition_2, geompy.ShapeType["EDGE"])
geompy.UnionIDs(Yes, [16, 4, 9])
#
geompy.addToStudy( Partition_2, 'Partition_2' )
geompy.addToStudyInFather( Partition_2, X03_02, 'X03_02' )
geompy.addToStudyInFather( Partition_2, X04_02, 'X04_02' )
geompy.addToStudyInFather( Partition_2, Inlet02, 'Inlet02' )
geompy.addToStudyInFather( Partition_2, Outlet02, 'Outlet02' )
geompy.addToStudyInFather( Partition_2, All02, 'All02' )
geompy.addToStudyInFather( Partition_2, Top02, 'Top02' )
geompy.addToStudyInFather( Partition_2, Bottom05, 'Bottom05' )
geompy.addToStudyInFather( Partition_2, Yes, 'Yes' )


###
### SMESH component
###

import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

smesh = smeshBuilder.New(theStudy)
Mesh08 = smesh.Mesh(Partition_1)
Regular_1D = Mesh08.Segment()
NoGeneral = Regular_1D.NumberOfSegments( 2 )
NoGeneral.SetDistrType( 0 )
Quadrangle_2D = Mesh08.Quadrangle(algo=smeshBuilder.QUADRANGLE)
#
NoX01 = smesh.CreateHypothesis('NumberOfSegments')
NoX01.SetNumberOfSegments( n01 )
NoX01.SetDistrType( 1 )
NoX01.SetScaleFactor( r01 )
NoX01.SetReversedEdges( [] )
NoX01.SetObjectEntry( "0:1:1:10" )
status = Mesh08.AddHypothesis(NoX01,X01)
#
Start_and_End_Length_1 = smesh.CreateHypothesis('StartEndLength')
Start_and_End_Length_1.SetStartLength( Segs01[-1] )
Start_and_End_Length_1.SetEndLength(   Segs03[ 0] )
Start_and_End_Length_1.SetReversedEdges( [] )
Start_and_End_Length_1.SetObjectEntry( 'Partition_1' )
status = Mesh08.AddHypothesis(Start_and_End_Length_1,X02)
#
NoX03 = smesh.CreateHypothesis('NumberOfSegments')
NoX03.SetNumberOfSegments( n02 )
NoX03.SetDistrType( 1 )
NoX03.SetScaleFactor( r03 )
NoX03.SetReversedEdges( [] )
NoX03.SetObjectEntry( "0:1:1:10" )
status = Mesh08.AddHypothesis(NoX03,X03)
#
NoX04 = smesh.CreateHypothesis('NumberOfSegments')
NoX04.SetNumberOfSegments( n02 )
NoX04.SetDistrType( 1 )
NoX04.SetScaleFactor( r04 )
NoX04.SetReversedEdges( [] )
NoX04.SetObjectEntry( "0:1:1:10" )
status = Mesh08.AddHypothesis(NoX04,X04)
#
NoX05 = smesh.CreateHypothesis('NumberOfSegments')
NoX05.SetNumberOfSegments( n02 )
NoX05.SetDistrType( 1 )
NoX05.SetScaleFactor( r05 )
NoX05.SetReversedEdges( [] )
NoX05.SetObjectEntry( "0:1:1:10" )
status = Mesh08.AddHypothesis(NoX05,X05)
#
NoYs = smesh.CreateHypothesis('NumberOfSegments')
NoYs.SetNumberOfSegments( n03 )
NoYs.SetScaleFactor( r06 )
NoYs.SetReversedEdges( [] )
NoYs.SetObjectEntry( "0:1:1:30" )
status = Mesh08.AddHypothesis(NoYs,Ys)
#
isDone = Mesh08.Compute()
isDone = Mesh08.SplitQuadObject( Mesh08, 1 )
#
Inlet_1 = Mesh08.GroupOnGeom(Inlet,'Inlet',SMESH.EDGE)
Outlet_1 = Mesh08.GroupOnGeom(Outlet,'Outlet',SMESH.EDGE)
Top01_1 = Mesh08.GroupOnGeom(Top01,'Top01',SMESH.EDGE)
Bottom01_1 = Mesh08.GroupOnGeom(Bottom01,'Bottom01',SMESH.EDGE)
Bottom02_1 = Mesh08.GroupOnGeom(Bottom02,'Bottom02',SMESH.EDGE)
Bottom03_1 = Mesh08.GroupOnGeom(Bottom03,'Bottom03',SMESH.EDGE)
Bottom04_1 = Mesh08.GroupOnGeom(Bottom04,'Bottom04',SMESH.EDGE)
All_1 = Mesh08.GroupOnGeom(All,'All',SMESH.FACE)
#
smesh.SetName(Inlet_1, 'Inlet')
smesh.SetName(Top01_1, 'Top01')
smesh.SetName(Outlet_1, 'Outlet')
smesh.SetName(Bottom02_1, 'Bottom02')
smesh.SetName(Bottom01_1, 'Bottom01')
smesh.SetName(Bottom04_1, 'Bottom04')
smesh.SetName(Bottom03_1, 'Bottom03')
smesh.SetName(All_1, 'All')

#
[ X01_1, X02_1, X03_1, X04_1, X05_1, Ys_1 ] = Mesh08.GetMesh().GetSubMeshes()
X01_1 = Mesh08.GetSubMesh( X01, 'X01' )
X02_1 = Mesh08.GetSubMesh( X02, 'X02' )
X03_1 = Mesh08.GetSubMesh( X03, 'X03' )
X04_1 = Mesh08.GetSubMesh( X04, 'X04' )
X05_1 = Mesh08.GetSubMesh( X05, 'X05' )
Ys_1  = Mesh08.GetSubMesh( Ys, 'Ys' )

## Set names of Mesh objects
smesh.SetName(Regular_1D.GetAlgorithm(), 'Regular_1D')
smesh.SetName(Quadrangle_2D.GetAlgorithm(), 'Quadrangle_2D')
smesh.SetName(NoX01, 'NoX01')
smesh.SetName(Start_and_End_Length_1, 'Start and End Length_1')
smesh.SetName(NoGeneral, 'NoGeneral')
smesh.SetName(NoX05, 'NoX05')
smesh.SetName(NoX03, 'NoX03')
smesh.SetName(NoX04, 'NoX04')
smesh.SetName(Mesh08.GetMesh(), 'Mesh08')
smesh.SetName(X03_1, 'X03')
smesh.SetName(X02_1, 'X02')
smesh.SetName(X01_1, 'X01')
smesh.SetName(X05_1, 'X05')
smesh.SetName(X04_1, 'X04')
smesh.SetName(Ys_1, 'Ys')

NoAll02 = smesh.CreateHypothesis('NumberOfSegments')
NoAll02.SetNumberOfSegments( n06 )
NoAll02.SetDistrType( 0 )
#
Mesh09 = smesh.Mesh(Partition_2)
status = Mesh09.AddHypothesis(NoAll02)
status = Mesh09.AddHypothesis(Regular_1D)
status = Mesh09.AddHypothesis(Quadrangle_2D)

status = Mesh09.AddHypothesis(NoX03,X03_02)
status = Mesh09.AddHypothesis(NoX04,X04_02)

NoYes = smesh.CreateHypothesis('NumberOfSegments')
NoYes.SetNumberOfSegments( n06 )
NoYes.SetDistrType( 1 )
NoYes.SetScaleFactor( 1.0/r06 )
NoYes.SetReversedEdges( [] )
NoYes.SetObjectEntry( "0:1:1:6" )
status = Mesh09.AddHypothesis(NoYes,Yes)

isDone = Mesh09.Compute()
isDone = Mesh09.SplitQuadObject( Mesh09, 1 )

smesh.SetName(Bottom05, 'Bottom05')
smesh.SetName(Mesh09.GetMesh(), 'Mesh09')

Inlet_2 = Mesh09.GroupOnGeom(Inlet02,'Inlet02',SMESH.EDGE)
Outlet02_1 = Mesh09.GroupOnGeom(Outlet02,'Outlet',SMESH.EDGE)
Top02_1 = Mesh09.GroupOnGeom(Top02,'Top02',SMESH.EDGE)
Bottom05_1 = Mesh09.GroupOnGeom(Bottom05,'Bottom05',SMESH.EDGE)
All02_1 = Mesh09.GroupOnGeom(All02,'All02',SMESH.FACE)
#
smesh.SetName(Inlet_2, 'Inlet_2')
smesh.SetName(Outlet02_1, 'Outlet02_1')
smesh.SetName(Top02_1, 'Top02_1')
smesh.SetName(Bottom05_1, 'Bottom05_1')
smesh.SetName(All02_1, 'All02_1')


### 
if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)
