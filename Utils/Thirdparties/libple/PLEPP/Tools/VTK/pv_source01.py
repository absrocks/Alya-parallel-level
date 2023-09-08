
#--------------------------------------------------------------------------||--#
"""
1) 
  Add Macro -> Add new macro
  Macros will be displayed in the macros menu and the macros toolbar
2)
  Select  vtkUnstructuredGrid or vtkPolyData 
3) 
  Check the directory, Alya files will be created 
  
Created: 
        miguel zavala, Oct 29, 2014, Barcelona

"""
#--------------------------------------------------------------------------||--#

Script03 = """
import numpy as np 
import vtk 
#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def Read_alya_geo( fname, INIT_KEY=''):
    #INIT_KEY = INIT_KEY.replace("_", "")
    #INIT_KEY = INIT_KEY.replace("-", "")
    #INIT_KEY = INIT_KEY.replace("&", "")
    END_KEY  = "END_" + INIT_KEY

    data  = open(fname, "r")
    lines = data.readlines()
    data.close()

    nline = len(lines)    

    IDs = { INIT_KEY:-1, END_KEY:-1 }
    for i in range(nline):
      ok   = False
      line = lines[i]
      if( line.find(INIT_KEY) >= 0 ):
        line    = line.replace(',', ' ')
        splited = line.split()
        n_split = len(splited)

        if(n_split == 1): 
          ok = True 
        if(n_split > 1): 
          print 
          print "WARNING: \'len( line.split() ) > 1\' ->" 
          print "                                       ", splited
          
          if(splited[1] == "NOTSORTED"):
            ok = True 
          print 

        if( ok ): 
          line = splited[0]
          if( line==INIT_KEY ): IDs[INIT_KEY] = i+1
          if( line==END_KEY  ): IDs[END_KEY ] = i+0

    print "IDs:", IDs, 

    n_XYZ = IDs[END_KEY]-IDs[INIT_KEY] 
    XYZ   = [ [] for i in range(n_XYZ) ]
    for i in range(IDs[INIT_KEY], IDs[END_KEY]): 
      line = lines[i]
      line = line.strip()
      line = line.split()
      idx  = eval(line[0]) - 1  # C style 
      rows  = [eval(val) for val in line[1:]]
      XYZ[idx] = rows 

    print ", n_lines", i 

    return np.array(XYZ) 
#--------------------------------------------------------------------------||--#


#--------------------------------------------------------------------------||--#
def Alya_cell_types():
    # col=1: kernel/domain/elmtyp.f90        <<==
    # col=2: kernel/elsest/elsest_geogid.f90 <<=
    alya_types = []
    alya_types.append( ("TRI03", 10) )
    alya_types.append( ("QUA04", 12) )
    alya_types.append( ("TET04", 30) )
    alya_types.append( ("PYR05", 32) )
    alya_types.append( ("PEN06", 34) )
    alya_types.append( ("HEX08", 37) )

    return alya_types

    #self.salome2vtk = {}
    #for salome_type, vtk_type in zip(self.salome_types, self.VTK.vtk_types):
    #  self.salome2vtk[salome_type] = vtk_type

def Vtk_cell_types():
    vtk_types = []    #(vtk_type, vtk_id, vtkCell) <= linear cell types found in VTK
    vtk_types.append( ("VTK_TRIANGLE",     5, vtk.vtkTriangle()   ) )
    vtk_types.append( ("VTK_QUAD",         9, vtk.vtkQuad()       ) )
    vtk_types.append( ("VTK_TETRA",       10, vtk.vtkTetra()      ) )
    vtk_types.append( ("VTK_PYRAMID",     14, vtk.vtkPyramid() ) )
    vtk_types.append( ("VTK_WEDGE",       13, vtk.vtkWedge()       ) )
    vtk_types.append( ("VTK_HEXAHEDRON",  12, vtk.vtkHexahedron() ) )

    return vtk_types


class Vtk2Alya_cell():
  def __init__(self):
    alya_types = Alya_cell_types()
    vtk_types  = Vtk_cell_types()

    self.alya2vtk = {}
    self.vtk2alya = {}
    for alya_type, vtk_type in zip(alya_types, vtk_types):
      alya_idx = alya_type[1]
      self.alya2vtk[ alya_idx ] = vtk_type

      vtk_idx  = vtk_type[1]
      self.vtk2alya[ vtk_idx ]  = alya_idx


  def get_vtk_cell(self, alya):
    vtk = self.alya2vtk.get(alya, "???")
    return vtk


  def get_alya_cell(self, vtk):
    alya = self.vtk2alya.get(vtk, "???")
    return alya


#--------------------------------------------------------------------------||--#
class VTU_CREATOR:
  def __init__(self):
    self.obj       = None
    self.vtk_pts   = vtk.vtkPoints()
    self.vtk_cells = vtk.vtkCellArray()

    self.vtk_props       = []
    self.vtk_cells_type  = []
    self.vtk_cells_id    = []

    self.n_pts   = 0
    self.n_cells = 0
    print "+Vtk_writer"


  def set_point(self, coord):
    self.vtk_pts.InsertNextPoint( coord[0], coord[1], coord[2] )
    self.n_pts += 1


  def set_cell(self, nodes, vtk_cell_type):
    self.__set_cell__(nodes, vtk_cell_type)


  def set_tetra(self, nodes):
    self.__set_cell__(nodes, vtk.vtkTetra())


  def set_vtk_unstructured(self, vtk_data=vtk.vtkUnstructuredGrid() ):
    vtk_data.SetPoints( self.vtk_pts )

    vtk_data.Allocate( self.n_cells, 1) #<-- paraview ?? 
    for cell_type, cell_ids in zip(self.vtk_cells_type, self.vtk_cells_id):
      vtk_data.InsertNextCell(cell_type, cell_ids)

    print "|_Vtk_writer.to_vtkunstructured_grid: ",
    print "pts/cells: %d %d" %( self.n_pts, self.n_cells )


  def __set_cell__(self, nodes, vtk_cell):
    i = 0
    vtk_list = vtk.vtkIdList()
    for node in nodes:
      if (node<0):
        print "Negative node!!"
        exit(1)

      vtk_list.InsertNextId(node)
      i += 1

    self.vtk_cells_id.append(   vtk_list               )
    self.vtk_cells_type.append( vtk_cell.GetCellType() )

    self.n_cells += 1
#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
import glob 
import sys
ALYA_KEYS_LIST = ["COORDINATES", "ELEMENTS"]

ALYA_PATH_FILE = "ALYA.txt"
ALYA_PATH_LIST = glob.glob(ALYA_PATH_FILE)
#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#



#--------------------------------------------------------------------------||--#
import argparse
parser = argparse.ArgumentParser( fromfile_prefix_chars='@' )
parser.add_argument('-P', action='store', dest='ALYA_PATH',
                    default='./', type=str, nargs=1,
                    help='ALYA_PATH',
                    )

parser.add_argument('-N', action='store', dest='ALYA_NAME',
                    default='---', type=str, nargs=1,
                    help='ALYA_NAME',
                    )

if( ALYA_PATH_LIST==[] ): 
  print "ERROR: \'%s\' doesnt find!!" % ALYA_PATH_FILE
  parser.print_help()
  print 
  print 
  sys.exit()    

Results = parser.parse_args( ['@%s' % ALYA_PATH_FILE] )

ALYA_PATH = Results.ALYA_PATH[0]
ALYA_NAME = Results.ALYA_NAME[0]

ALYA_CASE = "%s/%s" % (ALYA_PATH.strip(), ALYA_NAME.strip()) 
print "ALYA_CASE: \'%s\' " % ALYA_CASE


ALYA_FILES_LIST = glob.glob(ALYA_CASE+"*")
#print ALYA_FILES_LIST 

import re 

ALYA_GEO_FILES = {}
for f in ALYA_FILES_LIST:
    for line in open(f).readlines():
      for key in ALYA_KEYS_LIST: 
        #if re.match("END_%s" % key, line):
        if( line.find("END_%s" % key) >= 0 ): 
            if( ALYA_GEO_FILES.has_key( key ) ):
              ALYA_GEO_FILES.append( f )
            else: 
              ALYA_GEO_FILES[key] = [ f ] 

print "ALYA_GEO_FILES:", ALYA_GEO_FILES


ALYA_DATA = { }
for k,v in ALYA_GEO_FILES.iteritems(): 
  #print "ALYA_GEO_FILE:", k, v[0]
  if( not ALYA_DATA.has_key( k ) ):
    ALYA_DATA[ k ] = np.zeros(0)
    ALYA_DATA[ k ] = Read_alya_geo( fname=v[0], INIT_KEY=k)
    print ALYA_DATA[ k ].shape 

#print ALYA_DATA.keys()
#--------------------------------------------------------------------------||--#

#--------------------------------------------------------------------------||--#
VTU = VTU_CREATOR( )

for coord in ALYA_DATA['COORDINATES']: 
  VTU.set_point( coord ) 

VTK2ALYA = Vtk2Alya_cell()
for cell in ALYA_DATA['ELEMENTS']: 
  alya_cell = 30 
  vtk_cell  = VTK2ALYA.get_vtk_cell( alya_cell ) 

  cell -= 1 
  VTU.set_cell( cell, vtk_cell[2] )

VTU.set_vtk_unstructured( vtk_data=self.GetUnstructuredGridOutput() ) 
#--------------------------------------------------------------------------||--#
"""

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
import paraview.simple as Simple 


#-----------------------------------------------------| Programmable Source |--#
programmableSource1                          = Simple.ProgrammableSource()
programmableSource1.OutputDataSetType        = 'vtkUnstructuredGrid'
programmableSource1.Script                   = Script03 
programmableSource1.ScriptRequestInformation = ''  
programmableSource1.PythonPath               = ''
Simple.Show(programmableSource1) 



#--------------------------------------------------------------------------||--#
import os
print "Directory: \'%s\'" % os.getcwd()

print "OK!\n"


#--------------------------------------------------------------------------||--#
#
# http://www.paraview.org/Wiki/Python_Programmable_Filter 
#
#--------------------------------------------------------------------------||--#