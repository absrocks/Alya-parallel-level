import paraview.simple as Simple 
"""
1) 
  Add Macro -> Add new macro
  Macros will be displayed in the macros menu and the macros toolbar
2)
  Select  vtkUnstructuredGrid or vtkPolyData 
3) 
  Check the directory, Alya files will be created 
4)
  Enjoy it...
  
  Created: 
        miguel zavala, Sep 21, 2015, Aachen, Germany 

"""


Script03 = """
import sys 

def __vtp2vtu__(vtp):
    appendFilter = vtk.vtkAppendFilter() 
    if vtk.VTK_MAJOR_VERSION <= 5:
      appendFilter.SetInput( vtp )
    else:
      appendFilter.SetInputData( vtp )
    appendFilter.Update()

    unstructuredGrid = vtk.vtkUnstructuredGrid()
    unstructuredGrid.ShallowCopy( appendFilter.GetOutput() )
    #unstructuredGrid.Update()

    vtk_n_pts   = unstructuredGrid.GetNumberOfPoints()
    vtk_n_cells = unstructuredGrid.GetNumberOfCells()
    print "[vtp2vtu] n_pts: %d, n_cells: %d " % (vtk_n_pts, vtk_n_cells) 

    return unstructuredGrid


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
    #
    # undefined into Alya 
    alya_types.append( ("POLYH",-42) )
    alya_types.append( ("VOXEL",-11) )

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
    #
    # undefined into Alya   
    vtk_types.append( ("VTK_POLYHEDRON",  42, vtk.vtkPolyhedron() ) )
    vtk_types.append( ("VTK_VOXEL",       11, vtk.vtkVoxel() ) )

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
    if(alya=="???"):  
      print vtk 
      exit()
    return alya
#--------------------------------------------------------------------------||--#    

def __save_data_coords__(GetOutput, fname=""): 
    vtk_data    = GetOutput
    vtk_n_pts   = vtk_data.GetNumberOfPoints()
    vtk_n_cells = vtk_data.GetNumberOfCells()
    if(vtk_n_pts==0): 
      print "[save_data_coords] ERROR: n_pts==0"
      sys.exit() 

    fcoords = open(fname+"_COORDINATES.alya", "w")
    print>> fcoords, "COORDINATES"
    for idx in range(vtk_n_pts):
      print>> fcoords, idx+1, 
      pt = vtk_data.GetPoint(idx)
      n_pt = len(pt)
      for j in range(n_pt): print>> fcoords, pt[j],
      print>> fcoords, ""

    print>> fcoords, "END_COORDINATES"
    fcoords.close() 
    print "|_Vtk_writer.save_data_cell file: \'%s\' " % (fcoords.name)


#n_VTK_POLYHEDRON = 0  

def __save_data_cell__(GetOutput, fname="", GlobalIdx=None):
    vtk_data    = GetOutput
    vtk_n_pts   = vtk_data.GetNumberOfPoints()
    vtk_n_cells = vtk_data.GetNumberOfCells()
    if(vtk_n_cells==0): 
      print "[save_data_cell] ERROR: n_cells==0"
      sys.exit() 

    CelDic = Vtk2Alya_cell()

    ftype = open(fname+"_TYPES.alya", "w")
    felem = open(fname+"_ELEMENTS.alya", "w")
    print>> ftype, "TYPES"
    print>> felem, "ELEMENTS"

    alya_cells_found = {}  
    for idx in range(vtk_n_cells):
      cell_ids = vtk.vtkIdList()
      vtk_data.GetCellPoints(idx, cell_ids) 
      print>> felem, idx+1,
      n_cell_ids = cell_ids.GetNumberOfIds()  

      # ELEMENTS
      for j in range(n_cell_ids):
        cell_id = cell_ids.GetId(j) 
        if(GlobalIdx!=None): 
          print>> felem, GlobalIdx[cell_id]+1,
        else:
          print>> felem, cell_id+1,
      print>> felem, ""

      # TYPE 
      cell_type = vtk_data.GetCellType(idx)
      alya_cel  = CelDic.get_alya_cell(cell_type)

      if(alya_cells_found.get(alya_cel)==None): 
        alya_cells_found[alya_cel]  = 1 
      else:
        alya_cells_found[alya_cel] += 1

      print>> ftype, idx+1, alya_cel  


      # FACEs 
      global n_VTK_POLYHEDRON 
      if(cell_type == vtk.VTK_POLYHEDRON ): 
        ptIds = vtk.vtkIdList()
        vtk_data.GetFaceStream(idx, ptIds)

        n_faces = ptIds.GetId(0)  
        print "cell:%d, n_faces: %d " % (idx, n_faces)
#        print "%d) cell:%d, n_faces: %d " % (n_VTK_POLYHEDRON, idx, n_faces)  

        n_faceStream  =  ptIds.GetNumberOfIds()   
        faceStream    = [ ptIds.GetId(i) for i in range(n_faceStream) ]  
        #print faceStream 
        n_faces       = faceStream[0]  
        faceStream[0] = 0  

        face_init = [] 
        face_size = [] 

        j = 1 
        i = 0
        while(i < n_faces):
          k = faceStream[j]
          face_init.append(j+1)
          face_size.append(k)
          #print i, j, k
          j = j + k + 2 - 1
          i += 1

        #print "[face_init]", face_init
        #print "[face_size]", face_size

        face_end = [ face_init[i] + face_size[i] for i in range(n_faces) ]
        #print "[face_end]", face_end 

        for a,b in zip(face_init, face_end):
          idx = range(a,b)
          #print idx, 
          #print faceStream[a:b]

        #n_VTK_POLYHEDRON += 1
        #global n_VTK_POLYHEDRON  


    print>> ftype, "END_TYPES"
    print>> felem, "END_ELEMENTS"
    ftype.close()
    felem.close() 


    flist = open(fname+"_LIST.alya", "w")
    for key, val in alya_cells_found.items(): 
      print>> flist, key, val 
    flist.close()

    print "|_Vtk_writer.save_data_cell file: \'%s\' \'%s\' " % (ftype.name, felem.name)


def __point_locator__( VTU01, VTU02 ):
  #locator = vtk.vtkPointLocator() 
  locator = vtk.vtkOctreePointLocator() 
  locator.SetDataSet( VTU01 )
  locator.BuildLocator()

  n_Pts = VTU02.GetNumberOfPoints()
  Pts   = VTU02.GetPointData()
  
  Cell_idx = []
  for i in range(n_Pts): 
    Pt    = VTU02.GetPoint(i)
    PtId  = locator.FindClosestPoint(Pt)
    Cell_idx.append( PtId )

  print "[point_locator]", len(Cell_idx)
  return Cell_idx


def __get_pts_data__(GetOutput):
    Data    = GetOutput.GetPointData() 
    n_props = Data.GetNumberOfArrays()
    for i in range(n_props):
      property_name = Data.GetArrayName(i)
      print " |_%d) \'%s\'" % (i+1, property_name)
      
      Prop   = Data.GetArray(property_name)
      n_cols = Prop.GetNumberOfComponents()
      n_rows = Prop.GetNumberOfTuples() 
      print "\'%s\':" % Prop.GetName(),   

      prop_range = Prop.GetRange()   
      print "[%f,%f]" %( prop_range[0], prop_range[1] )


      import datetime
      now = datetime.datetime.now()
      property_name += "_%d" % now.year
      property_name += "_%d" % now.month
      property_name += "_%d" % now.hour
      property_name += "_%d" % now.minute
      Fout = open(property_name+".dat", "w")

      k=0
      if(n_cols>1): 
        for j in range(n_rows): 
          print>> Fout, j+1, 
          for i in range(n_cols): 
            val = Prop.GetValue(k)
            print>> Fout, val, 
            k+=1
          print>> Fout, ""
    
      if(n_cols==1): 
        for j in range(n_rows): 
          print>> Fout, j+1, 
          for i in range(n_cols): 
            val = Prop.GetValue(k) 
            print>> Fout, val, 
            k+=1
          print>> Fout, ""


def __get_cell_data__(GetOutput, property_name=None):
    Data    = GetOutput.GetCellData()
    n_props = Data.GetNumberOfArrays()
    for i in range(n_props):
      print " |_%d) \'%s\'" % (i+1, Data.GetArrayName(i) )

    DATA = None
    if(property_name!=None):
      print "[get_cell_data]",
      Prop   = Data.GetArray(property_name)
      n_cols = Prop.GetNumberOfComponents()
      n_rows = Prop.GetNumberOfTuples()
      print "\'%s\':" % Prop.GetName(),

      prop_range = Prop.GetRange()
      print "[%f,%f]" %( prop_range[0], prop_range[1] )

      k=0
      DATA = []
      if(n_cols>1):
        for j in range(n_rows):
          data = []
          for i in range(n_cols):
            data.append( Prop.GetValue(k) )
            k+=1
          DATA.append( data )
   
      if(n_cols==1):
        for j in range(n_rows):
          for i in range(n_cols):
            DATA.append( Prop.GetValue(k) )
            k+=1
   
    return DATA


def __save_data_faces__(GetOutput):
    vtk_data    = GetOutput
    vtk_n_pts   = vtk_data.GetNumberOfPoints()
    vtk_n_cells = vtk_data.GetNumberOfCells()
    vtk_faces   = vtk_data.GetFaces() 
    #print vtk_faces.GetName(), vtk_faces.GetRange()
    #print vtk_faces.GetNumberOfTuples(), vtk_faces.GetNumberOfComponents()
    #
    #for i in range(vtk_faces.GetNumberOfTuples()): 
    #  print vtk_faces.GetValue(i),
    #print  



def __get_cell_type__(GetOutput, Array01):    
    vtk_data    = GetOutput
    vtk_n_cells = vtk_data.GetNumberOfCells()
    if(vtk_n_cells==0): 
      print "[save_data_cell] ERROR: n_cells==0"
      sys.exit() 

    #Array01 = vtk.vtkUnsignedCharArray()
    #Array01 = vtk.vtkDoubleArray()
    Array01.SetName("CELL_TYPE")
    Array01.SetNumberOfTuples(vtk_n_cells)

    TYPES = {} 
    for idx in range(vtk_n_cells):
      cell_type = vtk_data.GetCellType(idx)
      Array01.SetValue(idx, cell_type)
      
      key = cell_type
      if(key in TYPES): 
        TYPES[key] += 1
      else:
        TYPES[key]  = 0

    vtk_types = {}
    vtk_types[ 5] = "VTK_TRIANGLE"
    vtk_types[ 9] = "VTK_QUAD"
    vtk_types[10] = "VTK_TETRA"
    vtk_types[14] = "VTK_PYRAMID"
    vtk_types[13] = "VTK_WEDGE" 
    vtk_types[12] = "VTK_HEXAHEDRON"
    #
    # undefined into Alya   
    vtk_types[42] = "VTK_POLYHEDRON"
    vtk_types[11] = "VTK_VOXEL"


    for key, value in TYPES.iteritems():
      name = vtk_types.get(key, "???")
      print " +'%s'[%d] = %d" % (name, key, value)

    return Array01 

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def func01(A, Array):
  __save_data_coords__(A, "A") 
  #__save_data_cell__(A, "A") 
  __get_pts_data__(A) 
  #__get_cell_data__(A)
  __get_cell_type__( A, Array) 


#--------------------------------------------------------------------------||--#
# Get the two inputs
A = self.GetInputDataObject(0, 0)
B = self.GetInputDataObject(0, 1)

print "A:", A.GetClassName()
if B is not None:
  print "B:", B.GetClassName()


DATA = vtk.vtkDoubleArray()
if( A.IsA("vtkUnstructuredGrid") ): 
  func01(A, DATA)
elif( A.IsA("vtkPolyData") ): 
  VTU = __vtp2vtu__(A)
  func01(VTU, DATA)
else:
  print "[Alya] ERROR: \'%s\' != \'vtkUnstructuredGrid\' or \'vtkPolyData\' " % A.GetClassName()
  sys.exit()


if B is not None:
  if( not B.IsA("vtkPolyData") ): 
    print "[Alya] ERROR: \'%s\' != \'vtkPolyData\' " % B.GetClassName()
    sys.exit()
  else:
    VTU = __vtp2vtu__(B)
    GlobalIdx   = __point_locator__(A,  VTU)
    #__save_data_coords__(VTU, "B") 
    __save_data_cell__(VTU, "B") 


#
output = self.GetOutput()
output.ShallowCopy( A )
output.GetCellData().AddArray( DATA )
#
"""


PF = Simple.ProgrammableFilter()
PF.Script = Script03
Simple.Show(PF)


import os
print "Directory: \'%s\'" % os.getcwd()

print "OK!\n"
