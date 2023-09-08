import paraview.simple as Simple 
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
      for j in range(n_cell_ids):
        cell_id = cell_ids.GetId(j) 
        if(GlobalIdx!=None): 
          print>> felem, GlobalIdx[cell_id]+1,
        else:
          print>> felem, cell_id+1,
      print>> felem, ""

      cell_type = vtk_data.GetCellType(idx)
      alya_cel  = CelDic.get_alya_cell(cell_type)

      if(alya_cells_found.get(alya_cel)==None): 
        alya_cells_found[alya_cel]  = 1 
      else:
        alya_cells_found[alya_cel] += 1

      print>> ftype, idx+1, alya_cel  

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

#--------------------------------------------------------------------------||--#
def __get_cell_data__(GetOutput, property_name=None):
  Data    = GetOutput.GetCellData()
  n_props = Data.GetNumberOfArrays()
  for i in range(n_props):
    print " |_%d) \'%s\'" % (i+1, Data.GetArrayName(i) )

    property_name = Data.GetArrayName(i) 
    if(property_name!=None):
      print "[get_cell_data]"
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


#--------------------------------------------------------------------------||--#

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def func01(A):
  #__save_data_coords__(A, "A") 
  #__save_data_cell__(A, "A") 
  __get_pts_data__(A)
  __get_cell_data__(A, "A")

#--------------------------------------------------------------------------||--#
# Get the two inputs
A = self.GetInputDataObject(0, 0)
B = self.GetInputDataObject(0, 1)

print "A:", A.GetClassName()
if B is not None:
  print "B:", B.GetClassName()

if( A.IsA("vtkUnstructuredGrid") ): 
  func01(A)
elif( A.IsA("vtkPolyData") ): 
  VTU = __vtp2vtu__(A)
  func01(VTU)
else:
  print "[Alya] ERROR: \'%s\' != \'vtkUnstructuredGrid\' or \'vtkPolyData\' " % A.GetClassName()
  #sys.exit()
  func01(A)


if B is not None:
  if( not B.IsA("vtkPolyData") ): 
    print "[Alya] ERROR: \'%s\' != \'vtkPolyData\' " % B.GetClassName()
    sys.exit()
  else:
    VTU = __vtp2vtu__(B)
    GlobalIdx   = __point_locator__(A,  VTU)
    #__save_data_coords__(VTU, "B") 
    #__save_data_cell__(VTU, "B") 


#
#output = self.GetOutput()
#output.ShallowCopy(VTU_SURFACE)
#output.GetPointData().AddArray(labels)
#
"""


PF = Simple.ProgrammableFilter()
PF.Script = Script03
Simple.Show(PF)


import os
print "Directory: \'%s\'" % os.getcwd()

print "OK!\n"


