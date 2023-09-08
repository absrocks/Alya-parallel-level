import paraview.simple as Simple 
#------------------------------------------------------------------------------#
"""
1) 
  Add Macro -> Add new macro
  Macros will be displayed in the macros menu and the macros toolbar
2)
  Select  vtkUnstructuredGrid or vtkPolyData 
3) 
  Check the directory, Alya files will be created 

Other:
  paraview --script=pv_lodi.py

Created: 
        miguel zavala, Nov 01, 2014, Barcelona
"""

#------------------------------------------------------------------------------#
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


def __point_locator__xxx( VTU01, VTU02 ):
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


def __point_locator__( A, B ):
  locator = vtk.vtkPointLocator()
  locator.InitPointInsertion(vtk.vtkPoints(), A.GetBounds())
  for i in range(A.GetNumberOfPoints()):
    locator.InsertNextPoint(A.GetPoints().GetPoint(i))

  ## Label the points in B that are also in A
  Cell_idx = []
  for i in range(B.GetNumberOfPoints()):
    point = B.GetPoints().GetPoint(i)
    if locator.IsInsertedPoint(point) != -1:
      Cell_idx.append(i)

  print "[point_locator]", len(Cell_idx)
  return Cell_idx


def __get_pts_data__(GetOutput, property_name=None):
    Data    = GetOutput.GetPointData() 
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

def __save_raw__(Data, fname, HEADER="", BOOTOM=""):
  Fout = open(fname, "w")
  print>> Fout, HEADER
  
  n_Data = len(Data)
  for i in range(n_Data):
    print>> Fout, i+1, 
    line   = Data[i]
    n_line =  1
    if(isinstance(line, list)): 
      n_line = len(line)
      for j in range(n_line):
        print>> Fout, line[j], 
    else:
      print>> Fout, line, 
    print>> Fout

  print>> Fout, BOOTOM

  Fout.close()


def __set_points_prop_scalar__(R, n_pts, prop_name=""):    
    array = vtk.vtkDoubleArray()
    array.SetName(prop_name)
    array.SetNumberOfComponents(1)
    array.SetNumberOfTuples( n_pts )
    for i in range(n_pts): array.SetValue(i, R[i])

    return array


def __set_points_prop_vector__(R, n_pts, prop_name=""):
    array = vtk.vtkDoubleArray()
    array.SetName(prop_name)
    array.SetNumberOfComponents(3)
    array.SetNumberOfTuples( n_pts )

    j = 0
    for i in range(n_pts):
      for k in range(3): array.SetValue(j+k, R[i,k])
      j += 3

    return array 


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __get_coords__(GetOutput):
    vtk_data    = GetOutput
    vtk_n_pts   = vtk_data.GetNumberOfPoints()
    vtk_n_cells = vtk_data.GetNumberOfCells()
    if(vtk_n_pts==0):
      print "[save_data_coords] ERROR: n_pts==0"
      sys.exit()

    Pts = []
    for idx in range(vtk_n_pts):
      pt = vtk_data.GetPoint(idx)
      Pts.append( pt )

    print "[get_coords] n_pts:", vtk_n_pts

    return Pts

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def Func01( VTU, Props = [] ):   
  import numpy as np 
  import sys 
  sys.path.append("/home/jmake/z2014/Alya/Runner/Couplings/VORTEX2D01/PYs")

 #-------------------------#
  Data = []

  Pts = __get_coords__( VTU ) 
  Pts =  np.array( Pts )

  X = Pts[:,0]
  Y = Pts[:,1]
  Z = Pts[:,2]

 #-------------------------#
  L       = np.max(Z) - np.min(Z)

  R       = np.sqrt(X**2 + Y**2)
  a       = np.max(R)
  R       = 1.0 - (R/a)**2
  D       = 2*a  # [m]

  T0     = 273.0  # [K]
  rho    = 1.2      # [kg/m3] 
  mu    = 2e-5    # [Pa.s]
  Re     = 2e3      # [m/s]        #<-- laminar 
  Vmax= Re*mu/rho/D

  #-------------------------#
  Re = (L*0.25)/0.06
  Vc = Re*mu/D/rho 
  g   = 4.0*mu*Vc/rho/(D/2.0)**2
  Vin=np.sqrt(  Vc**2 - 2*g*(L*1.00)  )
  
  print "L:  ",   L, "[m]"
  print "g:  ",   g, "[m/s2]" 
  #print "Vc:", Vc, "[m/s], Vmax:", Vmax, "Vin:", Vin
  print "Vin:    ", Vin, 
  print "Vinf:   ", Vc, 
  print "Vmax:", Vmax

  #-------------------------#
  Vel       = np.zeros( Pts.shape     ) 
  Vel[:,2] = R*Vc

  T          = np.zeros( Pts.shape[0] ) 
  T          = np.where( Z<=L*0.50,  T0, 3*T0) 

  #-------------------------#
  #Data.append( __set_points_prop_scalar__(R, R.shape[0], "Rxy") )
  Data.append( __set_points_prop_scalar__(T   ,     T.shape[0], "TEMPE") )
  Data.append( __set_points_prop_vector__(Vel, Vel.shape[0], "VELOC") )

  return Data 

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
A = self.GetInputDataObject(0, 0)
print "A:", A.GetClassName()

Arrays = [] 
if( A.IsA("vtkUnstructuredGrid") ): 
  Arrays = Func01( A ) 
elif( A.IsA("vtkPolyData") ): 
  VTU = __vtp2vtu__(A)
  Arrays = Func01( VTU ) 
else:
  print "[Alya] ERROR: \'%s\' != \'vtkUnstructuredGrid\' or \'vtkPolyData\' " % A.GetClassName()
  sys.exit()

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
output = self.GetOutput()
if( output.IsA("vtkUnstructuredGrid") or output.IsA("vtkPolyData") ): 

  Things = ["IDx"]
  Props  = []
  n_Props = output.GetPointData().GetNumberOfArrays()
  for i in range(n_Props):
    prop = output.GetPointData().GetArrayName(i)
    Props.append(prop)

  for think in Things:
    if(think in Props): Props.remove(think)

  for prop in Props:
    output.GetPointData().RemoveArray( prop )
    print "Deleted: '%s'" % prop

  output.ShallowCopy(A)
  for array in Arrays: output.GetPointData().AddArray(array)

else:
  print "[ERROR] Output Data Set Type must be chosen!!"
#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
"""

#------------------------------------------------------------------------------#
#
#VTU = Simple.OpenDataFile("../MESHES/Mesh01.vtu")
#Simple.Show(VTU)
#Simple.SetActiveSource(VTU)
#
PF = Simple.ProgrammableFilter()
PF.Script = Script03
PF.OutputDataSetType = "vtkUnstructuredGrid"
Simple.Show(PF)
#
#Writer = Simple.CreateWriter("xxx.vtu", PF)
#Writer.UpdatePipeline()
#
import os
print "Directory: \'%s\'" % os.getcwd()
print "OK!\n"
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
