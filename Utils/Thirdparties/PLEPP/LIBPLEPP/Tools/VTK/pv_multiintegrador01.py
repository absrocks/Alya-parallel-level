import os
import re
import glob
import sys

def Read_dat(fname):
    fname = glob.glob(fname)
    if(len(fname)<1):
      print "Error: \'%s\'! \n" % fname
      exit(1)
  
    fname = fname[-1]
    data = open(fname, "r")
    lines = data.readlines()
    data.close()

    nline = len(lines)
    XYZ = []
    for i in range(nline):
      line = lines[i]
      line = line.strip()
      line = line.split()
      XYZ.append([eval(val) for val in line[1:]])

    print "  |_No elements:", len(XYZ)
    return XYZ


TOOLS = """
import sys 
#-------------------------------------------------------------------------||---#
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
    #print>> fcoords, "COORDINATES"
    for idx in range(vtk_n_pts):
      print>> fcoords, idx+1, 
      pt = vtk_data.GetPoint(idx)
      n_pt = len(pt)
      for j in range(n_pt): print>> fcoords, pt[j],
      print>> fcoords, ""

    #print>> fcoords, "END_COORDINATES"
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


def __save_pts_data__(GetOutput, fname):
    Data    = GetOutput.GetPointData() 
    n_props = Data.GetNumberOfArrays()

    Fout = open(fname, "a")
    for i in range(n_props):
      property_name = Data.GetArrayName(i)
      print " |_%d) \'%s\'" % (i+1, property_name)
      
      Prop   = Data.GetArray(property_name)
      n_cols = Prop.GetNumberOfComponents()
      n_rows = Prop.GetNumberOfTuples() 
      print "\'%s\':" % Prop.GetName(),   

      prop_range = Prop.GetRange()   
      print "[%f,%f]" %( prop_range[0], prop_range[1] )

      k=0
      for j in range(n_rows): 
          for i in range(n_cols): 
            val = Prop.GetValue(k)
            print>> Fout, val, 
            k+=1
    print>> Fout, ""
    
"""

END01 = """
#--------------------------------------------------------------------------||--#
# Get the two inputs
A = self.GetInputDataObject(0, 0)
print "A:", A.GetClassName()

if( A.IsA("vtkUnstructuredGrid") ): 
  func01(A)
elif( A.IsA("vtkPolyData") ): 
  VTU = __vtp2vtu__(A)
  func01(VTU)
else:
  print "[Alya] ERROR: \'%s\' != \'vtkUnstructuredGrid\' or \'vtkPolyData\' " % A.GetClassName()
  sys.exit()
#
#output = self.GetOutput()
#output.ShallowCopy(VTU_SURFACE)
#output.GetPointData().AddArray(labels)
#
#-------------------------------------------------------------------------||---#
"""

#-------------------------------------------------------------------------||---#
Script00 = TOOLS + """
#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def func01(A):
  #__save_data_coords__(A, "XXX") 
  #__save_data_cell__(A, "A") 
  __get_pts_data__(A) 
#--------------------------------------------------------------------------||--#
""" + END01
#-------------------------------------------------------------------------||---#


#-------------------------------------------------------------------------||---#
Script01 = TOOLS + """
#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def func01(A):
  __save_data_coords__(A, "XXX") 
  #__save_data_cell__(A, "A") 
  #__get_pts_data__(A) 
#--------------------------------------------------------------------------||--#
""" + END01
#-------------------------------------------------------------------------||---#


#-------------------------------------------------------------------------||---#
Script02 = TOOLS + """
#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def func01(A):
      #__save_data_coords__(A, "XXX") 
      #__save_data_cell__(A, "A") 

      import datetime
      now = datetime.datetime.now()
      property_name  = "YYY"
      property_name += "_%d" % now.year
      property_name += "_%d" % now.month
      property_name += "_%d" % now.hour
      #property_name += "_%d" % now.minute
      #property_name += "_%d" % now.second
      property_name += ".dat"
      __save_pts_data__(A, property_name) 
#--------------------------------------------------------------------------||--#
""" + END01 
#-------------------------------------------------------------------------||---#



#-------------------------------------------------------------------------||---#
#-------------------------------------------------------------------------||---#
import paraview.simple as Simple 


#-------------------------------------------------------------------------||---#
PATH = ['/home/jmake/z2014/Alya/Runner/Couplings/TUBES03/GRAETZ01/xxx/xxx_0_0.vtu']
#
#InputData = Simple.XMLUnstructuredGridReader( 
#                        guiName="xxxx_0_0.vtu", 
#                        PointArrayStatus=['TEMPE', 'TU'], 
#                        FileName=PATH)
#
PATH = '/home/jmake/z2014/Alya/Runner/Couplings/TUBES03/GRAETZ01/internal01.ensi.case'
InputData = Simple.EnSightReader( 
                       guiName="internal01.ensi.case", 
                       PointArrays=['TEMPE', 'VELOC', "HEATF"], 
                       CaseFileName=PATH)
#-------------------------------------------------------------------------||---#

#-------------------------------------------------------------------------||---#
PlotOverLine02 = Simple.PlotOverLine( guiName="PlotOverLine02", Source="High Resolution Line Source" )
PlotOverLine02.Source.Resolution = 100
PlotOverLine02.Source.Point2 = [0.99, 0.0, 20.0]
PlotOverLine02.Source.Point1 = [0.99, 0.0, 0.0]
#
PF = Simple.ProgrammableFilter()
PF.Script = Script00
Simple.Show(PF)
#-------------------------------------------------------------------------||---#

#-------------------------------------------------------------------------||---#
PlotOverLine01 = Simple.PlotOverLine( guiName="PlotOverLine01", Source="High Resolution Line Source" )
PlotOverLine01.Source.Resolution = 100
PlotOverLine01.Source.Point2 = [0.00, 0.0, 20.0]
PlotOverLine01.Source.Point1 = [0.00, 0.0, 0.0]
#
PF = Simple.ProgrammableFilter()
PF.Script = Script01
Simple.Show(PF)
#-------------------------------------------------------------------------||---#

#-------------------------------------------------------------------------||---#
Coords = Read_dat("XXX_COORDINATES.alya")
for coord in Coords: 
  print coord, 
  Simple.SetActiveSource(InputData)
  Slice01 = Simple.Slice( guiName="Slice2", SliceOffsetValues=[0.0], SliceType="Plane" )
  Slice01.SliceType.Origin = coord
  Slice01.SliceType.Normal = [0.0, 0.0, 1.0]
  Simple.SetActiveSource(Slice01)

  Calculator01 = Simple.Calculator(guiName="Calculator1", 
                            Function='VELOC_Z*TEMPE', 
                            ResultArrayName='UT')
  Simple.SetActiveSource(Calculator01)

  IntegrateVariables01 = Simple.IntegrateVariables( guiName="IntegrateVariables2" )
  Simple.SetActiveSource(IntegrateVariables01)

  PF = Simple.ProgrammableFilter()
  PF.Script = Script02 
  Simple.Show(PF)
  #print "ok!"
#-------------------------------------------------------------------------||---#


#-------------------------------------------------------------------------||---#
#-------------------------------------------------------------------------||---#


#-------------------------------------------------------------------------||---#
import os
print "Directory: \'%s\'" % os.getcwd()

print "OK!\n"
#-------------------------------------------------------------------------||---#
