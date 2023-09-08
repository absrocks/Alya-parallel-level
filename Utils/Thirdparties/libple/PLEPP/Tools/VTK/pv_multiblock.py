import paraview.simple as Simple 

Script03 = """
def __cell_locator__(VTU01, VTU02, Cell_idx): 
  #locator = vtk.vtkOctreePointLocator() 
  #locator.SetDataSet( VTU01 )
  #locator.BuildLocator()

  #locator = vtk.vtkPointLocator()
  #locator.InitPointInsertion(vtk.vtkPoints(), VTU01.GetBounds())
  #for i in range(VTU01.GetNumberOfPoints()):
  #  locator.InsertNextPoint(VTU01.GetPoints().GetPoint(i))

  locator = vtk.vtkPointLocator() 
  #locator = vtk.vtkOctreePointLocator() 
  locator.SetDataSet( VTU01 )
  locator.BuildLocator()

  n_Cells = VTU02.GetNumberOfCells()
  for idx in range(n_Cells):  
    cell_ids = vtk.vtkIdList()
    VTU02.GetCellPoints(idx, cell_ids) 
    n_cell_ids = cell_ids.GetNumberOfIds()  

    for j in range(n_cell_ids): 
      i     = cell_ids.GetId(j)
      Pt    = VTU02.GetPoint(i)      
      PtId  = locator.FindClosestPoint(Pt)

      #if locator.IsInsertedPoint(Pt) != -1:
      #  Cell_idx[PtId] += 1
      #else:
      #  Cell_idx[PtId] = -1

  #return Cell_idx


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

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
#import sys
#sys.path.insert(0, "/home/jmake/z2014/Alya/Repository/Tools/VTK_Develoment")
#sys.path.insert("")
#
#import vtk_multiblock 
#import vtk_tools as tools 
#

# Get the two inputs
A = self.GetInputDataObject(0, 0)
B = self.GetInputDataObject(0, 1)

print "A:", A.GetClassName()
print "B:", B.GetClassName()

if( not A.IsA("vtkMultiBlockDataSet") ): 
  print "[Alya] ERROR: \'%s\' != \'vtkUnstructuredGrid\' " % A.GetClassName()
  #sys.exit()

if( not B.IsA("vtkUnstructuredGrid") ): 
  print "[Alya] ERROR: \'%s\' != \'vtkUnstructuredGrid\' " % B.GetClassName()
  #sys.exit()


# Create an array to store the labels
#labels = vtk.vtkUnsignedCharArray()
labels = vtk.vtkDoubleArray() 
labels.SetName("IDs")
#labels.SetNumberOfTuples(B.GetNumberOfCells())
labels.SetNumberOfTuples(B.GetNumberOfPoints())
#labels.SetNumberOfComponents( A.GetNumberOfBlocks() )

Cell_idx = [ [] ]*B.GetNumberOfPoints()

if A.IsA("vtkMultiBlockDataSet"):
    iter = A.NewIterator()
    iter.UnRegister(None)
    iter.InitTraversal()
    while not iter.IsDoneWithTraversal():
        curInput  = iter.GetCurrentDataObject()
        idx       = iter.GetCurrentFlatIndex() 
        blockName = iter.GetCurrentMetaData().Get( vtk.vtkCompositeDataSet.NAME() )
        print "%d) \'%s\' " % (idx, blockName.strip()), 
        
        curOutput = curInput.NewInstance()
        curOutput.UnRegister(None)

        #----------------------------------------------------------------------#
        #CurrentFlatIndex = vtk.vtkIntArray()
        #CurrentFlatIndex.SetName("CurrentFlatIndex")
        #CurrentFlatIndex.SetNumberOfComponents( 1 )
        #CurrentFlatIndex.SetNumberOfTuples( n_cells )
        #
        #for i in range( CurrentFlatIndex.GetNumberOfTuples() ): 
        #  CurrentFlatIndex.SetValue(i, idx-1)
        #curInput.GetCellData().AddArray(CurrentFlatIndex) 
        #----------------------------------------------------------------------#

        #----------------------------------------------------------------------#
        n_cells = curInput.GetNumberOfCells()   
        n_pts   = curInput.GetNumberOfPoints()
        pts     = curInput.GetPointData()
        print "n_pts, n_cells", n_pts, n_cells
        
        GlobalIdx = __point_locator__(B, curInput)
        for Gidx in GlobalIdx: labels.SetValue(Gidx, idx)
        #for Gidx in GlobalIdx: Cell_idx[Gidx].append( idx ) 

        #----------------------------------------------------------------------#
        iter.GoToNextItem()

#for i in range(B.GetNumberOfPoints()): labels.SetValue(i, len(Cell_idx[i]))

output = self.GetOutput()
output.ShallowCopy(B)
output.GetPointData().AddArray(labels)
#output.GetCellData().AddArray(labels)

PROP = __get_pts_data__(output, "IDs")
__save_raw__(PROP, "IDs.alya")
"""

PF = Simple.ProgrammableFilter()
PF.Script = Script03
Simple.Show(PF)

import os
print "Directory: \'%s\'" % os.getcwd()

print "OK!\n"
