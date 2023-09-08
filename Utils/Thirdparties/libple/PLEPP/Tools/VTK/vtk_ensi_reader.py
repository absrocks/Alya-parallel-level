import vtk 
import os
import sys 
#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __save_vtu__(filename, GetOutput):
    if(not GetOutput.IsA("vtkUnstructuredGrid")): 
      print "[__save_vtu__] ERROR: \'%s\' != \'vtkUnstructuredGrid\' " % GetOutput.GetClassName()
      sys.exit()

    filename, fileExtension = os.path.splitext( filename )
    filename += ".vtu"

    writer = vtk.vtkXMLUnstructuredGridWriter() 
    writer.SetFileName(filename)
    writer.SetInput( GetOutput )
    writer.Write()
    
    print "[save_vtu] \'%s\'"% filename

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __save_vtp__(filename, GetOutput):
    if(not GetOutput.IsA("vtkPolyData")): 
      print "[__save_vtp__] ERROR: \'%s\' != \'vtkPolyData\' " % GetOutput.GetClassName()
      sys.exit()

    filename, fileExtension = os.path.splitext( filename )
    filename += ".vtp"

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(filename)
    writer.SetInput( GetOutput )
    writer.Write()

    print "[save_vtp] \'%s\'"% filename


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __save_vtm__(filename, GetOutput):
    if(not GetOutput.IsA("vtkMultiBlockDataSet")): 
      print "[__save_vtp__] ERROR: \'%s\' != \'vtkMultiBlockDataSet\' " % GetOutput.GetClassName()
      sys.exit()

    filename, fileExtension = os.path.splitext( filename )
    filename  = filename.upper()
    filename += ".vtm"

    writer = vtk.vtkXMLMultiBlockDataWriter();
    writer.SetInput( GetOutput );
    writer.SetFileName(filename);
    #writer.SetDataModeToAscii();
    writer.Write();

    print "[save_vtm] \'%s\'"% filename


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __set_composite__(Obj, Merge, AppendFilter):
    if(not Obj.IsA("vtkMultiBlockDataSet")): 
      print "[__set_composite__] ERROR: \'%s\' != \'vtkMultiBlockDataSet\' " % Obj.GetClassName()
      sys.exit()
    
    meshMB = Obj

    if vtk.VTK_MAJOR_VERSION <= 5:
      Merge.AddInput( meshMB ) 
    else:
      Merge.AddInputData( meshMB ) 
    Merge.Update() 

    Input  = Obj
    Output = Input.NewInstance()    #<-- class vtkCompositeDataIterator
    Output.CopyStructure( Input ) 

    interfaceIterator = vtk.vtkCompositeDataIterator()
    interfaceIterator.SetDataSet( meshMB )
    interfaceIterator.VisitOnlyLeavesOn()
    interfaceIterator.SkipEmptyNodesOn()
    interfaceIterator.InitTraversal()
    interfaceIterator.GoToFirstItem()

    print "  [Multiblock]"
    N_Cells = []
    N_Pts   = [] 
    while( interfaceIterator.IsDoneWithTraversal() == 0 ):
      idx       = interfaceIterator.GetCurrentFlatIndex() 
      dObj      = interfaceIterator.GetCurrentDataObject() #<-- dObj.GetClassName()==vtkPolyData
      blockName = interfaceIterator.GetCurrentMetaData().Get( vtk.vtkCompositeDataSet.NAME() )

      Output.SetDataSet(interfaceIterator, dObj)      
      obj  = Output.GetDataSet(interfaceIterator) 
      grid = vtk.vtkDataSet.SafeDownCast( obj ) 
      AppendFilter.AddInput( grid )                        #<-- Unstructured??
      interfaceIterator.GoToNextItem()

      print "  |_%d) "    % idx, 
      print "n_pts: %d,"  % grid.GetNumberOfPoints(), 
      print "n_cells: %d" % grid.GetNumberOfCells(),
      #print grid.GetClassName(), 
      if(blockName!=None): 
        print "\'%s\' " % blockName.strip()
      else:
        print ""

      n_cells = grid.GetNumberOfCells()      
      n_pts   = grid.GetNumberOfPoints()
      pts     = grid.GetPointData()

      N_Cells.append(n_cells)
      N_Pts.append(n_pts)

      CurrentFlatIndex = vtk.vtkIntArray()
      CurrentFlatIndex.SetName("CurrentFlatIndex")
      CurrentFlatIndex.SetNumberOfComponents( 1 )
      CurrentFlatIndex.SetNumberOfTuples( n_cells )

      for i in range( CurrentFlatIndex.GetNumberOfTuples() ): 
        CurrentFlatIndex.SetValue(i, idx)
      grid.GetCellData().AddArray(CurrentFlatIndex) 

    AppendFilter.Update() 
    #
    #cleanFilter = vtk.vtkCleanPolyData() 
    #cleanFilter.SetInputConnection( AppendFilter.GetOutputPort() )
    #cleanFilter.Update() 
    #

    print "  |_Total: %d  %d " % (sum(N_Cells), sum(N_Pts) )

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
class Vtk_ensight_reader():
  def __init__(self):
    cdp = vtk.vtkCompositeDataPipeline()
    #self.eg = vtk.vtkEnSightGoldReader()
    self.eg = vtk.vtkGenericEnSightReader()
    self.eg.SetDefaultExecutivePrototype(cdp) 


  def init(self, filename=""):
    self.filename, self.fileExtension = os.path.splitext( filename )
    print "[Vtk_ensight_reader] <- \'%s\'" % (filename)

    self.eg.SetCaseFileName(filename) 
    self.eg.ReadAllVariablesOn() 
    self.eg.Update()

    self.nVariables = self.eg.GetNumberOfVariables() 
    for i in range(self.nVariables):
      print "%d)" % (i+1), 
      print "\'%s\'" % self.eg.GetPointArrayName(i) 


  def set_multiblock_merge(self):
    meshMB = self.eg.GetOutput()

    Merge = vtk.vtkMultiBlockMergeFilter()
    if vtk.VTK_MAJOR_VERSION <= 5:
      Merge.AddInput( meshMB ) 
    else:
      Merge.AddInputData( meshMB ) 
    Merge.Update() 

    self.Merge = Merge.GetOutput().GetBlock(0)

    self.PtsIds = vtk.vtkIntArray()
    self.PtsIds.SetName("PtsIds")
    self.PtsIds.SetNumberOfComponents( meshMB.GetNumberOfBlocks() )
    self.PtsIds.SetNumberOfTuples( self.Merge.GetNumberOfPoints() )

    k = 0
    for i in range( self.Merge.GetNumberOfPoints() ): 
      for j in range( meshMB.GetNumberOfBlocks() ): 
        self.PtsIds.SetValue(k, -1)
        k += 1

    self.CellIds = vtk.vtkIntArray()
    self.CellIds.SetName("CellIds")
    self.CellIds.SetNumberOfComponents( meshMB.GetNumberOfBlocks() )
    self.CellIds.SetNumberOfTuples( self.Merge.GetNumberOfCells() )

    k = 0
    for i in range( self.Merge.GetNumberOfCells() ): 
      for j in range( meshMB.GetNumberOfBlocks() ): 
        self.CellIds.SetValue(k, -(j+1))
        k += 1

    self.pointLocator = vtk.vtkPointLocator() 
    self.pointLocator.SetDataSet( self.Merge )
    self.pointLocator.BuildLocator()

    print "[set_multiblock_merge] "    


  def set_composite(self): 
    self.Merge        = vtk.vtkMultiBlockMergeFilter()
    self.AppendFilter = vtk.vtkAppendFilter()
    self.Obj = self.eg.GetOutput()
    __set_composite__(self.Obj, self.Merge, self.AppendFilter)


  def __set_composite__(self, idx=-1):
    Input  = self.eg.GetOutput()
    Output = Input.NewInstance()    #class vtkCompositeDataIterator
    Output.CopyStructure( Input ) 

    meshMB = self.eg.GetOutput()
    #for i in range(  meshMB.GetNumberOfBlocks() ): 
    #  print meshMB.HasMetaData( i )
    #  print meshMB.GetMetaData( i )

    AppendFilter = vtk.vtkAppendFilter()
    IDs = [[]]*self.Merge.GetNumberOfPoints()

    interfaceIterator = vtk.vtkCompositeDataIterator()
    interfaceIterator.SetDataSet( meshMB )
    interfaceIterator.VisitOnlyLeavesOn()
    interfaceIterator.SkipEmptyNodesOn()
    interfaceIterator.InitTraversal()
    interfaceIterator.GoToFirstItem()
    while( interfaceIterator.IsDoneWithTraversal() == 0 ):
      idx       = interfaceIterator.GetCurrentFlatIndex() 
      dObj      = interfaceIterator.GetCurrentDataObject();
      blockName = interfaceIterator.GetCurrentMetaData().Get( vtk.vtkCompositeDataSet.NAME() )

      Output.SetDataSet(interfaceIterator, dObj)      
      obj  = Output.GetDataSet(interfaceIterator)
      grid = vtk.vtkDataSet.SafeDownCast( obj )

      VTU = vtk.vtkUnstructuredGrid() 
      VTU.ShallowCopy( obj ) 
      AppendFilter.AddInputConnection( VTU.GetProducerPort() )

      n_pts = VTU.GetNumberOfPoints()
      pts   = VTU.GetPointData()

      n_found = 0
      for i in range(n_pts): 
        pt    = VTU.GetPoint(i)
        ptId  = self.pointLocator.FindClosestPoint(pt)
        IDs[ptId]  = IDs[ptId] + [idx]
        #if(len(IDs[ptId])>0): print ptId, IDs[ptId]

        #l = ptId * (meshMB.GetNumberOfBlocks())  + idx
        #self.PtsIds.SetValue(l, idx)

        cellIdList = vtk.vtkIdList()
        self.Merge.GetPointCells(ptId, cellIdList);
        n_cellIdList = cellIdList.GetNumberOfIds()  
        for j in range(n_cellIdList):
          #Ids = [cellIdList.GetId(j)  for j in range(n_cellIdList) ]
          #print ptId, n_cellIdList, Ids
          cellId = cellIdList.GetId(j)
          
          k = cellId*(meshMB.GetNumberOfBlocks())  + idx
          self.CellIds.SetValue(k, idx+1)

        n_found += 1

      print idx, dObj.GetClassName(), "\'%s\'" % blockName.strip(), 
      print grid.GetNumberOfPoints(), grid.GetNumberOfCells(), n_found
      interfaceIterator.GoToNextItem()

    AppendFilter.Update() 
    #cleanFilter = vtk.vtkCleanPolyData() 
    #cleanFilter.SetInputConnection( AppendFilter.GetOutputPort() )
    #cleanFilter.Update() 

    self.Output = AppendFilter.GetOutput() 


    self.Array01 = vtk.vtkIntArray()
    self.Array01.SetName("Array01")
    self.Array01.SetNumberOfComponents( 1 )
    self.Array01.SetNumberOfTuples( self.Merge.GetNumberOfPoints() )
    
    for i in range( self.Merge.GetNumberOfPoints() ): 
      n_IDs = len(IDs[i]) 
      #print IDs[i]
      self.Array01.SetValue(i, n_IDs )

      #pt    = VTU.GetPoint(i)
      #ptId  = self.pointLocator.FindClosestPoint(pt);        

      #l = ptId * (meshMB.GetNumberOfBlocks())  + idx
      #self.PtsIds.SetValue(l, idx)
    
    
    self.Merge.GetPointData().AddArray(self.Array01) 

    self.Merge.GetPointData().AddArray(self.PtsIds) 
    self.Merge.GetCellData().AddArray(self.CellIds) 

    filename  = self.filename
    filename += "_merge"
    filename += ".vtu"

    writer = vtk.vtkXMLUnstructuredGridWriter() 
    writer.SetFileName(filename);
    writer.SetInput( self.Merge )
    writer.Write()

    print "[set_composite] \'%s\'"% filename    


  def set_block(self, idx=0): 
    Vtu = self.eg.GetOutput().GetBlock(idx)
    self.vtu.SetInput( Vtu ) 
    self.vtu.Update() 


  def get_vtu_output_port(self):
    return self.vtu.GetOutputPort()


  def get_data(self, property_name=""):
    Vtu = self.eg.GetOutput().GetBlock(self.BLOCKid)
    Datai = Vtu.GetPointData() 
    
    Propi  = Datai.GetArray(property_name)
    nColsi = Propi.GetNumberOfComponents()
    nRowsi = Propi.GetNumberOfTuples() 
    print "\'%s\':" % Propi.GetName(),   
    print "(%d,%d)" % (nColsi, nRowsi)  


  def save_vtu(self):
    print "[Vtk_ensight_reader]", 
    __save_vtu__(self.filename, self.AppendFilter.GetOutput())


  def __save_vtu__(self, BLOCKid=0, filename=""):
    filename  = self.filename
    filename += "_"
    filename += str( BLOCKid ).zfill(3)
    filename += ".vtu"

    writer = vtk.vtkXMLUnstructuredGridWriter() 
    writer.SetFileName(filename);
    writer.SetInput( self.Output )

    #writer.SetInput( self.vtu.GetOutput() )
    writer.Write()
    
    print "[save_vtu] \'%s\'"% filename


  def save_multiblock(self):
    print "[Vtk_ensight_reader]", 
    __save_vtm__(self.filename, self.eg.GetOutput())


  def __save_multiblock__(self):
    filename  = self.filename
    filename  = filename.upper()
    filename += ".vtm"

    MB = self.eg.GetOutput()
    writer = vtk.vtkXMLMultiBlockDataWriter();
    writer.SetInput( MB );
    writer.SetFileName(filename);
    writer.SetDataModeToAscii();
    writer.Write();

    print "[save_multiblock] \'%s\'"% filename

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
class Vtk_multiblock_reader():
  def __init__(self):
    self.reader = vtk.vtkXMLMultiBlockDataReader() 


  def init(self, filename=""):
    self.filename, self.fileExtension = os.path.splitext( filename )
    print "[Vtk_multiblock_reader] <- \'%s\'" % (filename)

    self.reader.SetFileName(filename)
    self.reader.Update()
    self.Obj = self.reader.GetOutput()

  def set_composite(self):  
    self.Merge        = vtk.vtkMultiBlockMergeFilter()
    self.AppendFilter = vtk.vtkAppendFilter()
    __set_composite__(self.Obj, self.Merge, self.AppendFilter)

    #unstructuredGrid = vtk.vtkUnstructuredGrid() 
    #unstructuredGrid.ShallowCopy( self.Merge.GetOutput().GetBlock(0) )


  def save_vtp(self):
    print "[Vtk_multiblock_reader]", 
 
    filename  = self.filename
    filename += "_merge"
    __save_vtp__(filename, self.Merge.GetOutput().GetBlock(0) ) 


  def save_vtu(self):
    print "[Vtk_multiblock_reader]", 
    __save_vtu__(self.filename, self.AppendFilter.GetOutput() ) 


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
Fin  = "SURFACE.vtm"
Fout = "XXX"

MB = Vtk_multiblock_reader()
MB.init(Fin) 
MB.set_composite()
MB.save_vtp() 
MB.save_vtu() 


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
Fin  = "star.case"
Fout = "star"

ENSI = Vtk_ensight_reader()
ENSI.init(Fin) 

#ENSI.set_multiblock_merge()
ENSI.set_composite() 

ENSI.save_multiblock() 
ENSI.save_vtu() 

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
