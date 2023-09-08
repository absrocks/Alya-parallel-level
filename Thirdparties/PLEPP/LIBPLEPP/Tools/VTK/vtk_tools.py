import vtk 
#import paraview.vtk as vtk  
import os
import sys 
#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#

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



#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __ensight_reader__(filename, reader):
    filename, fileExtension = os.path.splitext( filename )
    print "[ensight_reader] <- \'%s\'" % (filename)

    cdp = vtk.vtkCompositeDataPipeline()
    eg  = vtk.vtkGenericEnSightReader() # vtk.vtkEnSightGoldReader()
    eg.SetDefaultExecutivePrototype(cdp) 

    eg.SetCaseFileName(filename) 
    eg.ReadAllVariablesOn() 
    eg.Update()
    
    reader = eg 


def __vtu_reader__(file_name):
    filename, fileExtension = os.path.splitext( file_name )
    print "[vtu_reader] <- \'%s\'" % (file_name)

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()

    #output = reader.GetOutput() 
    #n_pts  = output.GetNumberOfPoints()
    #print "n_pts:", n_pts

    return reader


def __vtp_reader__(file_name):
    filename, fileExtension = os.path.splitext( file_name )
    print "[vtu_reader] <- \'%s\'" % (file_name)

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(file_name)
    reader.Update()

    return reader


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
      #sys.exit()

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
def __save_mhd__(filename, GetOutput):
    if(not GetOutput.IsA("vtkImageData")): 
      print "[__save_mhd__] ERROR: \'%s\' != \'vtkImageData\' " % GetOutput.GetClassName()
      sys.exit()

    filename, fileExtension = os.path.splitext( filename )
    filename  = filename.upper()
    filename += ".mhd"

    writer = vtk.vtkMetaImageWriter();
    writer.SetInput( GetOutput );
    writer.SetFileName(filename);
    #writer.SetDataModeToAscii();
    writer.Write();

    print "[save_mhd] \'%s\'"% filename


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __save_vti__(filename, GetOutput):
  if(not GetOutput.IsA("vtkImageData")): 
      print "[__save_mhd__] ERROR: \'%s\' != \'vtkImageData\' " % GetOutput.GetClassName()
      sys.exit()

  filename, fileExtension = os.path.splitext( filename )
  filename  = filename.upper()
  filename += ".vti"
  
  writer = vtk.vtkXMLImageDataWriter()
  writer.SetFileName(filename)
  if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInputConnection(GetOutput.GetProducerPort());
  else:
    writer.SetInputData(GetOutput);
  writer.Write();


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __create_vti__(dim = [8,8,8]):
  imageData = vtk.vtkImageData() 
  #whiteImage.SetSpacing(spacing);  
  imageData.SetDimensions(dim);
 
  if vtk.VTK_MAJOR_VERSION <= 5:
    imageData.SetNumberOfScalarComponents(1);
    imageData.SetScalarTypeToDouble();
  else:
    imageData.AllocateScalars(VTK_DOUBLE, 1);

  return imageData


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __surface_filter__(Obj):
  print "[surface_filter] <-\'%s\' ->" % Obj.GetClassName(), 
  GetOutputPort = Obj.GetOutputPort() 
  surfaceFilter = vtk.vtkDataSetSurfaceFilter() 
  surfaceFilter.SetInputConnection( GetOutputPort ) 
  surfaceFilter.Update()

  __save_vtp__(Obj.GetClassName(), surfaceFilter.GetOutput() )

  return surfaceFilter

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __threshold__(Data, v_min, v_max, v_name):
  threshold = vtk.vtkThreshold()
  if vtk.VTK_MAJOR_VERSION <= 5:
    threshold.SetInput(Data);
  else:
    threshold.SetInputData(Data)

  threshold.ThresholdBetween(v_min, v_max) #ThresholdByLower 
  threshold.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, v_name);
  threshold.Update()

  thresholdedPolydata = threshold.GetOutput() #vtkUnstructuredGrid??
  print "[threshold] n_cells: %d \'%s\'" % (thresholdedPolydata.GetNumberOfCells(), v_name)

  return threshold


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __vtp2vtu__(vtp):
    appendFilter = vtk.vtkAppendFilter() 
    if vtk.VTK_MAJOR_VERSION <= 5:
      appendFilter.SetInput( vtp.GetOutput() )
    else:
      appendFilter.SetInputData( vtp.GetOutput() )
    appendFilter.Update()

    unstructuredGrid = vtk.vtkUnstructuredGrid()
    unstructuredGrid.ShallowCopy( appendFilter.GetOutput() )
    unstructuredGrid.Update()

    vtk_n_pts   = unstructuredGrid.GetNumberOfPoints()
    vtk_n_cells = unstructuredGrid.GetNumberOfCells()
    print "[vtp2vtu] n_pts: %d, n_cells: %d " % (vtk_n_pts, vtk_n_cells) 

    return unstructuredGrid


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __vtu2vtp__(vtu):
    PolyData = vtk.vtkPolyData()
    PolyData.ShallowCopy( vtu )
    PolyData.Update()
    return PolyData



#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __vtu2mhd__( vtp,   spacing = [1.0, 1.0, 1.0] ):
  from math import ceil 

  bounds = [0.0]*6;
  vtp.GetBounds(bounds);
  #print bounds

  whiteImage = vtk.vtkImageData() 
  whiteImage.SetSpacing(spacing);
 
  ## compute dimensions
  dim = [0.0, 0.0, 0.0]
  for i in range(3):
    #print (bounds[i * 2 + 1] - bounds[i * 2])/ spacing[i]
    dim[i] = ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i])

  print dim 
  
  whiteImage.SetDimensions(dim);
  whiteImage.SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

 
  origin = [0.0, 0.0, 0.0];
  origin[0] = bounds[0] + spacing[0] / 2;
  origin[1] = bounds[2] + spacing[1] / 2;
  origin[2] = bounds[4] + spacing[2] / 2;
  whiteImage.SetOrigin(origin);
 
  if vtk.VTK_MAJOR_VERSION <= 5:
    whiteImage.SetScalarTypeToUnsignedChar();
    whiteImage.AllocateScalars();
  else:
    whiteImage.AllocateScalars(vtk.VTK_UNSIGNED_CHAR,1);

  ## Fill the image with foreground voxels:
  inval = 255;
  outval = 0;
  count = whiteImage.GetNumberOfPoints()
  for i in range(count):
    whiteImage.GetPointData().GetScalars().SetTuple1(i, inval);


  ## Polygonal data --> image stencil:
  pol2stenc = vtk.vtkPolyDataToImageStencil() 
  if vtk.VTK_MAJOR_VERSION <= 5:
    pol2stenc.SetInput(vtp);
  else:
    pol2stenc.SetInputData(vtp);
  pol2stenc.SetOutputOrigin(origin);
  pol2stenc.SetOutputSpacing(spacing);
  pol2stenc.SetOutputWholeExtent( whiteImage.GetExtent() );
  pol2stenc.Update();

  ## Cut the corresponding white image and set the background:
  imgstenc = vtk.vtkImageStencil() 
  if vtk.VTK_MAJOR_VERSION <= 5:
    imgstenc.SetInput(whiteImage);
    imgstenc.SetStencil(pol2stenc.GetOutput());
  else:
    imgstenc.SetInputData(whiteImage);
    imgstenc.SetStencilConnection(pol2stenc.GetOutputPort());
 
  imgstenc.ReverseStencilOff();
  imgstenc.SetBackgroundValue(outval);
  imgstenc.Update();

  return pol2stenc


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __get_coords__(GetOutput): 
    vtk_data    = GetOutput
    vtk_n_pts   = vtk_data.GetNumberOfPoints()
    vtk_n_cells = vtk_data.GetNumberOfCells()
    if(vtk_n_pts==0): 
      print "[save_data_coords] ERROR: n_pts==0 \n"
      sys.exit() 

    Pts = []
    for idx in range(vtk_n_pts):
      pt = vtk_data.GetPoint(idx)
      Pts.append( pt )
      
    print "[get_coords] n_pts:", vtk_n_pts
    
    return Pts 


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __save_data_coords__(GetOutput, fname=""): 
    vtk_data    = GetOutput
    vtk_n_pts   = vtk_data.GetNumberOfPoints()
    vtk_n_cells = vtk_data.GetNumberOfCells()
    if(vtk_n_pts==0): 
      print "[save_data_coords] ERROR: n_pts==0 \n"
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
      print "[save_data_cell] ERROR: n_cells==0 \n"
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


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __cell_locator__(VTU01): 
  n_Pts = VTU02.GetNumberOfPoints()
  Pts   = VTU02.GetPointData()

  for i in range(n_Pts): 
    cellIdList = vtk.vtkIdList() 
    VTU02.GetPointCells(i,cellIdList)
    print cellIdList.GetNumberOfIds()




#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
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


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
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

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
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


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __cell_locator__xxx__( VTU01, VTU02 ): # no testiado!!
  #http://www.vtk.org/Wiki/VTK/Examples/Cxx/PolyData/CellLocator
  locator = vtk.vtkCellLocator() 
  locator.SetDataSet( VTU01 )
  locator.BuildLocator()

  n_Pts = VTU02.GetNumberOfPoints()
  Pts   = VTU02.GetPointData()
  
  closestPoint      = []*3
  closestPointDist2 = 0.0
  subId             = -1
  cellId            = -1

  n_found = 0
  for i in range(n_Pts): 
    Pt    = VTU02.GetPoint(i)
    locator.FindClosestPoint(Pt, closestPoint, cellId, subId, closestPointDist2);
    print i, PtId

    n_found += 1  

  print "[point_locator]", n_found 

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __set_composite02__(Obj, AppendFilter):
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
    blockIdx = 0
    while( interfaceIterator.IsDoneWithTraversal() == 0 ):
      idx       = interfaceIterator.GetCurrentFlatIndex() 
      dObj      = interfaceIterator.GetCurrentDataObject() #<-- dObj.GetClassName()==vtkPolyData
      blockName = interfaceIterator.GetCurrentMetaData().Get( vtk.vtkCompositeDataSet.NAME() )
      
      
      Output.SetDataSet(interfaceIterator, dObj)
      obj  = Output.GetDataSet(interfaceIterator) 
      grid = vtk.vtkDataSet.SafeDownCast( obj ) 
      
      if(grid.IsA("vtkUnstructuredGrid")):
        Poly = __vtu2vtp__(grid) 
        #AppendFilter.AddInputConnection( Poly.GetProducerPort() )        
        AppendFilter.AddInput( grid )                        #<-- Unstructured??

      else:
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
        CurrentFlatIndex.SetValue(i, idx-1)
      grid.GetCellData().AddArray(CurrentFlatIndex) 

    blockIdx += 1
    AppendFilter.Update() 
    #
    #cleanFilter = vtk.vtkCleanPolyData() 
    #cleanFilter.SetInputConnection( AppendFilter.GetOutputPort() )
    #cleanFilter.Update() 
    #

    print "  |_Total: %d  %d " % (sum(N_Cells), sum(N_Pts) )


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
    blockIdx = 0
    while( interfaceIterator.IsDoneWithTraversal() == 0 ):
      idx       = interfaceIterator.GetCurrentFlatIndex() 
      dObj      = interfaceIterator.GetCurrentDataObject() #<-- dObj.GetClassName()==vtkPolyData
      blockName = interfaceIterator.GetCurrentMetaData().Get( vtk.vtkCompositeDataSet.NAME() )
      
      
      Output.SetDataSet(interfaceIterator, dObj)
      obj  = Output.GetDataSet(interfaceIterator) 
      grid = vtk.vtkDataSet.SafeDownCast( obj ) 
      
      if(grid.IsA("vtkUnstructuredGrid")):
        Poly = __vtu2vtp__(grid) 
        #AppendFilter.AddInputConnection( Poly.GetProducerPort() )        
        AppendFilter.AddInput( grid )                        #<-- Unstructured??

      else:
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
        CurrentFlatIndex.SetValue(i, idx-1)
      grid.GetCellData().AddArray(CurrentFlatIndex) 

    blockIdx += 1
    AppendFilter.Update() 
    #
    #cleanFilter = vtk.vtkCleanPolyData() 
    #cleanFilter.SetInputConnection( AppendFilter.GetOutputPort() )
    #cleanFilter.Update() 
    #

    print "  |_Total: %d  %d " % (sum(N_Cells), sum(N_Pts) )

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __youngs_material_interface__( GetOutput ):
  # Reconstruct material interface
  # http://vtk.org/gitweb?p=VTK.git;a=blob;f=Examples/Graphics/Python/SingleYoungsMaterialInterface.py

  if( not GetOutput.IsA("vtkCompositeDataSet")):
      print "[__youngs_material_interface__] ERROR: \'%s\' != \'vtkCompositeDataSet\' " % GetOutput.GetClassName()
      sys.exit()

  meshMB = vtk.vtkMultiBlockDataSet()
  meshMB.SafeDownCast( GetOutput ) 
  meshMB.Update() 

  interface = vtk.vtkYoungsMaterialInterface()

  if vtk.VTK_MAJOR_VERSION <= 5:
    interface.SetInput(meshMB) 
  else:
    interface.SetInputData(meshMB)

  interface.SetNumberOfMaterials(1)
  interface.SetMaterialVolumeFractionArray(0, "CurrentFlatIndex")

  interface.SetVolumeFractionRange( .001, .999 );
  interface.FillMaterialOn()
  #interface.UseAllBlocksOff()
  #interface.UseAllBlocksOn()
  interface.Update()

  Merge = vtk.vtkMultiBlockMergeFilter()
  if vtk.VTK_MAJOR_VERSION <= 5:
      Merge.AddInput( interface.GetOutput() ) 
  else:
      Merge.AddInputData( interface.GetOutput() ) 
  Merge.Update() 

  #__save_vtp__("XXX", Merge.GetOutput() )
  __save_vtm__("XXX", Merge.GetOutput() )
  #__save_vtm__("XXX", interface.GetOutput() )

  return interface


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __vtkBlockIdScalars__( meshMB ):
  #BlockId = vtk.vtkBlockIdScalars() 
  #BlockId.SetInputData( meshMB )  
  #BlockId.SetInput( meshMB )
  #BlockId.Update() 
  
  #__get_cell_data__( BlockId ) 


  __save_vtm__("XXX", meshMB)
  
  return  BlockId


#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def __vtkStreamTracer__( GetOutputPort ):
  sphereSource = vtk.vtkSphereSource() 
  sphereSource.SetCenter(0.0, 0.0, 0.0)
  sphereSource.SetRadius(0.1)
  sphereSource.SetThetaResolution(10)
  sphereSource.SetPhiResolution (10)
  sphereSource.Update()

  LineSourceWidget0 = vtk.vtkLineSource()
  LineSourceWidget0.SetPoint1(0.0, 0.0,  1e-3)
  LineSourceWidget0.SetPoint2(0.0, 0.0, -1e-3)
  LineSourceWidget0.SetResolution(100)

  Stream = vtk.vtkStreamTracer()
  Stream.SetInput( LineSourceWidget0.GetOutput() )  
  Stream.SetInputConnection( GetOutputPort )
  Stream.SetStartPosition(0.0, 0.0, 0.0) 
  #Stream.SetMaximumPropagation 500
  Stream.SetIntegrationDirectionToBoth() 
  Stream.SetIntegratorTypeToRungeKutta4 ()
  Stream.SetComputeVorticity (False)
  
  
  return Stream

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
#def __vtkThreshold__( GetOutputPort ):



#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
def linear(x, x0, y0, x1, y1):
    return y0 + (y1-y0)*(x-x0)/(x1-x0)

#--------------------------------------------------------------------------||--#
#--------------------------------------------------------------------------||--#
