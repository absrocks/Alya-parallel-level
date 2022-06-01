import vtk 
import os
import sys 
#==============================================================================#
#==============================================================================#

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
class Vtk_creator:
  def __init__(self):
    self.obj       = None 
    self.num_pts   = 0
    self.vtk_pts   = vtk.vtkPoints() 
    self.vtk_cells = vtk.vtkCellArray()
    self.vtk_tetra = vtk.vtkTetra()

    self.vtk_props = [] 
    self.vtk_cells_type = []
    self.vtk_cells_id = []
    self.vtk_cells_props = [] 

    self.n_cells = 0
    self.n_pts = 0

    print "\n+Vtk_writer"
    

  def set_point(self, coord):
    self.vtk_pts.InsertNextPoint( coord[0], coord[1], coord[2] )
    self.num_pts += 1


  def set_points_prop_scalar(self, R, prop_name=""):    
    array = vtk.vtkDoubleArray()
    array.SetName(prop_name)
    array.SetNumberOfComponents(1)
    array.SetNumberOfTuples( self.num_pts )

    for i in range(self.num_pts): array.SetValue(i, R[i])
    #self.obj.GetPointData().AddArray(array)
    #self.obj.Update()

    self.vtk_props.append( array ) 
    print "|_Vtk_writer.set_vertices_prop: \'%s\' " % (prop_name)


  def set_points_prop_vector(self, R, prop_name=""):
    array = vtk.vtkDoubleArray()
    array.SetName(prop_name)
    array.SetNumberOfComponents(3)
    array.SetNumberOfTuples( self.num_pts )

    j = 0
    for i in range(self.num_pts): 
      for k in range(3): array.SetValue(j+k, R[i,k])
      j += 3
    #self.pdi.GetPointData().AddArray(array)

    self.vtk_props.append( array ) 
    print "|_Vtk_writer.set_vertices_prop_vec: \'%s\' " % (prop_name)


  def set_tetra(self, nodes):
    self.__set_cell__(nodes, vtk.vtkTetra())


  def set_triangle(self, nodes):
    self.__set_cell__(nodes, vtk.vtkTriangle())


  def set_cell(self, nodes, vtk_cell_type):
    self.__set_cell__(nodes, vtk_cell_type)


  def __set_cell__(self, nodes, vtk_cell):
    i = 0
    vtk_list = vtk.vtkIdList() 
    for node in nodes: 
      #vtk_cell.GetPointIds().SetId(i, node-1)
      if (node<0): 
        print "Negative node!!\n"
        exit(1)
        
      vtk_list.InsertNextId(node) 
      i += 1
    #self.vtk_cells.InsertNextCell(vtk_cell)
    self.vtk_cells_type.append( vtk_cell.GetCellType() ) 
    self.vtk_cells_id.append(vtk_list) 

    self.n_cells += 1


  def set_cell_prop_scalar(self, R, prop_name=""):
    array = vtk.vtkDoubleArray()
    array.SetName(prop_name)
    array.SetNumberOfComponents(1)
    array.SetNumberOfTuples( self.n_cells )

    for i in range(self.n_cells): array.SetValue(i, R[i])
    self.vtk_cells_props.append( array ) 


  def set_vtk_unstructured(self): 
    vtk_data = vtk.vtkUnstructuredGrid()
    vtk_data.SetPoints( self.vtk_pts ) 
    #vtk_data.GetClientSideObject().GetOutputDataObject(0)

    for cell_type, cell_ids in zip(self.vtk_cells_type, self.vtk_cells_id): 
      vtk_data.InsertNextCell(cell_type, cell_ids) 

    for array in self.vtk_props: 
      vtk_data.GetPointData().AddArray(array)
     #vtk_data.GetPointData().SetScalars(array)

    for array in self.vtk_cells_props: 
      vtk_data.GetCellData().AddArray(array)
      #vtk_data.Update() 

    
    #vtk_data.GetPointData().AddArray(array)
    #vtk_data.SetCells(self.vtk_cells)        
    #if vtk.VTK_MAJOR_VERSION <= 5:
    #vtk_data.Update()
    
    self.obj = vtk_data 
    print "|_Vtk_writer.to_vtkunstructured_grid: ", 
    print "pts/cells: %d %d" %(self.num_pts, self.n_cells) 


  def save_vtk_unstructured(self, file_name): 
    self.file_name = os.path.splitext(file_name)[0]

    #file_name_out  = self.file_name+".vtk"
    #writer = vtk.vtkUnstructuredGridWriter()

    file_name_out  = self.file_name+".vtu"
    writer = vtk.vtkXMLUnstructuredGridWriter() 
    writer.SetFileName(file_name_out)

    #if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(self.obj) 
    #else:
    #writer.SetInputData(self.obj)

    writer.Write() 
    print "|_Vtk_writer.writer file: \'%s\' " % (file_name_out)


  def psave(self, file_name, mpi_size, mpi_range):
    self.file_name = os.path.splitext(file_name)[0] 

    file_name_out  = self.file_name+"_%d" % mpi_range + ".vtu"
    writer = vtk.vtkXMLUnstructuredGridWriter() 
    writer.SetFileName(file_name_out)
    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput( self.obj ) #self.reader.GetOutput())
    else:
      writer.SetInputData( self.obj )#self.reader.GetOutput())
    writer.Write()

    if(mpi_range==0): 
      pwriter = vtk.vtkXMLPUnstructuredGridWriter() 
      pwriter.SetFileName(self.file_name+".pvtu")
      pwriter.SetNumberOfPieces(mpi_size)

      if vtk.VTK_MAJOR_VERSION <= 5:
        pwriter.SetInput( self.obj ) #self.reader.GetOutput())
      else:
        pwriter.SetInputData( self.obj ) #self.reader.GetOutput())

      pwriter.Write()
      print "|_Vtk_writer.writer file: \'%s\' " % (file_name_out)

#--------------------------------------------------------------------------||--#



#--------------------------------------------------------------------------||--#
class Vtk_structured_creator:
  def __init__(self):
    self.vtk_pts   = vtk.vtkPoints() 
    self.vtk_cells = vtk.vtkCellArray()


  def set_point(self, coord):
    self.vtk_pts.InsertNextPoint( coord[0], coord[1], coord[2] )


  def set_points_prop_scalar(self, R, prop_name=""): 
    array = vtk.vtkDoubleArray()
    array.SetName(prop_name)
    array.SetNumberOfComponents(1)
    array.SetNumberOfTuples( self.num_pts )

    for i in range(self.num_pts): array.SetValue(i, R[i])

    self.vtk_props.append( array ) 
    print "|_Vtk_writer.set_vertices_prop: \'%s\' " % (prop_name)


  def set_vtk_unstructured(self): 
    vts = vtk.vtkStructuredGrid()
    vts.SetDimensions(Nz,Ny,Nx)
    vts.SetPoints(pts)

    for array in self.vtk_props: 
      vts.GetPointData().AddArray(array)

    self.obj = vtk_data 
    print "|_Vtk_writer.to_vtkstructured"


  def save_vtk_unstructured(self, file_name): 
    self.file_name = os.path.splitext(file_name)[0]
    file_name_out  = self.file_name+".vtk"

    writer = vtk.vtkStructuredGridWriter()
    writer.SetFileName(file_name_out)

    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput( self.obj ) 
    else:
      writer.SetInputData( self.obj )

    writer.Write() 
    print "|_Vtk_writer.writer file: \'%s\' " % (file_name_out)

#--------------------------------------------------------------------------||--#

#--------------------------------------------------------------------------||--#
class Vtk_multiblock_reader:
  def __init__(self):
    self.vtk_pts   = None
    self.vtk_cells = None 
    self.prop_range = {}
    #print "\n+[Vtk_multiblock_reader]:"

  def load_vtk_structured(self, file_name): 
    self.fname = file_name
    #reader = vtk.vtkUnstructuredGridReader()
    reader = vtk.vtkXMLMultiBlockDataReader()
    reader.SetFileName(file_name)
    reader.Update()

    Merge = vtk.vtkMultiBlockMergeFilter() 
    Merge.AddInput( reader.GetOutput() )
    Merge.Update() 
    
    Mb = vtk.vtkStructuredGridOutlineFilter()
    if vtk.VTK_MAJOR_VERSION <= 5:
      Mb.SetInput( Merge.GetOutput() );
    else:
      Mb.SetInputData( Merge.GetOutput() );

    Geom01 = vtk.vtkCompositeDataGeometryFilter();
    Geom01.SetInputConnection( Merge.GetOutputPort());

    self.output      = Geom01.GetOutput()
    self.vtk_pts     = self.output.GetPointData()

    self.vtk_n_pts   = self.output.GetNumberOfPoints()
    self.vtk_n_cells = self.output.GetNumberOfCells()
    self.vtk_n_props = self.vtk_pts.GetNumberOfArrays()

    self.reader = reader
    #get_range = self.output.GetScalarRange() 
    #print get_range
    
    print " [Vtk_unstructured_reader]:"
    print "                           Cells: %d" % (self.vtk_n_cells)
    print "                           Point: %d" % (self.vtk_n_pts)
    print "                           Props: ", 
    for i in range(self.vtk_n_props): 
      print "\'%s\'" % (self.vtk_pts.GetArrayName(i)),  
    print 


    print "\n+[Vtk_multiblock_reader]: \'%s\'" % file_name 


#--------------------------------------------------------------------------||--#
class Vtk_unstructured_reader:
  def __init__(self):
    self.vtk_pts   = None
    self.vtk_cells = None 
    self.prop_range = {}
    print "\n+[Vtk_unstructured_reader]:"


  def load_vtk_structured(self, file_name): 
    self.fname = file_name
    #reader = vtk.vtkUnstructuredGridReader()
    reader = vtk.vtkXMLUnstructuredGridReader()

    reader.SetFileName(file_name)
    reader.Update()

    #reader = vtk.vtkOpenFOAMReader() 
    #numberOfCellArrays = reader.GetNumberOfCellArrays();
    #for i in range(numberOfCellArrays):
    #  print reader.GetCellArrayName(i)

    self.output      = reader.GetOutput()
    self.vtk_pts     = self.output.GetPointData()

    self.vtk_n_pts   = self.output.GetNumberOfPoints()
    self.vtk_n_cells = self.output.GetNumberOfCells()
    self.vtk_n_props = self.vtk_pts.GetNumberOfArrays()

    self.reader = reader
    #get_range = self.output.GetScalarRange() 
    #print get_range
    
    print " [Vtk_unstructured_reader]:"
    print "                           Cells: %d" % (self.vtk_n_cells)
    print "                           Point: %d" % (self.vtk_n_pts)
    print "                           Props: ", 
    for i in range(self.vtk_n_props): 
      print "\'%s\'" % (self.vtk_pts.GetArrayName(i)),  
    print 


  def to_polydata(self):
    geometryFilter = vtk.vtkGeometryFilter()
    if vtk.VTK_MAJOR_VERSION <= 5:
      geometryFilter.SetInput( self.output );
    else:
      geometryFilter.SetInputData( self.output );
    geometryFilter.Update(); 

    return geometryFilter #.GetOutput();


  def get_points_prop_scalar(self, prop_name=""):    
    DATA = self.vtk_pts.GetArray(prop_name)
    self.n_cols = DATA.GetNumberOfComponents()
    self.n_rows = DATA.GetNumberOfTuples()
  
    data = []
    for j in range(self.n_rows): data.append( DATA.GetValue(j) ) 

    prop_range = DATA.GetRange() 
    self.prop_range[prop_name] = DATA.GetRange()

    print " [Vtk_unstructured_reader]:", 
    print "\'%s\'=[%d,%d]" %( prop_name, prop_range[0], prop_range[1] )

    return data


  def get_points_prop_scalar_range(self, prop_name=""): 
    return self.prop_range[prop_name]


  def set_points_prop_scalar(self, R, prop_name=""):
    array = vtk.vtkDoubleArray()
    array.SetName(prop_name)
    array.SetNumberOfComponents(1)
    array.SetNumberOfTuples( self.vtk_n_pts )

    for i in range(self.vtk_n_pts): array.SetValue(i, R[i])
    #self.reader.GetOutput().GetPointData().SetScalars(array) 
    self.reader.GetOutput().GetPointData().AddArray(array) 
    self.reader.Update() 


  def set_cell_prop_scalar(self, R, prop_name=""):
    array = vtk.vtkDoubleArray()
    array.SetName(prop_name)
    array.SetNumberOfComponents(1)
    array.SetNumberOfTuples( self.vtk_n_cells )

    for i in range(self.vtk_n_cells): array.SetValue(i, R[i])
    self.reader.GetOutput().GetCellData().AddArray(array)
    self.reader.Update() 


  def save(self, file_name):
    self.file_name = os.path.splitext(file_name)[0]
    #file_name_out  = self.file_name+".vtk"
    #writer = vtk.vtkUnstructuredGridWriter()
    #
    file_name_out  = self.file_name+".vtu"
    writer = vtk.vtkXMLUnstructuredGridWriter() 
    writer.SetFileName(file_name_out)

    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(self.reader.GetOutput())
    else:
      writer.SetInputData(self.reader.GetOutput())

    writer.Write()
    print "|_Vtk_writer.writer file: \'%s\' " % (file_name_out)


  def psave(self, file_name, mpi_size, mpi_range):
    self.file_name = os.path.splitext(file_name)[0] 

    file_name_out  = self.file_name+"_%d" % mpi_range + ".vtu"
    writer = vtk.vtkXMLUnstructuredGridWriter() 
    writer.SetFileName(file_name_out)
    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(self.reader.GetOutput())
    else:
      writer.SetInputData(self.reader.GetOutput())
    writer.Write()

    if(mpi_range==0): 
      pwriter = vtk.vtkXMLPUnstructuredGridWriter() 
      pwriter.SetFileName(self.file_name+".pvtu")
      pwriter.SetNumberOfPieces(mpi_size)

      if vtk.VTK_MAJOR_VERSION <= 5:
        pwriter.SetInput(self.reader.GetOutput())
      else:
        pwriter.SetInputData(self.reader.GetOutput())

      pwriter.Write()
      print "|_Vtk_writer.writer file: \'%s\' " % (file_name_out)



  def get_cell(self):
    vtk_data = self.output

    for idx in range(self.vtk_n_cells):  
      cell_type = vtk_data.GetCellType(idx) 
      cell_ids  = vtk_data.GetCell(idx)


  def save_data_coords(self, fname=""): 
    vtk_data = self.output

    fcoords = open(fname+"_COORDINATES.alya", "w")
    print>> fcoords, "COORDINATES"
    for idx in range(self.vtk_n_pts):
      print>> fcoords, idx+1, 
      pt = vtk_data.GetPoint(idx)
      n_pt = len(pt)
      for j in range(n_pt): print>> fcoords, pt[j],
      print>> fcoords, ""

    print>> fcoords, "END_COORDINATES"
    fcoords.close() 
    print "|_Vtk_writer.save_data_cell file: \'%s\' " % (fcoords.name)


  def save_data_cell(self, fname=""):
    vtk_data = self.output
    CelDic = Vtk2Alya_cell()

    ftype = open(fname+"_TYPES.alya", "w")
    felem = open(fname+"_ELEMENTS.alya", "w")
    print>> ftype, "TYPES"
    print>> felem, "ELEMENTS"

    alya_cells_found = {}  
    for idx in range(self.vtk_n_cells):
      cell_ids = vtk.vtkIdList()
      vtk_data.GetCellPoints(idx, cell_ids) 
      print>> felem, idx+1,
      n_cell_ids = cell_ids.GetNumberOfIds()  
      for j in range(n_cell_ids): print>> felem, cell_ids.GetId(j),
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

class Vtk_glyph:
  def __init__(self): 
    self.glyphs = vtk.vtkGlyph3D()
    self.glyph.ScalingOn()
    self.glyph.SetScaleModeToScaleByScalar()
    self.glyph.ClampingOn()  
    self.glyph.Update()
    
    #self.glyph.SetVectorModeToUseVector()
    #self.glyph.OrientOn()
#--------------------------------------------------------------------------||--#

#--------------------------------------------------------------------------||--#
class Vtk_render:
   def __init__(self):
     self.vtk_window = vtk.vtkRenderWindow()
     self.vtk_window.SetSize(1000, 750)
     self.vtk_window.SetWindowName("Alya View")
  
     self.vtk_interactor = vtk.vtkRenderWindowInteractor()
     self.vtk_interactor.SetRenderWindow(self.vtk_window)

     self.actors = []

     print "[Vtk_plot_multiview]"


   def set_actor(self, vtk_cosa, coloring_by = 'RTData'):
    #def set_actor(self, vtk_cosa, range_cosa, coloring_by = 'RTData'):
    #mapper = vtk.vtkPolyDataMapper()
    mapper = vtk.vtkHierarchicalPolyDataMapper()
    mapper.SetInputConnection( vtk_cosa.GetOutputPort() )
    mapper.SetScalarModeToUsePointFieldData() 
    mapper.SetColorModeToMapScalars()
    mapper.ScalarVisibilityOn()
    #mapper.SetScalarRange( range_cosa ) 
    mapper.SelectColorArray(coloring_by)
 
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetRepresentationToWireframe()
    self.actors.append( actor )

 
   def __set_vtk_renderer__(self):
      renderer = vtk.vtkRenderer()
      renderer.SetBackground( [0,0,0] )
      for actor in self.actors: renderer.AddActor(actor)
  
      renderer.ResetCamera()
      self.vtk_window.AddRenderer( renderer )
 
 
   def render(self):
     self.__set_vtk_renderer__() 
     self.vtk_window.Render()
     self.vtk_interactor.Start()
 
     print "|_Vtk_plot_multiview.render ok"
 
 
   def render_animation(self):
     self.vtk_window.Render()
     self.vtk_interactor.Start()
 
     print "|_Vtk_plot_multiview.render_animation ok"
#--------------------------------------------------------------------------||--#


#--------------------------------------------------------------------------||--#
class Vtk_structured_reader:
  def __init__(self):
    self.vtk_pts   = None
    self.vtk_cells = None 
    self.prop_range = {}
    print "\n+[Vtk_structured_reader]:"


  def load_vtk_structured(self, file_name): 
    self.fname = file_name
    #reader = vtk.vtkUnstructuredGridReader()
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(file_name)
    reader.Update()

    self.output      = reader.GetOutput()
    self.vtk_n_pts   = self.output.GetNumberOfPoints()
    self.vtk_pts_dim = self.output.GetDataDimension() 
    self.vtk_pts     = self.output.GetPointData()
    self.vtk_n_props = self.vtk_pts.GetNumberOfArrays()


    self.reader = reader
    #get_range = self.output.GetScalarRange() 
    #print get_range

    
    print " [Vtk_structured_reader]:", 
    for i in range(self.vtk_n_props): 
      print "\'%s\'" % (self.vtk_pts.GetArrayName(i)),  
    print 


  def to_polydata(self):
    geometryFilter = vtk.vtkGeometryFilter()
    if vtk.VTK_MAJOR_VERSION <= 5:
      geometryFilter.SetInput( self.output );
    else:
      geometryFilter.SetInputData( self.output );
    geometryFilter.Update(); 

    return geometryFilter #.GetOutput();


  def get_dimensions(self): 
    print " [Vtk_structured_reader]:", 
    print "\'%s\'= [%d]" %("GetNumberOfPoints=", self.output.GetNumberOfPoints()) 
    return self.output.GetDimensions() 	  


  def get_points(self):
    Nx, Ny, Nz = self.output.GetDimensions() 

    Pts = []
    l = 0
    for i in range(Nx): 
      for j in range(Ny):
        for k in range(Nz):
          Pts.append( self.output.GetPoint(l) )
          l=l+1

    return Pts  
  

  def get_points_prop_scalar(self, prop_name=""):    
    DATA = self.vtk_pts.GetArray(prop_name)
    self.n_cols = DATA.GetNumberOfComponents()
    self.n_rows = DATA.GetNumberOfTuples()

    Data = []
    k=0

    if(self.n_cols>1): 
      for j in range(self.n_rows): 
        data = []
        for i in range(self.n_cols): 
          data.append( DATA.GetValue(k) ) 
          k+=1
        Data.append( data ) 
    
    if(self.n_cols==1): 
      for j in range(self.n_rows): 
        for i in range(self.n_cols): 
          Data.append( DATA.GetValue(k) ) 
          k+=1

    prop_range = DATA.GetRange( self.n_cols-1 ) 
    self.prop_range[prop_name] = DATA.GetRange()

    print " [Vtk_structured_reader]:", 
    print "\'%s\'= [%f,%f]" %( prop_name, prop_range[0], prop_range[1] )

    return Data


  def get_points_prop_scalar_range(self, prop_name=""): 
    return self.prop_range[prop_name]


  def set_points_prop_scalar(self, R, prop_name=""):
    array = vtk.vtkDoubleArray()
    array.SetName(prop_name)
    array.SetNumberOfComponents(1)
    array.SetNumberOfTuples( self.n_rows )

    for i in range(self.n_rows): array.SetValue(i, R[i])
    self.reader.GetOutput().GetPointData().SetScalars(array) 
    #self.reader.update() 

    print " [Vtk_structured_reader]:", 
    print "\'%s\'= [%f,%f]" %( prop_name+" range", min(R), max(R) ) 


  def remove_points_prop_scalar(self, prop_name=""):
    self.vtk_pts.RemoveArray(prop_name)  


  def save(self, file_name):
    self.file_name = os.path.splitext(file_name)[0]
    file_name_out  = self.file_name+".vtk"

    writer = vtk.vtkStructuredGridWriter()
    writer.SetFileName(file_name_out)

    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(self.reader.GetOutput())
    else:
      writer.SetInputData(self.reader.GetOutput())

    writer.Write()
    print "|_Vtk_writer.writer file: \'%s\' " % (file_name_out)


  def save_raw(self, file_name, X, Y):
    self.file_name = os.path.splitext(file_name)[0]
    file_name_out  = self.file_name+".dat"
    
    fout = open(file_name_out, "w") 
    for x, y in zip(X, Y):
      print>> fout, x, y
    fout.close()     

## Example 
#Vts = Vtk_structured_reader() 
#Vts.load_vtk_structured( Fin ) 
#Pts  = Vts.get_points()
#VELs = Vts.get_points_prop_scalar('VELOC')
#Vts.remove_points_prop_scalar('VISCO') 
#--------------------------------------------------------------------------||--#


#--------------------------------------------------------------------------||--#
class Vtk_ensight_reader():
  def __init__(self):
    cdp = vtk.vtkCompositeDataPipeline()
    self.eg = vtk.vtkEnSightGoldReader()
    #self.eg = vtk.vtkGenericEnSightReader()
    self.eg.SetDefaultExecutivePrototype(cdp) 


  def init(self, filename=""):
    self.eg.SetCaseFileName(filename) 
    self.eg.ReadAllVariablesOn() 
    self.eg.Update()

    self.nVariables = self.eg.GetNumberOfVariables() 
    for i in range(self.nVariables):
      print "%d)" % (i+1), 
      print "\'%s\'" % self.eg.GetPointArrayName(i) 

    if(self.eg.GetOutput().GetNumberOfBlocks()==1): self.BLOCKid=0
    else: sys.exit() 

    self.set_block() 


  def set_block(self, idx=0): 
    self.BLOCKid=idx
    Vtu = self.eg.GetOutput().GetBlock(self.BLOCKid)
    vtu = vtk.vtkAppendFilter()
    vtu.SetInput( Vtu ) 
    vtu.Update() 
    self.vtu = vtu


  def get_vtu_output_port(self):
    return self.vtu.GetOutputPort()


  def get_data(self, property_name=""):
    Vtu = self.eg.GetOutput().GetBlock(self.BLOCKid)
    Datai = Vtu.GetPointData() 
    
    Propi = Datai.GetArray(property_name)
    nColsi = Propi.GetNumberOfComponents()
    nRowsi = Propi.GetNumberOfTuples() 
    print "\'%s\':" % Propi.GetName(),   
    print "(%d,%d)" % (nColsi, nRowsi)  


  def save_vtu(self, filename=""):
    writer = vtk.vtkXMLUnstructuredGridWriter() 
    writer.SetFileName(filename);
    writer.SetInput( self.vtu.GetOutput() )
    writer.Write()
#--------------------------------------------------------------------------||--#





#--------------------------------------------------------------------------||--#
class Vtk_openfoam_reader:
  def __init__(self):
    self.vtk_pts   = None
    self.vtk_cells = None 
    self.prop_range = {}
    print "\n+[Vtk_unstructured_reader]:"


  def load_vtk_structured(self, file_name): 
    self.fname = file_name
    reader = vtk.vtkOpenFOAMReader() 
    reader.SetFileName(file_name)
    reader.Update()

    self.output      = reader.GetOutput().GetBlock(0) # vtkUnstructuredGrid
    self.vtk_n_pts   = self.output.GetNumberOfPoints()
    self.vtk_pts     = self.output.GetPointData()
    self.vtk_n_cells = self.output.GetNumberOfCells()
    self.vtk_n_props = self.vtk_pts.GetNumberOfArrays()
    self.vtk_n_props = reader.GetNumberOfCellArrays();

    self.reader = reader
    
    print " [Vtk_unstructured_reader]:"
    print "                           Cells: %d" % (self.vtk_n_cells)
    print "                           Point: %d" % (self.vtk_n_pts)
    print "                           Props: ", 
    for i in range(self.vtk_n_props):
      print "\'%s\'" % (reader.GetCellArrayName(i)),  
    print 

    self.vtk_n_patch = reader.GetNumberOfPatchArrays();
    print "                           Patch: ", 
    for i in range(self.vtk_n_patch):
      print "                                \'%s\'" % (reader.GetPatchArrayName(i))   

    Vtu = self.output 
    vtu = vtk.vtkAppendFilter()
    vtu.SetInput( Vtu ) 
    vtu.Update() 

    writer = vtk.vtkXMLUnstructuredGridWriter() 
    writer.SetFileName("xxx.vtu");
    writer.SetInput( vtu )
    writer.Write()

#--------------------------------------------------------------------------||--#


#==============================================================================#
#==============================================================================#
