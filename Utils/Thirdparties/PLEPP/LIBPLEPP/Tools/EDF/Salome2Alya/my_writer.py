import numpy as np 
import vtk 
import os
import sys 


class VTK_writer:
  def __init__(self):
    self.obj       = None 
    self.num_pts   = -1
    self.vtk_pts   = vtk.vtkPoints() 
    self.vtk_cells = vtk.vtkCellArray()
    self.vtk_tetra = vtk.vtkTetra()

    self.vtk_props = [] 
    self.vtk_cells_type = []
    self.vtk_cells_id = []

    self.vtk_types = []    #(vtk_type, vtk_id, vtkCell) <= linear cell types found in VTK 
    self.vtk_types.append( ("VTK_TRIANGLE",     5, vtk.vtkTriangle()   ) )
    self.vtk_types.append( ("VTK_QUAD",         9, vtk.vtkQuad()       ) )
    self.vtk_types.append( ("VTK_TETRA",       10, vtk.vtkTetra()      ) )
    self.vtk_types.append( ("VTK_PYRAMID",     14, vtk.vtkPyramid() ) )
    self.vtk_types.append( ("VTK_WEDGE",       13, vtk.vtkWedge()       ) ) 
    self.vtk_types.append( ("VTK_HEXAHEDRON",  12, vtk.vtkHexahedron() ) )

    print "\n+Vtk_writer"
    
    
  def writer(self, file_name): 
    self.file_name = os.path.splitext(file_name)[0]
    file_name_out  = self.file_name+".vtk"

    writer = vtk.vtkUnstructuredGridWriter()
    writer.SetFileName(file_name_out)

    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(self.obj) 
    else:
      writer.SetInputData(self.obj)
    #writer.Write() 
    #print "|_Vtk_writer.writer file: \'%s\' " % (file_name_out)
    self.writer_xml(file_name) 


  def writer_xml(self, file_name):
    self.file_name = os.path.splitext(file_name)[0]
    file_name_out  = self.file_name+".vtu"

    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(file_name_out)

    if vtk.VTK_MAJOR_VERSION <= 5:
      writer.SetInput(self.obj)
    else:
      writer.SetInputData(self.obj)
    writer.Write()
    print "|_Vtk_writer.writer_xml file: \'%s\' " % (file_name_out)


  def set_point(self, coord):
    self.vtk_pts.InsertNextPoint( coord[0], coord[1], coord[2] )


  def set_points_prop(self, R, prop_name=""):    
    array = vtk.vtkDoubleArray()
    array.SetName(prop_name)
    array.SetNumberOfComponents(1)
    array.SetNumberOfTuples( self.num_pts )

    for i in range(self.num_pts): array.SetValue(i, R[i])
    #self.obj.GetPointData().AddArray(array)
    #self.obj.Update()

    self.vtk_props.append( array ) 
    print "|_Vtk_writer.set_vertices_prop: \'%s\' " % (prop_name)


  def set_points_prop_vec(self, R, prop_name=""):
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


  def set_cell(self, nodes, vtk_cell):
    self.__set_cell__(nodes, vtk_cell[2])


  def __set_cell__(self, nodes, vtk_cell):
    i = 0
    vtk_list = vtk.vtkIdList() 
    for node in nodes: 
      #vtk_cell.GetPointIds().SetId(i, node-1)
      vtk_list.InsertNextId(node-1) 
      i += 1
    #self.vtk_cells.InsertNextCell(vtk_cell)
    self.vtk_cells_type.append( vtk_cell.GetCellType() ) 
    self.vtk_cells_id.append(vtk_list) 


  def to_vtkunstructured_grid(self): 
    vtk_data = vtk.vtkUnstructuredGrid()
    vtk_data.SetPoints(self.vtk_pts) 
    
    for cell_type, cell_ids in zip(self.vtk_cells_type, self.vtk_cells_id): 
      vtk_data.InsertNextCell(cell_type, cell_ids) 

    for array in self.vtk_props: 
      vtk_data.GetPointData().AddArray(array)
    
    #vtk_data.SetCells(self.vtk_cells)        
    #if vtk.VTK_MAJOR_VERSION <= 5:
    #vtk_data.Update()
    
    self.obj = vtk_data 
    print "|_Vtk_writer.to_vtkunstructured_grid "


  def get_tetra_vol(self, nodes, pts):
    a = np.array( pts[nodes[0]] ) #  \ 
    b = np.array( pts[nodes[1]] ) #  |__BASE anti-clockwise 
    c = np.array( pts[nodes[2]] ) # _|
    d = np.array( pts[nodes[3]] )
    
    # (b-a)x(c-a).(d-a)    
    n = np.cross(b-a, c-a)
    V = np.dot(d-a, n)/6.0
  
    # abd 
    n  = np.cross(b-a, d-a)
    A1 = np.dot(n, n)*0.5

    # adc
    n  = np.cross(d-a, c-a)
    A2 = np.dot(n, n)*0.5

    # acb
    n  = np.cross(b-a, c-a)
    A3 = np.dot(n, n)*0.5

    # bcd
    n  = np.cross(d-b, c-b)
    A4 = np.dot(n, n)*0.5

    ok = np.array([V, A1, A2, A3, A4]) 
    if((ok<0).any()): 
      print "ERROR (Negative value):", ok 
      print a, b, c, d
      print 
      sys.exit()

    return [V, A1, A2, A3, A4]
    
    



#============================================================================# 
