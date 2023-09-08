import vtk 


class Vtk_plot_cutter:
  def __init__(self):
    self.plane  = vtk.vtkPlane()
    self.cutter = vtk.vtkCutter()    
    self.mapper = vtk.vtkHierarchicalPolyDataMapper() 
    self.actor  = vtk.vtkActor()
    
    self.range = None
    
    print "+Vtk_plot_cutter"


  def set_plane(self, center, normal, ranges):
    self.plane.SetOrigin(center)
    self.plane.SetNormal(normal)
    #self.cutter.SetScalarRange(ranges)
    self.range = ranges 

    print "|_Vtk_plot_cutter.set_plane"


  def set_connection(self, obj):
    self.cutter.SetInputConnection(obj.GetOutputPort())
    self.cutter.SetCutFunction(self.plane)
    self.mapper.SetInputConnection(self.cutter.GetOutputPort())
    self.mapper.SetScalarRange(self.range)
    self.actor.SetMapper(self.mapper)

    print "|_Vtk_plot_cutter.set_connection"


  def get_actor(self):
    print "|_Vtk_plot_cutter.get_actor"
    return self.actor


## cutplane
#cutPlane = vtk.vtkPlane()
#cutPlane.SetOrigin(set_center)
#cutPlane.SetNormal(1,0,0)
#
#planeCut = vtk.vtkCutter()
#planeCut.SetInputConnection(reader.GetOutputPort())#SetInput(toRectilinearGrid.GetRectilinearGridOutput())
#planeCut.SetCutFunction(cutPlane)
#
#cut_mapper = vtk.vtkHierarchicalPolyDataMapper() #vtkDataSetMapper()
#cut_mapper.SetInputConnection(planeCut.GetOutputPort())
#cut_mapper.SetScalarRange(set_range)
#
#cut_actor = vtk.vtkActor()
#cut_actor.SetMapper(cut_mapper)

