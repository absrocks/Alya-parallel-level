import vtk 


class Vtk_plot_contour:
  def __init__(self):
    self.contour = vtk.vtkContourFilter()
    self.mapper  = vtk.vtkHierarchicalPolyDataMapper()
    self.actor   = vtk.vtkActor()
    self.values  = [] 

    print "+Vtk_plot_contour"


  def set_value(self, value): 
    self.values.append( value )

    print "|_Vtk_plot_contour.set_value", value 


  def set_connection(self, obj):  
    self.contour.SetInputConnection(obj.GetOutputPort())
    for i in range(len(self.values)): self.contour.SetValue(i, self.values[i])

    self.mapper.SetInputConnection(self.contour.GetOutputPort())
    self.mapper.SetImmediateModeRendering(1)
    self.mapper.SetScalarRange(0,1)
    self.mapper.SetScalarVisibility(1)

    self.actor.SetMapper(self.mapper)
    self.actor.GetProperty().SetInterpolationToGouraud()
    #self.actor.GetProperty().SetRepresentationToWireframe()

    print "|_Vtk_plot_contour.set_connection"


  def get_actor(self):
    print "|_Vtk_plot_contour.get_actor"
    return self.actor


## contour
#Contour0 = vtk.vtkContourFilter()
#Contour0.SetInputConnection(reader.GetOutputPort())
#Contour0.SetValue(0, set_average)
#
#contour_mapper = vtk.vtkHierarchicalPolyDataMapper()
#contour_mapper.SetInputConnection(Contour0.GetOutputPort())
#contour_mapper.SetImmediateModeRendering(1)
#contour_mapper.SetScalarRange(0,1)
#contour_mapper.SetScalarVisibility(1)
#
#contour_actor = vtk.vtkActor()
#contour_actor.SetMapper(contour_mapper)
#contour_actor.GetProperty().SetInterpolationToGouraud()
#contour_actor.GetProperty().SetRepresentationToWireframe()
#contour_actor.GetProperty().SetRepresentationToSurface()


