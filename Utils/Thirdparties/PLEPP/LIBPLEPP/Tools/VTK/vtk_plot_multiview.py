import vtk 
import numpy as np


#===================================================================| png |===#
class VtkAxis: 
  def __init__(self):
    axes   = vtk.vtkAxesActor() 
    self.widget = vtk.vtkOrientationMarkerWidget()
    self.widget.SetOutlineColor( 0.9300, 0.5700, 0.1300 )
    self.widget.SetOrientationMarker( axes )

  def set_actor(self, render): 
    self.widget.SetInteractor( render )
    self.widget.SetViewport( 0.0, 0.0, 0.2, 0.2 )
    self.widget.SetEnabled( 1 )
    #self.widget.InteractiveOn();


#===================================================================| png |===#
class VtkPNGWriter():
  def set_actor(self, render):
    self.render = render 
    
  def save(self, file_name=""):
    windowToImageFilter = vtk.vtkWindowToImageFilter()
    windowToImageFilter.SetInput(self.render);
    windowToImageFilter.SetMagnification(3) 
    windowToImageFilter.SetInputBufferTypeToRGBA() 
    windowToImageFilter.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(file_name) 
    writer.SetInputConnection(windowToImageFilter.GetOutputPort())
    writer.Write()
    
    print "|_save png %s" % file_name 


#===============================================================|animation|===#
class vtkTimerCallback():
  def __init__(self):
    self.timer_count = 0
    self.actor       = None 

  def set_actor(self, actor):
    self.actor = actor

  def execute(self, obj, event):
    print self.timer_count
    
    theta = self.timer_count/360.0*np.pi
    self.actor.SetPosition(np.sin(theta), np.cos(theta),0)

    iren = obj
    iren.GetRenderWindow().Render()
    self.timer_count += 1


#=========================================================|interactor style|===#
class Vtk_interactor_style(vtk.vtkInteractorStyleTrackballCamera):
  def __init__(self,parent=None):                
    self.AddObserver("LeftButtonPressEvent",   self.leftPressEvent)
    self.AddObserver("LeftButtonReleaseEvent", self.leftReleaseEvent)
    self.AddObserver("RightButtonPressEvent",  self.execute)
    print "+Vtk_interactor_style"

    self.timer_count = 0
    self.actor       = None 
 
  def set_actor(self, actor):
    self.actor = actor
    print "+Vtk_interactor_style.up"  

 
  def leftPressEvent(self, obj, event):
    self.OnLeftButtonDown() # OnLeftButtonDown Left/Right/Middle

    print self.actor.GetPosition()

    print "|_Vtk_interactor_style.down"
    return
 
  def leftReleaseEvent(self, obj, event):
    self.OnLeftButtonUp() 
    print "|_Vtk_interactor_style.up"
    return


  def execute(self, obj, event):
    print self.timer_count
    #self.actor.SetPosition(self.timer_count, self.timer_count,0)

    iren = obj
    iren.GetRenderWindow().Render()
    self.timer_count += 1


#===============================================================|multiview|===#
class Vtk_plot_multiview:
  def __init__(self):
    # rw 
    self.vtk_window = vtk.vtkRenderWindow()
    self.vtk_window.SetSize(1000, 750)
    self.vtk_window.SetWindowName("Alya View")

    self.vtk_interactor_style = Vtk_interactor_style()
    self.vtk_timer_callback   = vtkTimerCallback()

    # iren 
    self.vtk_interactor = vtk.vtkRenderWindowInteractor()
    self.vtk_interactor.SetRenderWindow(self.vtk_window)
    self.vtk_interactor.SetInteractorStyle( self.vtk_interactor_style ) 
    self.vtk_interactor.AddObserver('TimerEvent', self.vtk_timer_callback.execute)
    #self.vtk_interactor.CreateRepeatingTimer(100) 

    self.vtk_png = VtkPNGWriter()
    self.vtk_png.set_actor(self.vtk_window)
    
    self.vtk_axes = VtkAxis() 

    self.color_background = []
    self.color_background.append( (0.8,0.8,0.9) ) 
    self.color_background.append( (0.8,0.8,0.8) ) 
    self.color_background.append( (0.8,0.8,0.7) ) 
    
    x_middle = 0.75
    self.range_background = [] 
    self.range_background.append( [0.0,0.0,x_middle,1.0] ) 
    self.range_background.append( [x_middle,0.5,1.0,1.0] ) 
    self.range_background.append( [x_middle,0.0,1.0,0.5] ) 
    
    self.camara_position = [] 
    self.camara_position.append( [0.0,0.0,5.0] ) 
    self.camara_position.append( [0.0,0.0,5.0] ) 
    self.camara_position.append( [0.0,0.0,5.0] ) 

    self.actor2D = None
    print "+Vtk_plot_multiview"

  
  def set_camara_position(self, camara_pos):
    self.camara_position = camara_pos 
  

  def set_actors(self, actors): 
    for i in range(3):
      ranges = self.range_background[i]
      color  = self.color_background[i]
      actor  = actors[i]
      camara_pos = self.camara_position[i]
      self.__set_vtk_renderer__(ranges, color, actor, camara_pos)

    actor = actors[0]
    self.vtk_interactor_style.set_actor(actor[0])
    self.vtk_timer_callback.set_actor(actor[0])
    print "|_Vtk_plot_multiview.set_actors"


  def set_actor2D(self, actor):
    self.actor2D = actor


  def __set_vtk_renderer__(self, ranges, color, actors, camara_position):
    render = vtk.vtkRenderer()
    render.SetViewport(ranges)
    render.SetBackground(color)
    for actor in actors: render.AddActor(actor)
    if(self.actor2D!=None and len(actors)!=1): 
      render.AddActor2D(self.actor2D)
    
    render.ResetCamera()
    #render.GetActiveCamera().Zoom(1.25)      
    render.GetActiveCamera().SetPosition(camara_position)
   
    self.vtk_window.AddRenderer(render)


  def render(self, file_name="xxx_multiview.png"):
    #camera = vtk.vtkCamera() 
    #camera.SetPosition(3,0,0)
    #camera.SetFocalPoint(0,0,0)
    
    self.vtk_axes.set_actor(self.vtk_interactor) 
    self.vtk_window.Render()
    self.vtk_interactor.Start()
    
    self.vtk_png.save(file_name)
    
    print "|_Vtk_plot_multiview.render ok"


  def render_animation(self):
    self.vtk_window.Render()
    self.vtk_interactor.Initialize()
    #self.vtk_interactor.AddObserver('TimerEvent', self.vtk_timer_callback.execute)
    self.vtk_interactor.CreateRepeatingTimer(100) 
    self.vtk_interactor.Start()

    print "|_Vtk_plot_multiview.render_animation ok"


#===============================================================|example|===#
## EXAMPLE
"""
sphereSource = [] 
sphereSource.append( vtk.vtkSphereSource() ) 
sphereSource.append( vtk.vtkSphereSource() ) 
sphereSource.append( vtk.vtkSphereSource() ) 
    
sphereSource[0].SetCenter(0.0, 0.0, 0.0)
sphereSource[0].SetRadius(0.5)
sphereSource[1].SetCenter(0.0, 0.0, 0.0)
sphereSource[1].SetRadius(2.5)
sphereSource[2].SetCenter(0.0, 0.0, 0.0)
sphereSource[2].SetRadius(5.0)

color = []
color.append( (1.0,0.0,0.0) ) 
color.append( (0.0,1.0,0.0) ) 
color.append( (0.0,0.0,1.0) ) 

actors = [] 
for i in range(3):
  mapper = vtk.vtkPolyDataMapper()
  mapper.SetInputConnection(sphereSource[i].GetOutputPort())

  actor = vtk.vtkActor()
  actor.SetMapper(mapper)
  actor.GetProperty().SetColor(color[i])
  actor.GetProperty().SetRepresentationToWireframe()
  actors.append( [actor] ) 

MV = Vtk_plot_multiview()
MV.set_actors(actors)
MV.render_animation()
"""

