import vtk 
import os


class Vtk_plot_linesource:
  def __init__(self):
    self.source = vtk.vtkLineSource()
    self.probe  = vtk.vtkProbeFilter()
    print "+Vtk_plot_linesource"


  def set_source(self, init, end, resolution, ranges):
    self.source.SetPoint1(init)
    self.source.SetPoint2(end)
    self.source.SetResolution(resolution)
    
    self.range = ranges

    print "|_Vtk_plot_linesource.set_source"


  def set_connection(self, obj):
    self.probe.SetInputConnection(self.source.GetOutputPort())
    self.probe.SetSource(obj)
    self.probe.Update()
 
    self.data = self.probe.GetOutput()
     
    print "|_Vtk_plot_linesource.set_connection"


  def get_actor(self, diametro=0.001):
    tube = vtk.vtkTubeFilter()
    tube.SetInput(self.probe.GetPolyDataOutput())
    tube.SetNumberOfSides(5)
    tube.SetRadius(diametro)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())
    mapper.SetScalarRange(self.range)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    print "|_Vtk_plot_linesource.get_actor"
    return actor


  def get_data(self, key, step=0):
    #Y = [] 
    #for i in range(self.get_npts()): 
    #  Y.append(self.data.GetPointData().GetArray(key).GetValue(i))

    PTS_PROP = self.data.GetPointData().GetArray(key)
    NCOLS    = PTS_PROP.GetNumberOfComponents() 
    #for row in range(PTS_PROP.GetNumberOfTuples()):
    #  for col in range(NCOLS):
    #    print PTS_PROP.GetValue(NCOLS*row+col),  
    #  print 

    Y = [] 
    if(NCOLS==1): step = 0 
    for row in range(PTS_PROP.GetNumberOfTuples()):
      #print PTS_PROP.GetValue(NCOLS*row+step)  
      Y.append( PTS_PROP.GetValue(NCOLS*row+step) ) 

    print "|_Vtk_plot_linesource.get_data", key
    return Y 


  def get_npts(self):
    print "|_Vtk_plot_linesource.get_npts"
    return self.data.GetNumberOfPoints() 


  def save_data(self, file_name, key):
    file_name  = os.path.splitext(file_name)[0]    
    file_name += "_"+ key +".dat"
    
    #file_out = open(file_name, "w")
    #num = 1.0/(self.get_npts()-1.0) 
    #for i in range(self.get_npts()): 
    #  print >> file_out, i*num, 
    #  print >> file_out, self.data.GetPointData().GetArray(key).GetValue(i)
    #file_out.close() 


    PTS_PROP = self.data.GetPointData().GetArray(key)
    NCOLS    = PTS_PROP.GetNumberOfComponents()

    file_out = open(file_name, "w")
    num = 1.0/(self.get_npts()-1.0)

    for row in range(PTS_PROP.GetNumberOfTuples()):
      print >> file_out, row*num, 
      for col in range(NCOLS):
        print >> file_out, PTS_PROP.GetValue(NCOLS*row+col),  
      print >> file_out, "" 
    file_out.close() 


# line 
#probeLine = vtk.vtkLineSource()
#probeLine.SetPoint1(0.0,0.0, 0.0)
#probeLine.SetPoint2(0.0,0.0,30.0)
#probeLine.SetResolution(100)
#probe = vtk.vtkProbeFilter()
#probe.SetInputConnection(probeLine.GetOutputPort())
#probe.SetSource(DATAi)
#probe.Update()
#
#probeTube = vtk.vtkTubeFilter()
#probeTube.SetInput(probe.GetPolyDataOutput())
#probeTube.SetNumberOfSides(5)
#probeTube.SetRadius(.05)
#probeMapper = vtk.vtkPolyDataMapper()
#probeMapper.SetInputConnection(probeTube.GetOutputPort())
#probeMapper.SetScalarRange(set_range)
#probeActor = vtk.vtkActor()
#probeActor.SetMapper(probeMapper)


#from vtk_plot_table import Vtk_plot_table
#X = []
#Y = [] 

#numPoints = probe.GetOutput().GetNumberOfPoints() 
#for i in range(numPoints): X.append(i)
#for i in range(numPoints): Y.append(probe.GetOutput().GetPointData().GetArray(PROPERTY).GetValue(i))

