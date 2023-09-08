import vtk 
import numpy as np


class Vtk_plot:
  def __init__(self):
    self.i = 0 
    self.props_name = {} 

    self.table = vtk.vtkTable() 
    self.view  = vtk.vtkContextView() 
    self.view.GetRenderer().SetBackground(1.0, 1.0, 1.0);

    print "+Vtk_plot"


  def set_pts_size(self, numPoints=-1):
    self.numPoints = numPoints
    self.table.SetNumberOfRows(self.numPoints) 
    print "+Vtk_plot", numPoints


  def set_prop(self, key): 
    arrC = vtk.vtkFloatArray() 
    arrC.SetName(key)
    self.table.AddColumn(arrC)

    self.props_name[key] = self.i    
    self.i += 1


  def set_vals(self, vals, key):  
    col = self.props_name.get(key,-1)
    if(col==-1): 
      print "ERROR: there not exist \'key\' \n\n!!" % key
      sys.exit()

    self.table.SetNumberOfRows(self.numPoints) 
    for i in range(self.numPoints): self.table.SetValue(i, col, vals[i]) 

    print "+Vtk_plot \'%s\'" % key
    

  def plot(self): 
    for (key, val) in self.props_name.items(): 
      print key, val 
    
      chart = vtk.vtkChartXY()
      self.view.GetScene().AddItem(chart);
    
      line = chart.AddPlot(vtk.vtkChart.LINE)
      line.SetInput(self.table, 0, 1);
      line.SetColor(0, 255, 0, 255);
      line.SetWidth(5.0);
      line.GetPen().SetLineType(1) 

    self.view.GetRenderWindow().SetMultiSamples(0);
    self.view.GetInteractor().Initialize()
    self.view.GetInteractor().Start()

    print "+Vtk_plot \'%s\'"



numPoints = 69
inc = 7.5 / (numPoints-1)
X = []
Y = [] 
for i in range(numPoints):
  #table.SetValue(i, 0, i * inc);
  #table.SetValue(i, 1, np.cos(i * inc));
  X.append(i*inc)
  Y.append(np.cos(i * inc))


P = Vtk_plot()
P.set_pts_size(numPoints) 
P.set_prop("X")
P.set_prop("Y")
P.set_vals(X, "X")
P.set_vals(Y, "Y")
P.plot() 

#P.view.GetRenderWindow().SetMultiSamples(0);
#P.view.GetInteractor().Initialize()
#P.view.GetInteractor().Start()
print "OK! \n"


table = vtk.vtkTable()
 
arrX = vtk.vtkFloatArray() 
arrX.SetName("X Axis")

arrC = vtk.vtkFloatArray() 
arrC.SetName("Cosine")

table.AddColumn(arrX)
table.AddColumn(arrC)
 
numPoints = 69
inc = 7.5 / (numPoints-1)
table.SetNumberOfRows(numPoints);
for i in range(numPoints):
  table.SetValue(i, 0, i * inc);
  table.SetValue(i, 1, np.cos(i * inc));

view = vtk.vtkContextView() 
view.GetRenderer().SetBackground(1.0, 1.0, 1.0);
 
chart = vtk.vtkChartXY()
view.GetScene().AddItem(chart);
line = chart.AddPlot(vtk.vtkChart.LINE)
line.SetInput(table, 0, 1);
line.SetColor(0, 255, 0, 255);
line.SetWidth(5.0);
line.GetPen().SetLineType(1) 

#line = chart.AddPlot(vtk.vtkChart.LINE);
#line.SetInput(table, 0, 2);
#line.SetColor(255, 0, 0, 255);
#line.SetWidth(5.0);
#line.GetPen().SetLineType(2) 

#view.GetRenderWindow().SetMultiSamples(0);
#view.GetInteractor().Initialize()
#view.GetInteractor().Start()




