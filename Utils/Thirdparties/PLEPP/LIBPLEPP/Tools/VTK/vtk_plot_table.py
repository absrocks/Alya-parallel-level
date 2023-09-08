import vtk 
import numpy as np

"""
vtk V5.8.0 ok!
"""

class Vtk_plot_table:
  def __init__(self):
    self.i = 0 
    self.props_name = {} 

    self.table = vtk.vtkTable() 
    self.view  = vtk.vtkContextView() 
    self.view.GetRenderer().SetBackground(1.0, 1.0, 1.0)
    self.view.GetRenderWindow().SetSize(600, 600)

    print "+Vtk_plot \'%s\'" % vtk.vtkVersion().GetVTKVersion() 


  def set_pts_size(self, numPoints=-1):
    self.numPoints = numPoints
    self.table.SetNumberOfRows(self.numPoints) 
    print "+Vtk_plot.set_pts_size %d" % numPoints


  def set_prop(self, key): 
    arrC = vtk.vtkFloatArray() 
    arrC.SetName(key)
    self.table.AddColumn(arrC)

    self.props_name[key] = self.i    
    self.i += 1

    print "|_Vtk_plot.set_prop \'%s\'" % key


  def set_vals(self, vals, key):  
    col = self.props_name.get(key,-1)
    if(col==-1): 
      print "ERROR: there not exist \'key\' \n\n!!" % key
      sys.exit()

    #print vals 

    self.table.SetNumberOfRows(self.numPoints) 
    for i in range(self.numPoints): 
      self.table.SetValue(i, col, vals[i]) 
      #self.table.GetColumn(col).InsertNextValue( vals[i] )

    print "|_Vtk_plot.set_vals \'%s\'" % key
    

  def plot(self): 
    for (key, val) in self.props_name.items(): 
      if(val!=0): 
        print "|_Vtk_plot.plot Y:\'%s\'" % key
        chart = vtk.vtkChartXY()
        self.view.GetScene().AddItem(chart);
    
        line = chart.AddPlot(vtk.vtkChart.LINE)
        line.SetInput(self.table, 0, val);
        line.SetColor(0, 255, 0, 255);
        line.SetWidth(5.0);
        line.GetPen().SetLineType(1) 
      else:
        print "|_Vtk_plot.plot X:\'%s\'" % key

    self.view.GetRenderWindow().SetMultiSamples(0);
    self.view.GetInteractor().Initialize()
    self.view.GetInteractor().Start()

    print "|_Vtk_plot.plot ok!" 


## Example #################################################
"""
numPoints = 69
inc = 10.0 / (numPoints-1)
X = []
Y = [] 
Z = []
for i in range(numPoints):
  X.append(i*inc)
  Y.append(np.cos(i * inc))
  Z.append(np.sin(i * inc))

P = Vtk_plot_table()
P.set_pts_size(numPoints) 
P.set_prop("X")
P.set_prop("Y01")
P.set_prop("Y02")
P.set_vals(X, "X")
P.set_vals(Y, "Y02")
P.set_vals(Z, "Y01")
P.plot() 
"""

