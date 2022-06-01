import vtk 

"""
if XXX.SetInputData(       YYY.GetOutput() )
if XXX.SetInputConnection( YYY.GetOutputPort() )
"""

PROPERTY = "TEMPE"

cdp    = vtk.vtkCompositeDataPipeline()
reader = vtk.vtkGenericEnSightReader() 
reader.SetDefaultExecutivePrototype(cdp)
reader.SetCaseFileName('runner.ensi.case')
reader.SetTimeValue(2)
#reader.Update()

reader.ReadAllVariablesOff()
reader.SetPointArrayStatus(PROPERTY,1) # 1=actived, 0=unactived 
reader.Update()



if(reader.GetOutput().GetNumberOfBlocks()==1): BLOCKi=0
DATAi = reader.GetOutput().GetBlock(BLOCKi)

print reader.GetOutputPort() 

set_range  = DATAi.GetScalarRange() 
set_center = DATAi.GetCenter() 
print set_range


print DATAi.GetNumberOfCells() 
print DATAi.GetNumberOfPoints()  
print DATAi.GetPoint(0)
PTS      = DATAi.GetPointData()
print PTS.GetScalars().GetRange()
PTS_PROP = PTS.GetArray(PROPERTY)
print PTS_PROP.GetValue(0)


for i in range(reader.GetNumberOfVariables()):
  key = reader.GetPointArrayName(i)
  print key, 
  print reader.GetPointArrayStatus(key) 
  #print reader.GetDescription(i) 
  #print reader.GetVariableType(i)
#
#reader.SetCellArrayStatus("thickness",1)
#print reader.GetTimeSets()
#print reader.GetParticleCoordinatesByIndex()	
#print reader.GetOutputInformation()


streamer = vtk.vtkStreamLine()
#streamer.SetInputData(DATAi) 

#-----------------------------------------------------------------------------#
geom = vtk.vtkGeometryFilter()
geom.SetInputConnection(reader.GetOutputPort())
#
calc = vtk.vtkArrayCalculator()
calc.SetInputConnection(geom.GetOutputPort())
calc.SetAttributeModeToUsePointData()
calc.SetFunction(PROPERTY+"*2.0") 
calc.AddScalarArrayName(PROPERTY,0)
calc.SetResultArrayName("test")


# windows 
mapper = vtk.vtkHierarchicalPolyDataMapper()
mapper.SetInputConnection(calc.GetOutputPort()) # geom.GetOutputPort()
mapper.SetColorModeToMapScalars()
mapper.SetScalarModeToUsePointFieldData()
mapper.ColorByArrayComponent(PROPERTY,0) #"test"
mapper.SetScalarRange(set_range)
#
actor = vtk.vtkActor();
actor.SetMapper(mapper)
#actor.GetProperty().SetColor(0,0,0)



# The outline gives context to the original data.
outline = vtk.vtkOutlineFilter() # vtkStructuredGridOutlineFilter
outline.SetInputConnection(calc.GetOutputPort())
outlineMapper = vtk.vtkPolyDataMapper()
outlineMapper.SetInputConnection(outline.GetOutputPort())
#outlineMapper.SetInput(outline)
outlineActor = vtk.vtkActor()
outlineActor.SetMapper(outlineMapper)
outlineActor.GetProperty().SetColor(0, 0, 0)


#-------------------------------------- rendering 
renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
renderer.AddActor(outlineActor)

renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindow.StereoCapableWindowOn()
renderWindow.SetSize(1200, 600)

renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

renderWindow.Render()
renderWindowInteractor.Start()

reader.SetDefaultExecutivePrototype(None)
print "OK! \n"
#-----------------------------------------------------------------------------#

