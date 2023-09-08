import vtk 
import sys 
"""
if XXX.SetInputData(       YYY.GetOutput() )
if XXX.SetInputConnection( YYY.GetOutputPort() )
"""


PROPERTY = "PRESS"

cdp    = vtk.vtkCompositeDataPipeline()
reader = vtk.vtkGenericEnSightReader() 
reader.SetDefaultExecutivePrototype(cdp)
reader.SetCaseFileName('runner.ensi.case')
reader.SetTimeValue(10)
#reader.Update()

reader.ReadAllVariablesOff()
reader.SetPointArrayStatus(PROPERTY,1) # 1=actived, 0=unactived 
reader.Update()


if(reader.GetOutput().GetNumberOfBlocks()==1): BLOCKi=0
else: sys.exit()
DATAi = reader.GetOutput().GetBlock(BLOCKi)

print "(1)", reader.GetOutputPort() 

set_range  = DATAi.GetScalarRange() 
set_center = DATAi.GetCenter() 
set_average = 0.5*(set_range[0]+set_range[1])
print "(2)", set_range


print DATAi.GetNumberOfCells() 
print DATAi.GetNumberOfPoints()  
print DATAi.GetPoint(0)
PTS      = DATAi.GetPointData()
#print PTS.GetScalars().GetRange()
PTS_PROP = PTS.GetArray(PROPERTY)
print PTS_PROP.GetValue(0)


for i in range(reader.GetNumberOfVariables()):
  key = reader.GetPointArrayName(i)
  print key, 
  print reader.GetPointArrayStatus(key) 


#---------------------------------------------------------------| filters |----# 
## reader 
reader_geom = vtk.vtkGeometryFilter()
reader_geom.SetInputConnection(reader.GetOutputPort())
#
reader_mapper = vtk.vtkHierarchicalPolyDataMapper()
reader_mapper.SetInputConnection(reader_geom.GetOutputPort())
reader_mapper.SetScalarRange(set_range)
#
reader_actor = vtk.vtkActor();
reader_actor.SetMapper(reader_mapper)
reader_actor.GetProperty().SetRepresentationToWireframe()

## contour
Contour0 = vtk.vtkContourFilter()
Contour0.SetInputConnection(reader.GetOutputPort())
Contour0.SetValue(0,set_average)
#
contour_mapper = vtk.vtkHierarchicalPolyDataMapper()
contour_mapper.SetInputConnection(Contour0.GetOutputPort())
contour_mapper.SetImmediateModeRendering(1)
contour_mapper.SetScalarRange(0,1)
contour_mapper.SetScalarVisibility(1)
#
contour_actor = vtk.vtkActor()
contour_actor.SetMapper(contour_mapper)
contour_actor.GetProperty().SetInterpolationToGouraud()
contour_actor.GetProperty().SetRepresentationToWireframe()
#contour_actor.GetProperty().SetRepresentationToSurface()

## outline 
toRectilinearGrid = vtk.vtkCastToConcrete()
toRectilinearGrid.SetInput(DATAi)
toRectilinearGrid.Update()
#
outline = vtk.vtkStructuredGridOutlineFilter()
#outline.SetInputConnection(reader.GetOutputPort())
#outline.SetInputData(reader.GetOutput().GetBlock(0))
#outline.SetInput(reader.GetOutput().GetBlock(0))
outline.SetInput(toRectilinearGrid.GetRectilinearGridOutput())
#
outline_mapper = vtk.vtkHierarchicalPolyDataMapper()
outline_mapper.SetInputConnection(outline.GetOutputPort())
#
outline_actor = vtk.vtkActor()
outline_actor.SetMapper(outline_mapper)

## surface 
dss = vtk.vtkDataSetSurfaceFilter()
dss.SetInputConnection(reader.GetOutputPort())
#
dss_mapper = vtk.vtkHierarchicalPolyDataMapper()
dss_mapper.SetInputConnection(dss.GetOutputPort())
dss_mapper.SetColorModeToMapScalars()
#dss_mapper.SetScalarModeToUseCellFieldData()
dss_mapper.ColorByArrayComponent( 1*PROPERTY,0)
dss_mapper.ColorByArrayComponent(10*PROPERTY,1)
dss_mapper.SetScalarRange(set_range)
#
dss_actor = vtk.vtkActor()
dss_actor.SetMapper(dss_mapper)
dss_actor.GetProperty().SetRepresentationToWireframe()

## cutplane
cutPlane = vtk.vtkPlane()
cutPlane.SetOrigin(set_center)
cutPlane.SetNormal(1,0,0)
#
planeCut = vtk.vtkCutter()
planeCut.SetInputConnection(reader.GetOutputPort())#SetInput(toRectilinearGrid.GetRectilinearGridOutput())
planeCut.SetCutFunction(cutPlane)
#
cut_mapper = vtk.vtkHierarchicalPolyDataMapper() #vtkDataSetMapper()
cut_mapper.SetInputConnection(planeCut.GetOutputPort())
cut_mapper.SetScalarRange(set_range)
#
cut_actor = vtk.vtkActor()
cut_actor.SetMapper(cut_mapper)


## contour 
iso = vtk.vtkContourFilter()
iso.SetInputConnection(reader.GetOutputPort()) #XXX.SetInput(toRectilinearGrid.GetRectilinearGridOutput())
iso.SetValue(0, set_average)

normals = vtk.vtkPolyDataNormals() 
normals.SetInputConnection(iso.GetOutputPort())
normals.SetFeatureAngle(45)

isoMapper = vtk.vtkHierarchicalPolyDataMapper() # XXX.vtkPolyDataMapper()
isoMapper.SetInputConnection(normals.GetOutputPort())
isoMapper.ScalarVisibilityOff()

isoActor = vtk.vtkActor()
isoActor.SetMapper(isoMapper)
#isoActor.GetProperty().SetColor(GetRGBColor('bisque'))
isoActor.GetProperty().SetRepresentationToWireframe()


#--------------------------------------------------------------| rendering |----# 
renderer = vtk.vtkRenderer()
renderer.SetBackground((0.7,0.7,0.7))
renderer.AddActor(reader_actor)
renderer.AddActor(contour_actor)
renderer.AddActor(dss_actor)
#renderer.AddActor(outline_actor) # segmentation fault 
renderer.AddActor(cut_actor)
#renderer.AddActor(isoActor)


renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindow.StereoCapableWindowOn()
renderWindow.SetSize(600, 800)

renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

renderWindow.Render()
renderWindowInteractor.Start()

reader.SetDefaultExecutivePrototype(None)
print "OK! \n"
#-----------------------------------------------------------------------------#

