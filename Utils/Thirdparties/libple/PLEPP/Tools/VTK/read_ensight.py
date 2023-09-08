import vtk 
import sys 
"""
if XXX.SetInput[Data](       YYY.GetOutput() )
if XXX.SetInputConnection( YYY.GetOutputPort() )
"""


PROPERTY = "ZVELO"

cdp    = vtk.vtkCompositeDataPipeline()
reader = vtk.vtkGenericEnSightReader() 
reader.SetDefaultExecutivePrototype(cdp)
reader.SetCaseFileName('runner.ensi.case')
#reader.SetTimeValue(0.05)
reader.ReadAllVariablesOff() # ReadAllVariablesOn() 
reader.Update()


print "Times:",
print reader.GetTimeSets().GetItem(0).GetRange()

Times = [] 
for i in range(reader.GetTimeSets().GetItem(0).GetSize()): 
  Times.append( reader.GetTimeSets().GetItem(0).GetValue(i) ) 

time = Times[-1]
print time 

reader.SetTimeValue(time)
#reader.ReadAllVariablesOff() # ReadAllVariablesOn() 
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


#------------------------------------------------------------| plot table |----# 
# line 
probeLine = vtk.vtkLineSource()
probeLine.SetPoint1(1.0,0.0,-2.0)
probeLine.SetPoint2(1.0,0.0,20.0)
probeLine.SetResolution(100)
probe = vtk.vtkProbeFilter()
probe.SetInputConnection(probeLine.GetOutputPort())
probe.SetSource(DATAi)
probe.Update()
#
probeTube = vtk.vtkTubeFilter()
probeTube.SetInput(probe.GetPolyDataOutput())
probeTube.SetNumberOfSides(5)
probeTube.SetRadius(.05)
probeMapper = vtk.vtkPolyDataMapper()
probeMapper.SetInputConnection(probeTube.GetOutputPort())
probeMapper.SetScalarRange(set_range)
probeActor = vtk.vtkActor()
probeActor.SetMapper(probeMapper)


from vtk_plot_table import Vtk_plot_table
X = []
Y = [] 

numPoints = probe.GetOutput().GetNumberOfPoints() 
for i in range(numPoints): X.append(i)
for i in range(numPoints): Y.append(probe.GetOutput().GetPointData().GetArray(PROPERTY).GetValue(i))

print "LEN:", 
print len(X), 
print len(Y)

P = Vtk_plot_table()
P.set_pts_size(numPoints) 
P.set_prop("Nodes")
P.set_prop(PROPERTY)
P.set_vals(X, "Nodes")
P.set_vals(Y, PROPERTY)

#P.plot() 


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
#
## scalarbar 
scalarBar = vtk.vtkScalarBarActor()
scalarBar.SetLookupTable(reader_mapper.GetLookupTable()) 
scalarBar.SetTitle(PROPERTY)
scalarBar.SetNumberOfLabels(5)
scalarBar.SetOrientationToHorizontal()
scalarBar.SetWidth(0.50)
scalarBar.SetHeight(0.1)
scalarBar.GetPositionCoordinate().SetValue(0.05,0.01)


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


## 


#--------------------------------------------------------------| rendering |----# 
renderer = vtk.vtkRenderer()
renderer.SetBackground((0.7,0.7,0.7))
renderer.AddActor(reader_actor)
renderer.AddActor(contour_actor)
renderer.AddActor(dss_actor)
renderer.AddActor(cut_actor)
renderer.AddActor(isoActor)
renderer.AddActor(probeActor)
renderer.AddActor2D(scalarBar)


#renderer.AddActor(outline_actor) # segmentation fault 


renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindow.StereoCapableWindowOn()
renderWindow.SetSize(600, 800)

renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

#renderWindow.Render()
#renderWindowInteractor.Start()

#reader.SetDefaultExecutivePrototype(None)


from vtk_plot_linesource import Vtk_plot_linesource

numPoints = 250 
H = 20.0/1000

vtkL = [] 
vtkL.append( Vtk_plot_linesource() ) 
vtkL.append( Vtk_plot_linesource() ) 
vtkL.append( Vtk_plot_linesource() ) 
vtkL.append( Vtk_plot_linesource() ) 

vtkL[0].set_source([0.01,0.0,0.0], [0.01,0.05,0.0], numPoints, set_range) 
vtkL[0].set_connection(DATAi)
vtkL[0].save_data("xxx01.dat", PROPERTY) 

vtkL[1].set_source([0.01,0.0,2*H], [0.01,0.05,2*H], numPoints, set_range) 
vtkL[1].set_connection(DATAi)
vtkL[1].save_data("xxx02.dat", PROPERTY) 

vtkL[2].set_source([0.01,0.0,4*H], [0.01,0.05,4*H], numPoints, set_range) 
vtkL[2].set_connection(DATAi)
vtkL[2].save_data("xxx03.dat", PROPERTY) 

vtkL[3].set_source([0.01,0.025,-0.02], [0.01,0.025,0.2], numPoints, set_range) 
vtkL[3].set_connection(DATAi)
vtkL[3].save_data("xxx04.dat", PROPERTY) 



X = range(numPoints)
Y = vtkL[0].get_data(PROPERTY) 

P = Vtk_plot_table()
P.set_pts_size(numPoints) 
P.set_prop(PROPERTY)
P.set_prop("Nodes")

P.set_vals(X, "Nodes")
P.set_vals(Y, PROPERTY)
#P.plot() 


from vtk_plot_cutter import Vtk_plot_cutter

vtkXY01 = Vtk_plot_cutter()
vtkXY01.set_plane([0.01,0.025,0.0], [0,0,1], set_range) 
vtkXY01.set_connection(reader) 

vtkXY02 = Vtk_plot_cutter()
vtkXY02.set_plane([0.01,0.025,0.01], [0,0,1], set_range) 
vtkXY02.set_connection(reader) 

vtkYZ01 = Vtk_plot_cutter()
vtkYZ01.set_plane([0.01,0.025,0.0], [0,1,0], set_range) 
vtkYZ01.set_connection(reader) 


from vtk_plot_contour import Vtk_plot_contour

vtkC = Vtk_plot_contour()
vtkC.set_value(set_average*0.5) 
vtkC.set_value(set_average*1.0) 
vtkC.set_value(set_average*1.5) 
vtkC.set_connection(reader)


from vtk_plot_multiview import Vtk_plot_multiview

subactors = []
subactors.append( reader_actor ) 
subactors.append( vtkC.get_actor() ) 
subactors.append( vtkYZ01.get_actor() ) 
subactors.append( vtkL[0].get_actor() ) 
subactors.append( vtkL[1].get_actor() ) 
subactors.append( vtkL[2].get_actor() ) 
#subactors.append( probeActor ) 


actors = [] 
actors.append( subactors ) 
actors.append( [vtkXY01.get_actor()] ) 
actors.append( [vtkXY02.get_actor()] ) 

MV = Vtk_plot_multiview()
MV.set_actor2D(scalarBar)
MV.set_actors(actors)
MV.render()


print "OK! \n"
#-----------------------------------------------------------------------------#

