import vtk
 
def main():
    rw = vtk.vtkRenderWindow()
    rw.SetSize(1000, 600)
    rw.SetWindowName("Alya")

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(rw)

    #Create a sphere
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
    
    color_background = []
    color_background.append( (0.8,0.8,0.9) ) 
    color_background.append( (0.8,0.8,0.8) ) 
    color_background.append( (0.8,0.8,0.7) ) 

    xmins=[0.0, 0.5, 0.5]
    ymins=[0.0, 0.5, 0.0]
    xmaxs=[0.5, 1.0, 1.0]
    ymaxs=[1.0, 1.0, 0.5]

    iren_list = []    
    for i in range(3):
        ren = vtk.vtkRenderer()
        rw.AddRenderer(ren)
        ren.SetViewport(xmins[i],ymins[i],xmaxs[i],ymaxs[i])
 
        mapper = vtk.vtkPolyDataMapper()
        actor  = vtk.vtkActor()
        
        mapper.SetInputConnection(sphereSource[i].GetOutputPort())
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(color[i])
        actor.GetProperty().SetRepresentationToWireframe()
        
        ren.AddActor(actor)
        ren.SetBackground(color_background[i])
        ren.ResetCamera()
 
    rw.Render()
    iren.Start()
 
 
if __name__ == '__main__':
    main()
    
