import vtk_creator as Creator


fin  = "banc_interpolated.vtu"
fout = "banc_alya_format"


VTU = Creator.Vtk_unstructured_reader() 
VTU.load_vtk_structured(fin)
V1 = VTU.get_points_prop_scalar("Velocityi")
V2 = VTU.get_points_prop_scalar("Velocityj")
V3 = VTU.get_points_prop_scalar("Velocityk")


falya = open("VELOC.alya", "w")

i = 0
for v1, v2 in zip(V1, V2):
  print>> falya, i+1, v1, v2  
  i += 1
falya.close() 

#VTU.save_data_coords(fout)
#VTU.save_data_cell(fout) 


print "OK!!"
print 
