import vtk_tools as tools 


VTP = tools.__vtp_reader__("YMI/surface01.vtp") 

VTU = []
VTU.append( tools.__vtu_reader__("YMI/volume01.vtu")  )
VTU.append( tools.__vtp2vtu__(VTP) )

GlobalIdx = tools.__point_locator__( VTU[0].GetOutput(), VTU[1] )

tools.__save_data_coords__(VTU[0].GetOutput(), "zzz") 
tools.__save_data_cell__(VTU[0].GetOutput(), "zzz") 

tools.__save_data_coords__(VTU[1], "yyy") 
tools.__save_data_cell__(VTU[1], "yyy", GlobalIdx) 

PROP = tools.__get_cell_data__(VTU[1], "CurrentFlatIndex")
tools.__save_raw__(PROP, "yyy_ON_BOUNDARIES.alya", "ON_BOUNDARIES", "END_ON_BOUNDARIES")

