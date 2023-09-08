import vtk_multiblock 
import vtk_tools as tools 

Fin  = "star.case"

ENSI = vtk_multiblock.Vtk_ensight_reader()

ENSI.init(Fin) 

ENSI.set_composite() 
ENSI.BlockIdScalars()

ENSI.save_multiblock() 
#ENSI.save_vtp() 
ENSI.save_vtu() 


#VTU_COSA = vtk_multiblock.__vtp2vtu__(ENSI.AppendFilter)

COSA_IN = ENSI.AppendFilter.GetOutput()

VTP_THRESHOLD = tools.__threshold__(COSA_IN , 1, 1, "CurrentFlatIndex")
VTU_THRESHOLD = tools.__vtp2vtu__(VTP_THRESHOLD)
tools.__save_data_coords__(VTU_THRESHOLD , "zzz") 
tools.__save_data_cell__(VTU_THRESHOLD, "zzz") 


#VTP_SURFACE = tools.__surface_filter__( ENSI.AppendFilter )
VTP_SURFACE = tools.__surface_filter__( VTP_THRESHOLD )
VTU_SURFACE = tools.__vtp2vtu__(VTP_SURFACE)

GlobalIdx   = tools.__point_locator__( VTU_THRESHOLD, VTU_SURFACE )
tools.__save_data_coords__(VTU_SURFACE, "yyy") 
tools.__save_data_cell__(VTU_SURFACE, "yyy", GlobalIdx) 

PROP = tools.__get_cell_data__(VTU_SURFACE, "CurrentFlatIndex")
tools.__save_raw__(PROP, "yyy_ON_BOUNDARIES.alya", "ON_BOUNDARIES", "END_ON_BOUNDARIES")

