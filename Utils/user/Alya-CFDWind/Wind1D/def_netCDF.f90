module def_netcdf

  
  integer(4)                     :: ncID
  integer(4)                     :: nz_ID,nt_ID,site_ID,z_ID,time_ID,lat_ID,lon_ID
  integer(4)                     :: u_ID,v_ID,vel_ID,dir_ID,th_ID
  integer(4)                     :: uadv_ID,vadv_ID,thadv_ID
  integer(4)                     :: tke_ID,eps_ID,mut_ID,lm_ID,Rf_ID
  integer(4)                     :: prdm_ID,prdt_ID,hfx_ID,strs_ID
  integer(4)                     :: ust_ID,sflux_ID,LMO_ID,T2_ID
  integer(4)                     :: istat
  integer(4)                     :: nt_netCDF = 0
  integer(4)                     :: iout = 0       !output number 

  character(len=256)           :: fnc_out,rm
  logical                         :: there
  !Attributes
  character(len=25 ) :: attr_units                
  character(len=100) :: attr_desc,attr_coor
  !netCDF file varialbles names
  character(len=256)           :: time_name = 'time' !time in days since 01/01/0000
  character(len=256)           :: z_name = 'z'       !height levels
  character(len=256)           :: site_name = 'site' !location
  character(len=256)           :: fc_name = 'fc'     !Coriolis term
  character(len=256)           :: ug_name = 'Ug'     !Ug = Press. grad. RHS term (divided by fc) V-component
  character(len=256)           :: vg_name = 'Vg'     !Vg = - Press. grad. RHS term (divided by fc) U-component
  character(len=256)           :: uadv_name = 'Uadv' !Advective momentum RHS term (divided by fc) U-component
  character(len=256)           :: vadv_name = 'Vadv' !Advective momentum RHS term (divided by fc) V-component
  character(len=256)           :: thadv_name = 'Thadv' !Potential temperature advective RHS term
  character(len=256)           :: T2_name = 'T2'     !2m temperature
  character(len=256)           :: tsk_name = 'TSK'   !Skin temperature (Surface temperature)
  character(len=256)           :: u_name = "U"       !U velocity component
  character(len=256)           :: v_name = "V"       !V velocity component
  character(len=256)           :: vel_name = "VEL"   !Velocity module
  character(len=256)           :: dir_name = "DIR"   !Meteorological direction  
  character(len=256)           :: th_name = "Th"     !Potential temperature
  character(len=256)           :: tke_name = "TKE"   !Turbulent kinetic energy
  character(len=256)           :: eps_name = "eps"   !Diss4ation
  character(len=256)           :: mut_name = "mut"   !Turbulent viscosity
  character(len=256)           :: lm_name = "lm"     !Mixing lenght
  character(len=256)           :: Rf_name = "Rf"     !Flux Richardson Number
  character(len=256)           :: prdm_name = "prdm" !TKE Mechanical Production
  character(len=256)           :: prdt_name = "prdt" !TKE Thermal Production
  character(len=256)           :: hfx_name = "hfx"   !Heat flux
  character(len=256)           :: strs_name = "strs" !Shear stress
  !NEW 2018 for Meso-Micro Challange
  character(len=256)           :: ust_name = "us"     !Friction Velocity
  character(len=256)           :: sflux_name = "wt"  !Kinematic Heat Flux
  character(len=256)           :: LMO_name = "L"     !Obukhov Lenght

  
  real(8),allocatable            :: time_out      
  real(8),allocatable            :: u_out    (:)   !nz
  real(8),allocatable            :: v_out    (:)   !nz
  real(8),allocatable            :: vel_out  (:)   !nz
  real(8),allocatable            :: dir_out  (:)   !nz
  real(8),allocatable            :: th_out   (:)   !nz
  real(8),allocatable            :: uadv_out (:)   !nz
  real(8),allocatable            :: vadv_out (:)   !nz
  real(8),allocatable            :: thadv_out(:)   !nz
  real(8),allocatable            :: tke_out  (:)   !nz
  real(8),allocatable            :: eps_out  (:)   !nz
  real(8),allocatable            :: mut_out  (:)   !nz
  real(8),allocatable            :: lm_out   (:)   !nz
  real(8),allocatable            :: Rf_out   (:)   !nz
  real(8),allocatable            :: prdm_out (:)   !nz
  real(8),allocatable            :: prdt_out (:)   !nz
  real(8),allocatable            :: hfx_out  (:)   !nz
  real(8),allocatable            :: strs_out (:)   !nz
  real(8),allocatable            :: ust_out        
  real(8),allocatable            :: sflux_out      
  real(8),allocatable            :: LMO_out        
  real(8),allocatable            :: T2_out         

end module def_netcdf
