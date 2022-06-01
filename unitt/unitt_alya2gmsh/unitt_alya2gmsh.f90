program unitt_alya2gmsh
  !
  ! Test copy, extract and merge operators
  !
  use def_elmtyp
  use def_kintyp_basic
  use def_kintyp_mesh_basic
  use mod_alya2gmsh,                only : alya2gmsh_initialization, alya2gmsh_remeshing
  use mod_elmgeo
  use def_domain,       only : ndime
  
  implicit none

  type(mesh_type_basic) :: mesh
  type(mesh_type_basic) :: mesh_new
  integer(ip)           :: ipoin,idime,ielem,inode
  real(rp)              :: min_size_amr
  real(rp)              :: max_size_amr
  real(rp),    pointer  :: mesh_size(:)
  !
  ! Initialize problem dimension value
  !
#ifdef NDIMEPAR
#else  
  ndime = 2_ip;
#endif  


  
  call elmgeo_element_type_initialization()

  call mesh     % init('MY_MESH')
  call mesh_new % init('MY_MESH')
  
  mesh % nelem = 9
  mesh % npoin = 16
  mesh % ndime = 2
  mesh % mnode = 4
  call mesh % alloca()

  mesh % lnods(:,1)  = (/  12 ,7 ,5 ,10  /)
  mesh % lnods(:,2)  = (/  14 ,9 ,7 ,12   /)
  mesh % lnods(:,3)  = (/  16, 15, 9 ,14  /) 
  mesh % lnods(:,4)  = (/  7, 4 ,2 ,5   /)
  mesh % lnods(:,5)  = (/  9 ,8 ,4 ,7  /) 
  mesh % lnods(:,6)  = (/  15, 13, 8, 9 /)  
  mesh % lnods(:,7)  = (/  4 ,3 ,1 ,2  /) 
  mesh % lnods(:,8)  = (/  8 ,6 ,3 ,4  /) 
  mesh % lnods(:,9)  = (/  13 ,11, 6, 8  /)

  mesh % ltype(1:9)  = QUA04

  mesh % coord(:,1 ) = (/   0.000000e+00_rp ,  1.000000e+00_rp /) 
  mesh % coord(:,2 ) = (/   0.000000e+00_rp ,  6.666667e-01_rp /) 
  mesh % coord(:,3 ) = (/   3.333333e-01_rp ,  1.000000e+00_rp /) 
  mesh % coord(:,4 ) = (/   3.333333e-01_rp ,  6.666667e-01_rp /) 
  mesh % coord(:,5 ) = (/   0.000000e+00_rp ,  3.333333e-01_rp /) 
  mesh % coord(:,6 ) = (/   6.666667e-01_rp ,  1.000000e+00_rp /) 
  mesh % coord(:,7 ) = (/   3.333333e-01_rp ,  3.333333e-01_rp /) 
  mesh % coord(:,8 ) = (/   6.666667e-01_rp ,  6.666667e-01_rp /) 
  mesh % coord(:,9 ) = (/   6.666667e-01_rp ,  3.333333e-01_rp /) 
  mesh % coord(:,10) = (/   0.000000e+00_rp ,  0.000000e+00_rp /) 
  mesh % coord(:,11) = (/   1.000000e+00_rp ,  1.000000e+00_rp /) 
  mesh % coord(:,12) = (/   3.333333e-01_rp ,  0.000000e+00_rp /) 
  mesh % coord(:,13) = (/   1.000000e+00_rp ,  6.666667e-01_rp /) 
  mesh % coord(:,14) = (/   6.666667e-01_rp ,  0.000000e+00_rp /) 
  mesh % coord(:,15) = (/   1.000000e+00_rp ,  3.333333e-01_rp /) 
  mesh % coord(:,16) = (/   1.000000e+00_rp ,  0.000000e+00_rp /) 
  !
  ! Initialize alya2gmsh permutation arrays between gmsh and Alya element type identifiers
  !
  call alya2gmsh_initialization()

  allocate(mesh_size(1))
  mesh_size(1) = 0.0_rp
  min_size_amr = 0.0_rp
  max_size_amr = 0.0_rp 
  call alya2gmsh_remeshing(mesh_new,mesh,mesh_size,min_size_amr,max_size_amr)

  do ielem = 1,mesh % nelem
     do inode = 1,mesh % mnode
        if( mesh_new % lnods(inode,ielem) /= mesh_new % lnods(inode,ielem) ) then
           print*,'Wrong lnods'
           stop 1
        end if
     end do
  end do    
  do ipoin = 1,mesh % npoin
     do idime = 1,mesh % ndime
        if( abs(mesh_new % coord(idime,ipoin)-mesh % coord(idime,ipoin)) > 1.0e-12_rp ) then           
           print*,'Wrong coord'
           stop 1
        end if
     end do
  end do

  
end program unitt_alya2gmsh
