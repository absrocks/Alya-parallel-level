!---------------------------------------------------------------------!
!This subroutine writes a VTK file for Paraview. 
!Written by S. Gopalakrishnan on 4/10 
!and James F. Kelly, NPS, Monterey (CA)
! 
!  To write one file every n-time steps use this
!  call:
!  
!        !
!        !Write vtk every n-steps
!        !
!        if(ittim == 0 .or. mod(ittim,irestart) == 0) then
!           write(fnp1,'(i3)')ittim
!           iloop=2 - int(log10(real(ittim)))
!           do j=1,iloop
!              fnp1(j:j)='0'
!           end do
!           fnp=trim('OUTVTK_alya') // '_' // trim(fnp1) // '.vtk'
!           call nsa_outvtk(fnp)
!        end if
!
!
!---------------------------------------------------------------------!
subroutine nsa_outvtk(fname)  
  
  use      def_master
  use      def_domain
  use      mod_postpr
  use      def_nastal

  implicit none
  
  !global
  real      :: theta, rho
  character :: fname*72
  integer   :: ie, i, j, k, nglm1, nglm13
  integer   :: ncells, nsize, icomp
  
  nglm13 = 1                   !Number of cells per element
  ncells = nelem*nglm13        !Total number of cells
  nsize = 5*ncells             !Number of points required to represent cells
  icomp=min(TIME_N,ncomp_nsa)
  
  !Open VTK file
  open(1,file=trim(fname))
  
  write(1,'(a)')"# vtk DataFile Version 2.0"
  write(1,'(a)')'EULER2D data'
  write(1,'(a)')'ASCII'
  write(1,'(a)')'DATASET UNSTRUCTURED_GRID'
  
  !Write out the coordinates
  write(1,*)'POINTS',npoin,'float'
  do i=1,npoin
     write(1,*) coord(1,i),coord(2,i), 0.0
  end do
  
  !Write out the connectivity
  write(1,*)'CELLS',ncells,nsize
  do ie =1,nelem
     write(1, *)'4', &
          (lnods(1,ie)-1), &
          (lnods(2,ie)-1), &
          (lnods(3,ie)-1), &
          (lnods(4,ie)-1)
  end do
  
  write(1,*)'CELL_TYPES',ncells
  do i=1,ncells
     write(1,*) '9'
  end do

  write(1,*)"CELL_DATA",ncells
  write(1,*)"POINT_DATA",npoin
  
  write(1,*)"SCALARS DENSI-HY double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) densi(i,icomp) - rekee_nsa(ndime+1, i)
  end do
  write(1,*)"SCALARS DENSI double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) densi(i,icomp)
  end do
  write(1,*)"SCALARS PRESS-HY double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) press(i,icomp) - rekee_nsa(ndime+3, i)
  end do
  write(1,*)"SCALARS PRESS double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) press(i,icomp)
  end do
  write(1,*)"SCALARS U-velo double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) veloc(1,i,icomp)
  end do
  write(1,*)"SCALARS W-velo double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) veloc(ndime,i,icomp)
  end do
  write(1,*)"SCALARS Theta-HY double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*)  tempe(i,icomp) - rekee_nsa(ndime+2, i)
  end do
  write(1,*)"SCALARS Theta double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) tempe(i,icomp)
  end do
  
  write(1,*)"VECTORS Velocity double"
  do i=1,npoin
     write(1,*) veloc(1,i,icomp), veloc(2,i,icomp), 0.0
  end do

  close(1)
  
end subroutine nsa_outvtk

subroutine nsa_outvtk_kessler(fname,q)
  
  use      def_master
  use      def_domain
  use      mod_postpr
  use      def_nastal
  
  implicit none
  
  !global
  real(rp)    :: q(nvar_nsa,npoin), theta, rho
  character   :: fname*72
  integer(ip) :: ie, i, j, k, icol,nglm1, nglm13
  integer(ip) :: ncells, nsize, icomp
  
  nglm13 = 1                   !Number of cells per element
  ncells = nelem*nglm13        !Total number of cells
  nsize = 5*ncells             !Number of points required to represent cells
  icomp=min(TIME_N,ncomp_nsa)
  
  !Open VTK file
  open(1,file=fname)
  
  write(1,'(a)')"# vtk DataFile Version 2.0"
  write(1,'(a)')'EULER2D data'
  write(1,'(a)')'ASCII'
  write(1,'(a)')'DATASET UNSTRUCTURED_GRID'
  write(1,'(a)')'FIELD FieldData 2'
  write(1,'(a)')'TIME 1 1 double'
  write(1,*) cutim
  write(1,'(a)')'CYCLE 1 1 int'
  write(1,*) ittim
  !Write out the coordinates
  write(1,*)'POINTS',npoin,'float'
  do i=1,npoin
     write(1,*) coord(1,i),coord(2,i), 0.0
  end do
  
  !Write out the connectivity
  write(1,*)'CELLS',ncells,nsize
  do ie =1,nelem
     write(1, *)'4', &
          (lnods(1,ie)-1), &
          (lnods(2,ie)-1), &
          (lnods(3,ie)-1), &
          (lnods(4,ie)-1)
  end do
  
  write(1,*)'CELL_TYPES',ncells
  do i=1,ncells
     write(1,*) '9'
  end do

  write(1,*)"CELL_DATA",ncells
  write(1,*)"POINT_DATA",npoin
  
  write(1,*)"SCALARS DENSI-HY double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) q(1,i)
  end do
  write(1,*)"SCALARS DENSI double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) q(1,i) + rekee_nsa(ndime+1, i)
  end do
  write(1,*)"SCALARS Qv double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) q(5,i)
  end do
  write(1,*)"SCALARS Qc double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) q(6,i)
  end do
  write(1,*)"SCALARS Qr double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) q(7,i)
  end do
  write(1,*)"SCALARS U-velo double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) q(2,i)
  end do
  write(1,*)"SCALARS W-velo double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) q(3,i)
  end do
  write(1,*)"SCALARS Theta-HY double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*)  q(4,i)
  end do
  write(1,*)"SCALARS Theta double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) q(4,i) + rekee_nsa(ndime+2, i)
  end do
  write(1,*)"VECTORS Velocity double"
  do i=1,npoin
     write(1,*) q(2,i), q(3,i), 0.0_rp
  end do
  
  close(1)
 
end subroutine nsa_outvtk_kessler

subroutine nsa_outvtk_moist(fname_moist,q,q_ref)  
  
  use      def_master
  use      def_domain
  use      mod_postpr
  use      def_nastal
  
  implicit none
  
  !global
  real(rp)      :: q(7,npoin),q_ref(7,npoin), theta, rho
  character     :: fname_moist*72
  integer(ip)   :: ie, i, j, k, nglm1, nglm13
  integer(ip)   :: ncells, nsize, icomp
  
  nglm13 = 1                   !Number of cells per element
  ncells = nelem*nglm13        !Total number of cells
  nsize = 5*ncells             !Number of points required to represent cells
  icomp=min(TIME_N,ncomp_nsa)
  
  !Open VTK file
  open(1,file=fname_moist)
  
  write(1,'(a)')"# vtk DataFile Version 2.0"
  write(1,'(a)')'EULER2D data'
  write(1,'(a)')'ASCII'
  write(1,'(a)')'DATASET UNSTRUCTURED_GRID'
  
  !Write out the coordinates
  write(1,*)'POINTS',npoin,'float'
  do i=1,npoin
     write(1,*) coord(1,i),coord(2,i), 0.0
  end do
  
  !Write out the connectivity
  write(1,*)'CELLS',ncells,nsize
  do ie =1,nelem
     write(1, *)'4', &
          (lnods(1,ie) -1), &
          (lnods(2,ie)-1), &
          (lnods(3,ie)-1), &
          (lnods(4,ie)-1)
  end do
  
  write(1,*)'CELL_TYPES',ncells
  do i=1,ncells
     write(1,*) '9'
  end do

  write(1,*)"CELL_DATA",ncells
  write(1,*)"POINT_DATA",npoin
!!$  
!!$  write(1,*)"SCALARS DENSI double"
!!$  write(1,*) "LOOKUP_TABLE default"
!!$  do i=1,npoin
!!$     write(1,*) q(1,i)
!!$  end do
!!$  write(1,*)"SCALARS DENSI_REF double"
!!$  write(1,*) "LOOKUP_TABLE default"
!!$  do i=1,npoin
!!$     write(1,*)  rekee_nsa(3,i)
!!$  end do
!!$  write(1,*)"SCALARS PRESS-HY double"
!!$  write(1,*) "LOOKUP_TABLE default"
!!$  do i=1,npoin
!!$     write(1,*) press(i,icomp) - rekee_nsa(1, i)
!!$  end do
!!$  write(1,*)"SCALARS PRESS double"
!!$  write(1,*) "LOOKUP_TABLE default"
!!$  do i=1,npoin
!!$     write(1,*) press(i,icomp)
!!$  end do
!!$
!!$  write(1,*)"SCALARS U-velo double"
!!$  write(1,*) "LOOKUP_TABLE default"
!!$  do i=1,npoin
!!$     write(1,*) veloc(1,i,icomp)
!!$  end do
!!$  write(1,*)"SCALARS W-velo double"
!!$  write(1,*) "LOOKUP_TABLE default"
!!$  do i=1,npoin
!!$     write(1,*) veloc(2,i,icomp)
!!$  end do
!!$  write(1,*)"SCALARS theta double"
!!$  write(1,*) "LOOKUP_TABLE default"
!!$  do i=1,npoin
!!$     write(1,*)  q(4,i)
!!$  end do
!!$ 
  write(1,*)"SCALARS qv double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) q(5,i)
  end do
  write(1,*)"SCALARS qc double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) q(6,i)
  end do
  write(1,*)"SCALARS qr double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) q(7,i)
  end do
  
  close(1)
  
end subroutine nsa_outvtk_moist

subroutine nsa_outvtk_water(fname)  
  
  use      def_master
  use      def_domain
  use      mod_postpr
  use      def_nastal

  implicit none
  
  !global
  real(rp)      :: theta, rho
  character     :: fname*72
  integer(ip)   :: ie, i, j, k, nglm1, nglm13
  integer(ip)   :: ncells, nsize, icomp
  
  nglm13 = 1                   !Number of cells per element
  ncells = nelem*nglm13        !Total number of cells
  nsize = 5*ncells             !Number of points required to represent cells
  icomp=min(TIME_N,ncomp_nsa)
  
  !Open VTK file
  open(1,file=fname)
  
  write(1,'(a)')"# vtk DataFile Version 2.0"
  write(1,'(a)')'EULER2D data'
  write(1,'(a)')'ASCII'
  write(1,'(a)')'DATASET UNSTRUCTURED_GRID'
  
  !Write out the coordinates
  write(1,*)'POINTS',npoin,'float'
  do i=1,npoin
     write(1,*) coord(1,i),coord(2,i), 0.0
  end do
  
  !Write out the connectivity
  write(1,*)'CELLS',ncells,nsize
  do ie =1,nelem
     write(1, *)'4', &
          (lnods(1,ie)-1), &
          (lnods(2,ie)-1), &
          (lnods(3,ie)-1), &
          (lnods(4,ie)-1)
  end do
  
  write(1,*)'CELL_TYPES',ncells
  do i=1,ncells
     write(1,*) '9'
  end do

  write(1,*)"CELL_DATA",ncells
  write(1,*)"POINT_DATA",npoin
  
  write(1,*)"SCALARS CONCE1 double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) conce(i,1,3)
  end do
  write(1,*)"SCALARS CONCE2 double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) conce(i,2,3)
  end do
  write(1,*)"SCALARS CONCE3 double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) conce(i,3,3)
  end do
  close(1)
  
end subroutine nsa_outvtk_Water

subroutine nsa_outvtk_pert(fname)
  
  use      def_master
  use      def_domain
  use      mod_postpr
  use      def_nastal
  
  implicit none
  
  !global
  real theta, rho, q(7,npoin), q_ref(7,npoin)
  character fname*72
  integer ie, i, j, k, nglm1, nglm13
  integer ncells, nsize, icomp
  
  nglm13 = 1                   !Number of cells per element
  ncells = nelem*nglm13        !Total number of cells
  nsize = 5*ncells             !Number of points required to represent cells
  icomp=min(TIME_N,ncomp_nsa)
  
  !Open VTK file
  open(1,file=fname)
  
  write(1,'(a)')"# vtk DataFile Version 2.0"
  write(1,'(a)')'EULER2D data'
  write(1,'(a)')'ASCII'
  write(1,'(a)')'DATASET UNSTRUCTURED_GRID'
  
  !Write out the coordinates
  write(1,*)'POINTS',npoin,'float'
  do i=1,npoin
     write(1,*) coord(1,i),coord(2,i), 0.0
  end do
  
  !Write out the connectivity
  write(1,*)'CELLS',ncells,nsize
  do ie =1,nelem
     write(1, *)'4', &
          (lnods(1,ie) -1), &
          (lnods(2,ie)-1), &
          (lnods(3,ie)-1), &
          (lnods(4,ie)-1)
  end do
  
  write(1,*)'CELL_TYPES',ncells
  do i=1,ncells
     write(1,*) '9'
  end do

  write(1,*)"CELL_DATA",ncells
  write(1,*)"POINT_DATA",npoin
  
  write(1,*)"SCALARS DPERT double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) densi(i,icomp)
  end do
  write(1,*)"SCALARS DENSI double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) densi(i,icomp) + rekee_nsa(3, i)
  end do
  !
  ! Pressure is initialized as total regardless
  ! of the set of variables used. This is why 
  ! we write differently wrt densi and tempe.
  !
  write(1,*)"SCALARS PPERT double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) press(i,icomp) - rekee_nsa(1, i)
  end do
  write(1,*)"SCALARS PRESS double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) press(i,icomp)
  end do

  write(1,*)"SCALARS U-velo double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) veloc(1,i,icomp)
  end do
  write(1,*)"SCALARS W-velo double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) veloc(2,i,icomp)
  end do
  write(1,*)"SCALARS TPERT double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*)  tempe(i,icomp)
  end do
  write(1,*)"SCALARS Theta double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) tempe(i,icomp) + rekee_nsa(2, i)
  end do
  write(1,*)"VECTORS Velocity double"
  do i=1,npoin
     write(1,*) veloc(1,i,icomp), veloc(2,i,icomp), 0.0
  end do

  close(1)
  
end subroutine nsa_outvtk_pert


subroutine nsa_outmatlab(fname,q)
  
  use      def_master
  use      def_domain
  
  implicit none
  
  !global
  real(rp)    :: q(npoin)
  character   :: fname*72
  integer(ip) :: i
      
  !Open dat file
  open(1,file=fname)
    do i=1,npoin
     write(1,*) q(i)
  end do
  
  close(1)
 
end subroutine nsa_outmatlab


subroutine nsa_outvtk_3d(fname)  
  
  use      def_master
  use      def_domain
  use      mod_postpr
  use      def_nastal

  implicit none
  
  !global
  real      :: theta, rho
  character :: fname*72
  integer   :: ie, i, j, k, nglm1, nglm13
  integer   :: ncells, nsize, icomp
  
  nglm13 = 1                   !Number of cells per element
  ncells = nelem*nglm13        !Total number of cells
  nsize = 9*ncells             !Number of points required to represent cells
  icomp=min(TIME_N,ncomp_nsa)
  
  !Open VTK file
  open(1,file=fname)
  
  write(1,'(a)')"# vtk DataFile Version 2.0"
  write(1,'(a)')'EULER2D data'
  write(1,'(a)')'ASCII'
  write(1,'(a)')'DATASET UNSTRUCTURED_GRID'
  
  !Write out the coordinates
  write(1,*)'POINTS',npoin,'float'
  do i=1,npoin
     write(1,*) coord(1,i),coord(2,i), coord(3,i)
  end do
  
  !Write out the connectivity
  write(1,*)'CELLS',ncells,nsize
  do ie =1,nelem
     write(1, *)'8', &
          (lnods(1,ie)-1), &
          (lnods(2,ie)-1), &
          (lnods(3,ie)-1), &
          (lnods(4,ie)-1), &
          (lnods(5,ie)-1), &
          (lnods(6,ie)-1), &
          (lnods(7,ie)-1), &
          (lnods(8,ie)-1)
  end do
  
  write(1,*)'CELL_TYPES',ncells
  do i=1,ncells
     write(1,*) '12'
  end do

  write(1,*)"CELL_DATA",ncells
  write(1,*)"POINT_DATA",npoin
  
  write(1,*)"SCALARS DENSI-HY double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) densi(i,icomp) - rekee_nsa(ndime+1, i)
  end do
  write(1,*)"SCALARS DENSI double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) densi(i,icomp)
  end do
  write(1,*)"SCALARS PRESS-HY double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) press(i,icomp) - rekee_nsa(ndime+3, i)
  end do
  write(1,*)"SCALARS PRESS double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) press(i,icomp)
  end do
  write(1,*)"SCALARS U-velo double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) veloc(1,i,icomp)
  end do
  write(1,*)"SCALARS V-velo double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) veloc(2,i,icomp)
  end do
  write(1,*)"SCALARS W-velo double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) veloc(ndime,i,icomp)
  end do
  write(1,*)"SCALARS Theta-HY double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*)  tempe(i,icomp) - rekee_nsa(ndime+2, i)
  end do
  write(1,*)"SCALARS Theta double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) tempe(i,icomp)
  end do
  
  write(1,*)"VECTORS Velocity double"
  do i=1,npoin
     write(1,*) veloc(1,i,icomp), veloc(2,i,icomp), veloc(ndime,i,icomp)
  end do

  close(1)
  
end subroutine nsa_outvtk_3d

subroutine nsa_outvtk_3d_q(fname,q)  
  
  use      def_master
  use      def_domain
  use      mod_postpr
  use      def_nastal

  implicit none
  
  !global
  real      :: theta, rho, q(ndofn_nsa+3,npoin)
  character :: fname*72
  integer   :: ie, i, j, k, nglm1, nglm13
  integer   :: ncells, nsize, icomp
  
  nglm13 = 1                   !Number of cells per element
  ncells = nelem*nglm13        !Total number of cells
  nsize = 9*ncells             !Number of points required to represent cells
  icomp=min(TIME_N,ncomp_nsa)
  
  !Open VTK file
  open(1,file=fname)
  
  write(1,'(a)')"# vtk DataFile Version 2.0"
  write(1,'(a)')'EULER2D data'
  write(1,'(a)')'ASCII'
  write(1,'(a)')'DATASET UNSTRUCTURED_GRID'
  
  !Write out the coordinates
  write(1,*)'POINTS',npoin,'float'
  do i=1,npoin
     write(1,*) coord(1,i),coord(2,i), coord(3,i)
  end do
  
  !Write out the connectivity
  write(1,*)'CELLS',ncells,nsize
  do ie =1,nelem
     write(1, *)'8', &
          (lnods(1,ie)-1), &
          (lnods(2,ie)-1), &
          (lnods(3,ie)-1), &
          (lnods(4,ie)-1), &
          (lnods(5,ie)-1), &
          (lnods(6,ie)-1), &
          (lnods(7,ie)-1), &
          (lnods(8,ie)-1)
  end do
  
  write(1,*)'CELL_TYPES',ncells
  do i=1,ncells
     write(1,*) '12'
  end do

  write(1,*)"CELL_DATA",ncells
  write(1,*)"POINT_DATA",npoin
  
  write(1,*)"SCALARS DENSI-HY double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) densi(i,icomp) - rekee_nsa(ndime+1, i)
  end do
  write(1,*)"SCALARS DENSI double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) densi(i,icomp)
  end do
  write(1,*)"SCALARS PRESS-HY double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) press(i,icomp) - rekee_nsa(ndime+3, i)
  end do
  write(1,*)"SCALARS PRESS double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) press(i,icomp)
  end do
  write(1,*)"SCALARS U-velo double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) veloc(1,i,icomp)
  end do
  write(1,*)"SCALARS V-velo double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) veloc(2,i,icomp)
  end do
  write(1,*)"SCALARS W-velo double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) veloc(ndime,i,icomp)
  end do
  write(1,*)"SCALARS Theta-HY double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*)  tempe(i,icomp) - rekee_nsa(ndime+2, i)
  end do
  write(1,*)"SCALARS Theta double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) tempe(i,icomp)
  end do
  
  write(1,*)"VECTORS Velocity double"
  do i=1,npoin
     write(1,*) veloc(1,i,icomp), veloc(2,i,icomp), veloc(ndime,i,icomp)
  end do

  close(1)
  
end subroutine nsa_outvtk_3d_q
