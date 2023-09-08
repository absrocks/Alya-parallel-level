subroutine chm_outvtk_moist(fname_moist,q,q_ref)  
  
  use      def_master
  use      def_domain
  use      mod_postpr
!  use      def_nastal
  use      def_chemic
  
  implicit none
  
  !global
  real(rp)      :: q(7,npoin),q_ref(7,npoin), theta, rho
  character     :: fname_moist*72
  integer(ip)   :: ie, i, j, k, nglm1, nglm13
  integer(ip)   :: ncells, nsize, icomp
  
  nglm13 = 1                   !Number of cells per element
  ncells = nelem*nglm13        !Total number of cells
  nsize = 5*ncells             !Number of points required to represent cells
    
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
  
end subroutine chm_outvtk_moist


subroutine chm_outvtk(fname_moist)
  
  use      def_master
  use      def_domain
  use      mod_postpr
  use      def_chemic
  
  implicit none
  
  !global
  real(rp)      :: theta, rho
  character     :: fname_moist*72
  integer(ip)   :: ie, i, j, k, nglm1, nglm13
  integer(ip)   :: ncells, nsize, icomp
  
  nglm13 = 1                   !Number of cells per element
  ncells = nelem*nglm13        !Total number of cells
  nsize = 5*ncells             !Number of points required to represent cells
    
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
  write(1,*)"SCALARS qv double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) conce(i,1,3)
  end do
  write(1,*)"SCALARS qc double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) conce(i,2,3)
  end do
  write(1,*)"SCALARS qr double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) conce(i,3,3)
  end do
  
  close(1)
  
end subroutine chm_outvtk

subroutine chm_outvtk_3D(fname_moist)
  
  use      def_master
  use      def_domain
  use      mod_postpr
  use      def_chemic
  
  implicit none
  
  !global
  real(rp)      :: theta, rho
  character     :: fname_moist*72
  integer(ip)   :: ie, i, j, k, nglm1, nglm13
  integer(ip)   :: ncells, nsize, icomp
  
  nglm13 = 1                   !Number of cells per element
  ncells = nelem*nglm13        !Total number of cells
  nsize = 9*ncells             !Number of points required to represent cells
    
  !Open VTK file
  open(1,file=fname_moist)
  
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
  write(1,*)"SCALARS qv double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) conce(i,1,3)
  end do
  write(1,*)"SCALARS qc double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) conce(i,2,3)
  end do
  write(1,*)"SCALARS qr double"
  write(1,*) "LOOKUP_TABLE default"
  do i=1,npoin
     write(1,*) conce(i,3,3)
  end do
  
  close(1)

end subroutine chm_outvtk_3D
