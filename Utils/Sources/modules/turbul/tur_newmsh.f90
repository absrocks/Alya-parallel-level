subroutine tur_newmsh
!-----------------------------------------------------------------------
!****f* Turbul/tur_newmsh
! NAME 
!    tur_newmsh
! DESCRIPTION
!    This routine projects data from the old to the new mesh
! USES
!    tur_reabcs
!    tur_memall
! USED BY
!    Turbul
!***
!-----------------------------------------------------------------------
  use      def_master
  implicit none

  if(kfl_algor_msh==1) then                  ! SHEAR_SLIP


  elseif(kfl_algor_msh==2) then              ! GLOBAL REFINEMENT

     ! Save data and deallocate vectors

     ! Project data onto new mesh

     ! Read new mesh data (bc's)


     ! Read the boundary conditions
     ! call tur_reabcs
     ! Allocate memory (veloc, press and solver)
     ! call tur_memall

  end if
  
end subroutine tur_newmsh
