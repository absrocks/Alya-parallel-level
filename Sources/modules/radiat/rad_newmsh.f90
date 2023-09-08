subroutine rad_newmsh
!-----------------------------------------------------------------------
!****f* Radiat/rad_newmsh
! NAME 
!    rad_newmsh
! DESCRIPTION
!    This routine projects data from the old to the new mesh
! USES
!    rad_openfi
!    rad_reabcs
!    rad_memall
! USED BY
!    Radiat
!***
!-----------------------------------------------------------------------
  use      def_master
  implicit none

  if(kfl_algor_msh==1) then                  ! SHEAR_SLIP


  elseif(kfl_algor_msh==2) then              ! GLOBAL REFINEMENT

     ! Save data and deallocate vectors

     ! Project data onto new mesh

     ! Read new mesh data (bc's)

     ! Open files
     ! call rad_openfi
     ! Read the boundary conditions
     ! call rad_reabcs
     ! Allocate memory (veloc, press and solver)
     ! call rad_memall

  end if
  
end subroutine rad_newmsh
