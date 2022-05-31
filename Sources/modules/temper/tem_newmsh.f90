subroutine tem_newmsh
!-----------------------------------------------------------------------
!****f* Temper/tem_newmsh
! NAME 
!    tem_newmsh
! DESCRIPTION
!    This routine projects data from the old to the new mesh
! USES
!    tem_openfi
!    tem_reabcs
!    tem_memall
! USED BY
!    Temper
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
     ! call tem_openfi
     ! Read the boundary conditions
     ! call tem_reabcs
     ! Allocate memory (veloc, press and solver)
     ! call tem_memall

  end if
  
end subroutine tem_newmsh
