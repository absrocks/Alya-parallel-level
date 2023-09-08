subroutine nsi_newmsh
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_newmsh
  ! NAME 
  !    nsi_newmsh
  ! DESCRIPTION
  !    This routine set NST data on the new mesh
  ! USES
  !    nsi_updmsh
  !    nsi_promsh
  !    nsi_openfi
  !    nsi_reabcs
  !    nsi_memall
  ! USED BY
  !    Nastin
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_solver 
  implicit none

  if(kfl_algor_msh==1) then                  ! SHEAR_SLIP
     ! Save data and deallocate vectors
     call nsi_updmsh(one)
     !  
     ! Allocate memory for the solver
     !
     !call nsi_memsol
     !
     ! Allocate memory for all the unknowns of the problem and the
     ! coefficients.
     !
     !call memunk


  elseif(kfl_algor_msh==2) then              ! GLOBAL REFINEMENT

     ! Save data and deallocate vectors
     call nsi_updmsh(one)

     ! Project data onto new mesh
     call nsi_promsh

     ! Deallocate old mesh data
     call nsi_updmsh(two)

     ! Read new mesh data (bc's)
     ! Open data file
     !call nsi_openda
     ! Read the boundary conditions
     !call nsi_reabcs
     ! Allocate memory (veloc, press and solver)
     !call nsi_memall

  end if

end subroutine nsi_newmsh
