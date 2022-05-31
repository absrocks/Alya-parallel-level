!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_begrun.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Beginning the run... 
!> @details Beginning the run... we can compute matrices!
!> @}
!-----------------------------------------------------------------------

subroutine nsi_begrun()

  use def_master
  use def_nastin
  use def_domain
  use mod_nsi_fractional_step,  only : nsi_fractional_step_matrices
  use mod_nsi_multi_step_fs,    only : nsi_multi_step_fs_matrices
  use mod_nsi_semi_implicit,    only : nsi_semi_implicit_matrices
  use mod_nsi_schur_complement, only : nsi_schur_complement_matrices
  implicit none
  
  !-------------------------------------------------------------------
  !
  ! Solution strategies
  !
  !-------------------------------------------------------------------
  
  if( NSI_FRACTIONAL_STEP ) then
     !
     ! Graction astep
     !
     if(kfl_tisch_nsi == 4) then
        call nsi_multi_step_fs_matrices()
     else
        call nsi_fractional_step_matrices()       
     end if
  else if(NSI_SEMI_IMPLICIT ) then
     call nsi_semi_implicit_matrices()       
     
  else if( NSI_SCHUR_COMPLEMENT ) then
     !
     ! Schur complement
     !
     call nsi_schur_complement_matrices()       
  end if
  !
  ! If pressure matrix comes from a Schur complement
  ! Must be here because amatr must be allocated
  !
  if( NSI_SCHUR_COMPLEMENT ) then
     solve(2) % A1       => amatr(poapp_nsi:)
     solve(2) % A2       => amatr(poapu_nsi:)
     solve(2) % A3       => amatr(poauu_nsi:)
     solve(2) % A4       => amatr(poaup_nsi:)
     solve(2) % ndofn_A3 =  ndime
     nullify(solve(2) % invA3)
  end if
  
end subroutine nsi_begrun
