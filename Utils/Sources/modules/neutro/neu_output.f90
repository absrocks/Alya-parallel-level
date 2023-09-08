!-----------------------------------------------------------------------
!> @addtogroup Neutro
!> @{
!> @file    neu_output.f90
!> @date    01/04/2016
!> @author  Guillaume Houzeaux
!> @brief   Finalize coupling iteration
!> @details Check convergence 
!> @} 
!-----------------------------------------------------------------------

subroutine neu_output()

  use def_master
  use def_domain
  use def_neutro
  use mod_matrix, only : matrix_output_gid_format
  implicit none
  integer(ip) :: ivari,ivarp
  !
  ! Initial solution, end of a time step and end of run
  !
  do ivarp = 1,nvarp
     ivari = ivarp
     call posdef(11_ip,ivari)
     call neu_outvar(ivari)
  end do

  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! Calculations on sets
     !
     !call neu_outset()
     !
     ! Calculations on witness points
     !
     !call neu_outwit()

  else if( ittyp == ITASK_ENDRUN ) then
     !
     ! End of the run
     !
     !call pspltm(&
     !     npoin,npoin,1_ip,0_ip,c_dom,r_dom,amatr,&
     !     trim(title)//': '//namod(modul),0_ip,18.0_rp,'cm',&
     !     0_ip,0_ip,2_ip,99_ip)
     !call matrix_output_gid_format(npoin,1_ip,r_dom,c_dom,amatr)

  end if


end subroutine neu_output
