subroutine hlm_reanut()

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_reanut.f90
  ! NAME 
  !    hlm_reanut
  ! DESCRIPTION
  !    This routine reads the numerical problem for 'helmoz' module
  ! USES
  !    ecoute
  ! USED BY
  !    hlm_turnon
  !-----------------------------------------------------------------------

  use def_parame
  use def_inpout
  use def_master
  use def_solver
  use def_helmoz
  use def_domain
  use mod_ecoute, only :  ecoute

  implicit none

  cotol_hlm     = 1.0e-12_rp  ! Convergence tolerance

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     solve_sol => solve  
     !
     ! Reach the section
     !
     call ecoute('hlm_reanut')
     do while (words(1) /= 'NUMER')
        call ecoute('hlm_reanut')
     enddo
     !
     ! Begin to read data
     !
     do while (words(1) /= 'ENDNU')
        call ecoute('hlm_reanut')

        if (words(1) == 'ALGEB') then
           !
           ! Algebraic solver
           !
           call reasol(1_ip)
           !else if (words(1) == 'PRECO') then 
           !  call reasol(2_ip)

        endif
     end do
  end if

end subroutine hlm_reanut
