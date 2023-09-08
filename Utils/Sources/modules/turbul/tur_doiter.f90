subroutine tur_doiter
!-----------------------------------------------------------------------
!****f* Turbul/tur_doiter
! NAME 
!    tur_doiter
! DESCRIPTION
!    This routine controls the internal loop of the temperature equation.
! USES
!    tur_begite
!    tur_solite
!    tur_endite
! USED BY
!    Turbul
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_turbul
  implicit none
  
  
  if( kfl_stead_tur == 0 ) then
     call tur_begite()

     do while( kfl_goite_tur == 1 )
        if( kfl_algor_tur == 1 ) then
           do iunkn_tur = 1,nturb_tur
              do itera_tur =1,niter_tur
                 call tur_solite()
                 call tur_endite(one)
              end do              
              ! acualizes untur(:,2) with relaxation factor
              call tur_updunk(12_ip)
           end do
        else
           call tur_solite()
           call tur_endite(one)
        end if
     end do
     call tur_endite(two)
  end if

end subroutine tur_doiter
