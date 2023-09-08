subroutine wav_doiter
!-----------------------------------------------------------------------
!****f* Wavequ/wav_doiter
! NAME 
!    wav_doiter
! DESCRIPTION
!    This routine controls the internal loop of the temperature equation.
! USES
!    wav_begite
!    wav_solite
!    wav_endite
! USED BY
!    Temper
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_wavequ
  implicit none

  if(kfl_stead_wav==0) then
     call wav_begite
     do while(kfl_goite_wav==1)
        call wav_solite
        call wav_endite(one)
     end do
     call wav_endite(two)
  end if

end subroutine wav_doiter
