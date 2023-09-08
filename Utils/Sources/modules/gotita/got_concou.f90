subroutine got_concou
!-----------------------------------------------------------------------
!****f* Gotita/got_concou
! NAME 
!    got_concou
! DESCRIPTION
!    This routine checks the Gotita convergence of the run and
!    set the general convergence flags.
! USED BY
!    Gotita
!***
!-----------------------------------------------------------------------
  use      def_master
  use      def_gotita
  use mod_outfor, only : outfor
  implicit none
  !
  ! Check convergence
  !
  if(kfl_conve(modul)==1) then
     if(resiv_got>cotol_got) kfl_gocou = 1
  end if 
  glres(modul) = resiv_got
  !
  ! Output residuals
  !
  coutp(1)='DROPLET VELOCITY'
  routp(1)=resiv_got
  call outfor(9_ip,lun_outpu,' ')
  coutp(1)='WATER VOL. FRAC.'
  routp(1)=resic_got
  call outfor(9_ip,lun_outpu,' ')

end subroutine got_concou
