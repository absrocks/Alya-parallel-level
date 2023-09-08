subroutine wav_endite(itask)
!-----------------------------------------------------------------------
!****f* Wavequ/wav_endite
! NAME 
!    wav_endite
! DESCRIPTION
!    This routine checks convergence and performs updates of the
!    temperature  at:
!    - itask=1 The end of an internal iteration
!    - itask=2 The end of the internal loop iteration
! USES
!    wav_cvgunk
!    wav_updunk
! USED BY
!    wav_doiter
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_wavequ
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: itask

  select case(itask)

  case(1)
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||) and update unknowns:
     !  u(n,i,j-1) <-- u(n,i,j) 
     !
     call wav_cvgunk(1_ip)
     call wav_updunk(3_ip)

  case(2)
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||) and update unknowns:
     !  u(n,i-1,*) <-- u(n,i,*) 
     !
     call livinf(16_ip,' ',itinn(modul))
     call wav_cvgunk(2_ip)
     call wav_updunk(4_ip)

  end select

end subroutine wav_endite
