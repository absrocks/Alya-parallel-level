subroutine rad_solsgs()
  !-----------------------------------------------------------------------
  !****f* Nastin/rad_solsgs
  ! NAME 
  !    rad_solsgs
  ! DESCRIPTION
  !    This routine solves the SGS equation
  ! USES
  ! USED BY
  !    rad_endite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  use def_radiat
  implicit none
  integer(ip) :: ipoin
  real(rp)    :: time1,time2

  if(  kfl_ortho_rad >= 1 ) then
     !
     ! Initialization
     !
     call cputim(time1)
     resgs_rad(1) = 0.0_rp
     resgs_rad(2) = 0.0_rp
     !
     ! Residual projections
     !
     if( INOTMASTER .and. kfl_ortho_rad /= 0 ) then
        do ipoin = 1,npoin
           rhsid(ipoin) = 0.0_rp
        end do
        !
        ! Update SGS
        !
        call rad_elmope(4_ip)
        !
        ! Residual projections
        !
        call rhsmod(1_ip,rhsid)
        do ipoin = 1,npoin
           tepro_rad(ipoin) = rhsid(ipoin) / vmass(ipoin)
        end do
     end if
     !
     ! Output SGS convergence
     !
     !call rad_cvgsgs()

     call cputim(time2)
     !cputi_rad(3) = cputi_rad(3) + time2 - time1

  end if

end subroutine rad_solsgs
