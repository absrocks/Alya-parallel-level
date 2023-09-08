subroutine chm_solsgs()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_solsgs
  ! NAME 
  !    chm_solsgs
  ! DESCRIPTION
  !    This routine solves the SGS equation
  ! USES
  ! USED BY
  !    chm_endite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  use def_chemic
  implicit none
  integer(ip) :: ipoin,iclas
  real(rp)    :: time1,time2

  if( kfl_sgsti_chm /= 0 .or. kfl_stabi_chm >= 1 ) then
     !
     ! Initialization
     !
     call cputim(time1)
     !
     ! Residual projections
     !
     if( INOTMASTER ) then    
        do iclas = iclai_chm,iclaf_chm
           do ipoin = 1,npoin
              rhsid(ipoin) = 0.0_rp
           end do
           !
           ! Update projection
           ! 
           call runend('NOT PROGRAMMED FOR THIS MODEL')
           !
           ! Residual projections
           !
           call rhsmod(1_ip,rhsid)
           do ipoin = 1,npoin
              proje_chm(ipoin,iclas) = rhsid(ipoin) / vmass(ipoin)
           end do
        end do
     end if

     call cputim(time2)

  end if

end subroutine chm_solsgs
