subroutine setgts(order)
  !-----------------------------------------------------------------------
  !****f* master/setgts
  ! NAME
  !    setgts
  ! DESCRIPTION
  !    This routine computes the time step
  !    ITASK     = 1 ... Initialization (Turnon)
  !              = 2 ... On each time step (Timste)
  !
  !    KFL_TIMCO = 0 ... Prescribed time step
  !              = 1 ... From critical time step
  !              = 2 ... local/per module time step
  ! USES
  ! USED BY
  !    Turnon, Timste
  !***
  !-----------------------------------------------------------------------
  use def_master
  use mod_outfor, only : outfor
  implicit none
  integer(ip), intent(in) :: order

  select case (order)

  case (1_ip)
     !
     ! Compute DTINV and the time for the first time step, as well as the
     ! maximum number of time steps allowed
     !
     dtinv = 0.0_rp
     cutim = timei

     !
     ! Only execute if time is prescribed
     !
     if( kfl_timco == 0 ) call timfun()

     if( dtime >  0.0_rp ) dtinv = 1.0_rp / dtime

     !
     ! This has no sense as time_f can't be smaller than time_i
     !
     if( timef <= timei ) then
        dtime = 0.0_rp
        dtinv = 0.0_rp
     end if

     dtold(1) = dtime
     dtinv_old(1) = 1.0_rp / max(dtold(1),zeror)
     oltim    = 0.0_rp

  case (2_ip)
     !
     ! Use minimum of critical time steps
     !
     dtold(4) = dtold(3)
     dtold(3) = dtold(2)
     dtold(2) = dtold(1)
     dtold(1) = dtime

     !
     ! Only execute if time is prescribed
     !
     if( kfl_timco == 0 ) call timfun()

     !
     ! From critical time step
     !
     if( kfl_timco == 1 ) then
        if( dtinv /= 0.0_rp ) dtime = 1.0_rp / dtinv
     end if

     if( dtold(1) == 0.0_rp ) dtold(1) = dtime
     if( dtold(2) == 0.0_rp ) dtold(2) = dtime
     if( dtold(3) == 0.0_rp ) dtold(3) = dtime
     if( dtold(4) == 0.0_rp ) dtold(4) = dtime
     dtinv_old(1:10) = 1.0_rp / max(dtold(1:10),zeror)
     oltim           = cutim
     cutim           = cutim + dtime
     ioutp(1)        = ittim
     routp(1)        = dtime
     routp(2)        = cutim
     call outfor(15_ip,lun_outpu,' ')

  end select

end subroutine setgts
