!-----------------------------------------------------------------------
!> @addtogroup Master
!> @{
!> @file    Restart_new.f90
!> @author  houzeaux
!> @date    2019-11-26
!> @brief   Restart
!> @details Restart main subroutine
!> @} 
!-----------------------------------------------------------------------

subroutine Restar(itask)

  use def_parame
  use def_master
  use def_domain
  use def_postpr
  use mod_iofile
  use mod_communications, only : PAR_BROADCAST
  use mod_messages,       only : messages_live
  use mod_moduls,         only : moduls
  implicit none
  integer(ip), intent(in) :: itask

  real(rp)                :: dtime_false

  read_restart  = .false.
  write_restart = .false.

  !----------------------------------------------------------------------
  !
  ! Read continue restart file
  !
  !----------------------------------------------------------------------

  if( itask == ITASK_READ_RESTART ) then

     if( kfl_rstar >= 1 ) then
        read_restart  = .true.
        kfl_reawr     = 1
     end if

     if( read_restart ) then

        call messages_live('READ RESTART FILES','START SECTION')
        !
        ! Reading restart files
        !
        do modul = 0,mmodu
           call moddef(7_ip)
        end do
        modul = 0

        if( kfl_rstar == 2 ) then

           if( INOTSLAVE ) then
              read(lun_rstar) ittim,cutim,dtinv,dtime_false  ! time things
              read(lun_rstar) dpthe,prthe(1:4)               ! for lowmach flows
              read(lun_rstar) dtold(1:10)
              read(lun_rstar) dtinv_old(1:10)
           end if

           if( IPARALL ) then
              call PAR_BROADCAST(ittim)
              call PAR_BROADCAST(cutim)
              call PAR_BROADCAST(dtinv)
              call PAR_BROADCAST(dtime_false)
              call PAR_BROADCAST(dpthe)
              call PAR_BROADCAST( 4_ip,prthe)
              call PAR_BROADCAST(10_ip,dtold)
              call PAR_BROADCAST(10_ip,dtinv_old)
           end if

           if( cutim >= abs(timef) ) kfl_gotim = 0
           if( ittim >= mitim      ) kfl_gotim = 0

        end if

        call Kermod(ITASK_READ_RESTART)
        do iblok = 1,nblok
           call moduls(ITASK_READ_RESTART)
        end do

        do modul = 0,mmodu
           call moddef(6_ip)
        end do
        modul = 0

        call messages_live('READ RESTART FILES','END SECTION')

     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Write restart file
  !
  !----------------------------------------------------------------------

  if( itask == ITASK_WRITE_RESTART ) then

     if(  kfl_preli        == 1                      .and. ( &
          mod(ittim,nprit) == 0                      .or.    &
          cutim            >= abs(timef)-1.0e-10_rp  .or.    &
          ittim            >= mitim                  .or.    &
          kfl_timei        == 0                      .or.    &
          kfl_stop         /= 0                            ) ) then
        kfl_reawr     = 2
        write_restart = .true.
     end if

     if( write_restart ) then

        call messages_live('WRITE RESTART FILES','START SECTION')
        do modul = 0,mmodu
           call moddef(8_ip)
        end do
        modul = 0

        if( INOTSLAVE ) then

           if( kfl_timei == 0 ) then
              write(lun_rstar) 0_ip,0.0_rp,0.0_rp,0.0_rp
           else
              write(lun_rstar) ittim,cutim,dtinv,dtime
           end if
           write(lun_rstar) dpthe,prthe(1:4)
           write(lun_rstar) dtold(1:10)
           write(lun_rstar) dtinv_old(1:10)
        end if

        call Kermod(ITASK_WRITE_RESTART)
        do iblok = 1,nblok
           call moduls(ITASK_WRITE_RESTART)
        end do

        do modul = 0,mmodu
           call moddef(6_ip)
        end do
        modul = 0

        call messages_live('WRITE RESTART FILES','END SECTION')

     end if

  end if
  !
  ! Close restart file
  !
  kfl_reawr = 0             
  call moddef(6_ip)

end subroutine Restar
