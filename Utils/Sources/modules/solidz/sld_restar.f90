!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_restar.f90
!> @author  Solidz Team
!> @date    September, 2017
!>          Adds state dependent variables for any model
!> @brief   This routine reads/writes values required for a restart
!> @details Displacement is  always read and written. Velocity and
!>          and acceleration only for dynamic problems.
!>
!>          \verbatim
!>          ITASK = 1 ... Reads the initial values from the restart file
!>                  2 ... Writes restart file
!>          \endverbatim
!> @}
!------------------------------------------------------------------------

subroutine sld_restar(itask)

  use def_kintyp,        only : ip,rp
  use def_master,        only : INOTMASTER
  use def_master,        only : ITER_K
  use def_master,        only : cutim, ittim
  use def_master,        only : displ, postp, gevec, ger3p, READ_RESTART_FILE
  use mod_postpr,        only : postpr
  use def_solidz,        only : nprev_sld
  use def_solidz,        only : kfl_timei_sld, SLD_DYNAMIC_PROBLEM
  use def_solidz,        only : veloc_sld, accel_sld
  use def_solidz,        only : kfl_sdvar_sld, svegm_sld
  use mod_sld_fe2

  implicit none

  integer(ip), intent(in)    :: itask                    !< What to do
  integer(ip)                :: icomp, iwopo, kfl_gores

  !
  ! Check if restart file should be read or written
  !
  call respre(itask,kfl_gores)
  if ( kfl_gores == 0 ) return

  if ( itask == READ_RESTART_FILE ) then
     icomp = nprev_sld
  else
     icomp = ITER_K
  end if
  !----------------------------------------------------------------------
  !
  ! Displacement
  !
  !----------------------------------------------------------------------

  iwopo = 1
  gevec => displ(:,:,icomp)
  call postpr(gevec,postp(1)%wopos(1:3,iwopo),ittim,cutim)

  !----------------------------------------------------------------------
  !
  ! Velocity
  !
  !----------------------------------------------------------------------

  if ( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then

     iwopo = 2
     gevec => veloc_sld(:,:,icomp)
     call postpr(gevec,postp(1)%wopos(1:3,iwopo),ittim,cutim)

  end if

  !----------------------------------------------------------------------
  !
  ! Acceleration
  !
  !----------------------------------------------------------------------

  if ( kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then

     iwopo = 3
     gevec => accel_sld(:,:,icomp)
     call postpr(gevec,postp(1)%wopos(1:3,iwopo),ittim,cutim)

  end if

  !----------------------------------------------------------------------
  !
  ! FE2 read and write restart in the sub-scale
  !
  !----------------------------------------------------------------------

  if (itask == READ_RESTART_FILE .and. INOTMASTER) then
          call fe2_read_restart()
  else if (INOTMASTER) then
          call fe2_write_restart()
  end if

  !----------------------------------------------------------------------
  !
  ! Cauchy stresses
  !
  !----------------------------------------------------------------------
  !<GGU> Adds corresponding flag when they are required
  !iwopo = 4
  !gevec => caust_sld
  !call postpr(gevec,postp(1)%wopos(1:3,iwopo),ittim,cutim)

  !----------------------------------------------------------------------
  !
  ! State Dependent Variables (SDVs)
  !
  !----------------------------------------------------------------------

  if ( kfl_sdvar_sld == 1_ip ) then

     iwopo = 60
     ger3p => svegm_sld(:)
     call postpr(ger3p,postp(1)%wopos(1:3,iwopo),ittim,cutim)

     if ( itask == READ_RESTART_FILE .and. INOTMASTER ) then
        ! Update after reading restart
        call sld_updunk(7_ip)
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Finish
  !
  !----------------------------------------------------------------------

  call respre(3_ip,kfl_gores)

end subroutine sld_restar
