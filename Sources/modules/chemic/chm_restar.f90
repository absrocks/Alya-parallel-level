subroutine chm_restar(itask)
  !------------------------------------------------------------------------
  !****f* Chemic/chm_restar
  ! NAME 
  !    chm_restar
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... Reads the initial values from the restart file
  !            2 ... Writes restart file
  ! USES
  ! USED BY
  !    lev_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_chm_arrays
  use mod_communications, only : PAR_BROADCAST
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: icomp,iwopo,kfl_gores,iclas
  integer(ip)             :: dummi(10), nreal, ipoin
  character(5)            :: wopos(3), outna

  !----------------------------------------------------------------------
  !
  ! Arrays
  !
  !----------------------------------------------------------------------

  if( itask == ITASK_READ_RESTART ) then
     call chm_arrays('READ RESTART')
  else if( itask == ITASK_WRITE_RESTART ) then
     call chm_arrays('WRITE RESTART')
  end if

  !----------------------------------------------------------------------
  !
  ! Variables
  !
  !---------------------------------------------------------------------- 

  if( itask == ITASK_READ_RESTART ) then
     !
     ! Read restart
     !
     if( INOTSLAVE ) then
        read(momod(modul) % lun_rstar) avtim_chm
     end if
     if( IPARALL ) then
        call PAR_BROADCAST(avtim_chm)
     end if

  else if( itask == ITASK_WRITE_RESTART ) then 
     !
     ! Write restart file
     !
     if( INOTSLAVE ) then
        write(momod(modul) % lun_rstar) avtim_chm
     end if
     
  end if
  
end subroutine chm_restar
 
