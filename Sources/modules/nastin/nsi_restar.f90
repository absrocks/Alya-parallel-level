subroutine nsi_restar(itask)
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_restar
  ! NAME 
  !    nsi_restar
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... Reads the initial values from the restart file
  !            2 ... Writes restart file
  ! USES
  ! USED BY
  !    nsi_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_postpr
  use mod_memchk
  use mod_nsi_arrays,     only : nsi_arrays
  use mod_communications, only : PAR_BROADCAST
  use def_kermod,         only : avta1_nsw_ker
  use def_kermod,         only : fact_nsw_ker
  use def_kermod,         only : kfl_wlaav_ker
  use def_kermod,         only : kfl_noslw_ker
  use def_kermod,         only : kount_nsw_ele_ker
  implicit none
  integer(ip), intent(in) :: itask
  
  !----------------------------------------------------------------------
  !
  ! Primary arrays
  !
  !----------------------------------------------------------------------

  if( itask == ITASK_READ_RESTART ) then
     call nsi_arrays('READ RESTART')
  else if( itask == ITASK_WRITE_RESTART ) then
     call nsi_arrays('WRITE RESTART')
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
        read(momod(modul) % lun_rstar) grnor_nsi,ubpre_nsi,avtim_nsi
     end if
     if( IPARALL ) then
        call PAR_BROADCAST(grnor_nsi)
        call PAR_BROADCAST(ubpre_nsi)
        call PAR_BROADCAST(avtim_nsi)
     end if

  else if( itask == ITASK_WRITE_RESTART ) then 
     !
     ! Write restart file
     !
     if( INOTSLAVE ) then
        write(momod(modul) % lun_rstar) grnor_nsi,ubpre_nsi,avtim_nsi
     end if
     
  end if
  
end subroutine nsi_restar
 
