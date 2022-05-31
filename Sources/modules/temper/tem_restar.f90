subroutine tem_restar(itask)
  !------------------------------------------------------------------------
  !****f* Temper/tem_restar
  ! NAME 
  !    tem_restar
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... Reads the initial values from the restart file
  !            2 ... Writes restart file
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_kermod,         only : kfl_adj_prob
  use mod_tem_arrays,     only : tem_arrays
  use mod_communications, only : PAR_BROADCAST
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,icomp,kfl_gores
  
  !----------------------------------------------------------------------
  !
  ! Primary arrays
  !
  !----------------------------------------------------------------------

  if( itask == ITASK_READ_RESTART ) then
     call tem_arrays('READ RESTART')
  else if( itask == ITASK_WRITE_RESTART ) then
     call tem_arrays('WRITE RESTART')
  end if

  !----------------------------------------------------------------------
  !
  ! Restart file tem.rst
  !
  !----------------------------------------------------------------------

  if( itask == ITASK_READ_RESTART ) then
     !
     ! Read restart
     !
     if( INOTSLAVE ) then
        read(momod(modul) % lun_rstar) vinvt_tem(1:4)
        read(momod(modul) % lun_rstar) prthe(1:4)
        read(momod(modul) % lun_rstar) avtim_tem
     end if
     if( IPARALL ) then
        call PAR_BROADCAST(4_ip,vinvt_tem)
        call PAR_BROADCAST(4_ip,prthe)
        call PAR_BROADCAST(avtim_tem)
     end if

  else if( itask == ITASK_WRITE_RESTART) then 
     !
     ! Write restart file
     !
     if( INOTSLAVE ) then
        write(momod(modul) % lun_rstar) vinvt_tem(1:4)
        write(momod(modul) % lun_rstar) prthe(1:4)
        write(momod(modul) % lun_rstar) avtim_tem
     end if
     
  end if

  !----------------------------------------------------------------------
  !
  ! Assign constant tempe forward values for adjoint
  !
  !----------------------------------------------------------------------  
   
  if( itask == ITASK_READ_RESTART ) then
     icomp = min(3_ip,ncomp_tem)  
  else
     icomp = 1 
  end if
  if( itask == ITASK_READ_RESTART .and. kfl_adj_prob == 1 ) then
    do ipoin = 1, npoin
        tempe_forw(ipoin,1) = tempe(ipoin,icomp)    
        tempe(ipoin,icomp)  = 0.0_rp
    end do    
 endif
 
 !----------------------------------------------------------------------
 !
 ! Compute dependent variables not accounted ofr in the restart
 !
 !----------------------------------------------------------------------  
 
 if( itask == ITASK_READ_RESTART ) then
    !
    ! Compute temperature from enthalpy if needed (:,1), even in the case of restart.
    ! This is needed to compute properties
    !
    !call tem_temperature_from_enthalpy()
 end if

end subroutine tem_restar
 
