subroutine tur_restar(itask)
  !------------------------------------------------------------------------
  !****f* turbul/tur_restar
  ! NAME 
  !    tur_restar
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
  use def_turbul
  use mod_postpr
  use mod_memchk
  use def_kermod, only     :  kfl_adj_prob
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: icomp,iwopo,kfl_gores,iturb,ipoin
  !
  ! Check if restrt file should be read or written
  !
  call respre(itask,kfl_gores)
  if( kfl_gores == 0 ) return

  if( itask == READ_RESTART_FILE ) then
     icomp = 3
  else
     icomp = 1
  end if

  !----------------------------------------------------------------------
  !
  ! Turbulent variables and turbulent viscosity
  !
  !----------------------------------------------------------------------

  iwopo = 24
  do iturb = 1,nturb_tur
     iwopo = iwopo + 1
     if( INOTMASTER ) gesca => untur(iturb,1:npoin,icomp)
     call postpr(gesca,postp(1)%wopos(1:2,iwopo),ittim,cutim)
  end do
 
  iwopo = 6
  call postpr(turmu,postp(1)%wopos(1:2,iwopo),ittim,cutim)
    
  !----------------------------------------------------------------------
  !
  ! Distance to the wall!  in version 796 this did not make sense. wopos(1) is empty so it resulted in CANNOT OPEN RESTART FILE: pet-.rst
  ! Moreover even in previous version I do not understand we we were reading it from a restart file instead of recalculating it since one might want to change it 
  ! during a restart 
  !
  !----------------------------------------------------------------------

!  if( kfl_walld_tur == 1 .and. kfl_walgo_tur == 1 ) then
!     iwopo = 10
!     if( IMASTER ) walld_tur => nul1r
!     call postpr(walld_tur,postp(1)%wopos(1:2,iwopo),ittim,cutim)    
!  end if

  !----------------------------------------------------------------------
  !
  ! Projections for OSS stabilization
  !
  !----------------------------------------------------------------------

  if( kfl_ortho_tur >= 1 ) then

     iwopo = 28
     do iturb = 1,nturb_tur
        iwopo = iwopo + 1
        if( INOTMASTER ) gesca => unpro_tur(iturb,1:npoin)
        call postpr(gesca,postp(1)%wopos(1:2,iwopo),ittim,cutim)
     end do
     if (kfl_ortho_tur==2) then
        do iturb = 1,nturb_tur
           iwopo = iwopo + 1
           if( INOTMASTER ) gesca => unprr_tur(iturb,1:npoin)
           call postpr(gesca,postp(1)%wopos(1:2,iwopo),ittim,cutim)
        end do
     end if
 
  end if
!!$  if (kfl_shock_tur/=0) then
!!$     do iturb = 1,nturb_tur
!!$        iwopo = iwopo + 1
!!$        if( INOTMASTER ) gesca => unpgr_tur(iturb,1:npoin)
!!$        call postpr(gesca,postp(1)%wopos(1:2,iwopo),ittim,cutim)
!!$     end do    
!!$  end if

!!$  
  
  !----------------------------------------------------------------------
  !
  ! assign constant tempe forward values for adjoint
  !
  !----------------------------------------------------------------------  
   
!   if( itask == READ_RESTART_FILE .and. kfl_adj_prob == 1 ) then
!     do iturb = 1,nturb_tur
!       do ipoin = 1, npoin
! 	untur_forw(iturb,ipoin,1) = untur(iturb,ipoin,icomp)
! 	untur_forw(iturb,ipoin,2) = untur(iturb,ipoin,icomp)
! 	untur(iturb,ipoin,icomp) = 0.0_rp
!       end do    
!     enddo
!   endif
  
  !----------------------------------------------------------------------
  !
  ! Finish
  !
  !----------------------------------------------------------------------

  call respre(3_ip,kfl_gores)

end subroutine tur_restar
 
