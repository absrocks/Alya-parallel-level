!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Starts an iteration
!! @file    pts_begite.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!! @brief   This routine starts an iteration
!! @details The different tasks carried out are:
!>          - Initialize the array DEPOE_PTS : deposited particles
!>            per element
!> @}
!------------------------------------------------------------------------

subroutine pts_begite()
  use def_kintyp
  use def_master
  use def_domain
  use def_kermod, only : kfl_vefun
  use def_partis
  use mod_gradie
  implicit none
  integer(ip) :: ielem,itype,ilagr
  

  !
  ! Initialize deposition map
  !
  if( INOTMASTER .and. kfl_depos_pts /= 0 ) then
     do ielem = 1,nelem
        do itype = 1,ntyla_pts
           if( parttyp(itype) % kfl_exist /= 0 ) then
              depoe_pts(itype,ielem) = 0.0_rp
           end if
        end do
     end do
  end if
  !
  ! Take off particles that have dissapeared in order to free the space
  ! to save more particles in the future
  !
  if( INOTMASTER ) then
     do ilagr = 1,mlagr
        if( lagrtyp(ilagr) % kfl_exist <= -2 ) then
           lagrtyp(ilagr) % kfl_exist = 0
        end if
     end do
  end if
  !
  ! Saffman force: compute velocity deformation tensor
  !
  if( INOTMASTER ) then
     itype_loop: do itype = 1,ntyla_pts
        if( parttyp(itype) % kfl_exist /= 0 .and. parttyp(itype) % kfl_saffm /= 0 ) then
           call graten(advec,defor_pts)
           if (kfl_vefun /= 0_ip ) then
              call vortic_advec(-1_ip)
           else
              call vortic(-1_ip)
           end if
           exit itype_loop
        end if
     end do itype_loop
  end if
  !
  ! Momentum is accumulated for nastin
  !
  if( kfl_coupl(ID_NASTIN,ID_PARTIS) /= 0 .and. INOTMASTER ) then
     momen = 0.0_rp
  end if


end subroutine pts_begite

