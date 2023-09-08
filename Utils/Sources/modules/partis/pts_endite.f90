!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Partis inner iteration
!! @file    pts_doiter.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!! @brief   This routine ends a time step for particles
!! @details Update residence time
!> @} 
!------------------------------------------------------------------------

subroutine pts_endite()
  use def_kintyp
  use def_master
  use def_domain
  use def_partis
  use mod_communications, only : PAR_FROM_GHOST_ELEMENT_EXCHANGE
  use mod_communications, only : PAR_FROM_GHOST_BOUNDARY_EXCHANGE
  use mod_pts_transport,  only : pts_transport_finalize
  implicit none
  integer(ip) :: ielem,iboun
  !
  ! Compute some numbers and output some information
  !
  if( kfl_version_pts == 1 ) call pts_transport_finalize()
  !
  ! Residence: pass information from ghost element
  ! Pass information computed on my ghost element and send it to my neighbors
  !
  if( INOTMASTER .and. kfl_resid_pts /= 0 ) then
     call PAR_FROM_GHOST_ELEMENT_EXCHANGE(resid_pts,'SUM','IN MY CODE')
     do ielem = nelem+1,nelem_2
        resid_pts(1:ntyla_pts,ielem) = 0.0_rp
     end do 
  end if
  !
  ! Boundary deposition on boundaries
  !
  if( INOTMASTER .and. kfl_depos_pts /= 0 ) then
     call PAR_FROM_GHOST_BOUNDARY_EXCHANGE(depob_pts,'SUM','IN MY CODE')
     do iboun = nboun+1,nboun_2
        depob_pts(1:ntyla_pts,iboun) = 0.0_rp
     end do 
  end if

end subroutine pts_endite
