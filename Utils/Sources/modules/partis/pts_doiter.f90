!------------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @name    Partis inner iteration
!> @file    pts_doiter.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   This routine solves a time step
!> @details This routine solves a time step
!>
!>          call pts_begite
!>          call pts_transport_particles
!>             do communication loops
!>                do particles
!>                  call pts_transport_single_particle
!>                end do
!>                call pts_parallelization_migration
!>             end do
!>          call pts_endite
!> @} 
!------------------------------------------------------------------------

subroutine pts_doiter()
  use def_partis
  use mod_pts_transport
  implicit none
  
  call pts_begite()
  if( kfl_version_pts == 1 ) then
     call pts_transport_particles()
  else
     call pts_solite()
  end if
  call pts_endite()

end subroutine pts_doiter
