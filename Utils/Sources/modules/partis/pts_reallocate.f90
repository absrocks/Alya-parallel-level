!------------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_injec.f90
!> @author  Guillaume Houzeaux
!> @brief   Reallocate particle type
!> @details Reallocate particle type LAGRTYP
!> @} 
!------------------------------------------------------------------------

subroutine pts_reallocate(klagr)
  use def_master,   only : ip
  use def_partis,   only : latyp
  use def_partis,   only : lagrtyp
  use def_partis,   only : mlagr
  use def_partis,   only : lagrtyp_init
  use def_partis,   only : lagrtyp_tmp
  use mod_messages, only : livinf
  implicit none  
  integer(ip), intent(in) :: klagr
  integer(ip) :: memorysize_for_particles
  type(latyp) :: particle_structure   !used to calculate the amount of memory allocated
  integer(4)  :: istat
  if( size(lagrtyp) > klagr ) then
     !
     ! Required size is lower than actual
     !
     return

  else
     !
     ! Reallocate
     !
     allocate( lagrtyp_tmp(klagr), stat=istat )
     if(istat/=0) then
        memorysize_for_particles = 0_8!sizeof(particle_structure)*klagr
        call livinf(-7_ip,"Partis failed to reallocate storage for particles (lagrtyp). Requested memory in bytes",memorysize_for_particles)
        call livinf(-7_ip,"Number of requested particles",klagr)
        call runend("pts_reallocate: Particle storage allocation failure")
     end if

     lagrtyp_tmp(1:mlagr)       =  lagrtyp(1:mlagr)
     lagrtyp_tmp(mlagr+1:klagr) =  lagrtyp_init
     mlagr                      =  klagr
     deallocate(lagrtyp)
     lagrtyp                    => lagrtyp_tmp

  end if

end subroutine pts_reallocate
