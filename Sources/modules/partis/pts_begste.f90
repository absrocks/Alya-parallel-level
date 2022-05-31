!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_begste.f90
!> @author  houzeaux
!> @date    2018-09-22
!> @brief   Start a time step
!> @details Inject particles if required
!> @} 
!-----------------------------------------------------------------------

subroutine pts_begste()
  
  use def_parame
  use def_master
  use def_domain
  use def_partis
  use mod_physics,       only : physics_set_liquid_pressure
  use mod_pts_injection, only : pts_injection_injectors
  use mod_pts_thermodynamic
  implicit none
  integer(ip) :: idime,kk,itype
  real(rp) :: xx

  !
  ! Update pressure dependendent particle properties if thermodynamic pressure changes
  !
  do itype = 1,mtyla
     if (( parttyp(itype) % kfl_exist == 1 ) .and. &
         ( parttyp(itype) % kfl_therm /= 0 ) .and. &
         ( kfl_prthe /= 0 )                 ) then
        
        call physics_set_liquid_pressure(parttyp(itype) % liq, prthe(1))

     endif
  enddo


  if( kfl_timco == 1 ) dtinv = max(dtinv,1.0_rp/dtime)
  !
  ! Inject particles
  !
  nlagr = nlacc_pts
  call pts_injection_injectors()
  !
  ! injetion happens before timestep is updated, so save injection into the old file still
  !
  ! If pts.res splitting is requested, open the new file
  !
  if(kfl_ptsres_split>0) then
     if( mod(ittim, kfl_oufre_pts) == 0 ) then
        call pts_openfi(2_ip)
     end if
  end if

end subroutine pts_begste
