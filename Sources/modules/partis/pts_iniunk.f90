!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Partis initial solution
!! @file    pts_iniunk.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!! @brief   This routine determines the initial solution
!! @details Initial solution, allocate memory and read restart file
!>          \verbatim
!>          NTYLA_PTS ... Maximum particle type number
!>          DEPOE_PTS ... Number of deposited particle per element
!>          LEDEP_PTS ... .TRUE. if an eement hosts a deposited particle 
!>          DEPOS_PTS ... Smoothed number of particles per node 
!>                        (computed from DEPOE_PTS in pts_endite)
!>          \endverbatim
!> @} 
!------------------------------------------------------------------------

subroutine pts_iniunk()
  use def_parame
  use def_domain
  use def_master 
  use def_kermod
  use def_partis
  use mod_memory
  use mod_elmgeo, only : elmgeo_element_characteristic_length
  use mod_random, only : random_initialization
  use mod_messages
  implicit none  
  integer(ip) :: ivari,itype,ielem,pelty,pnode,ptopo,inode,ipoin
  real(rp)    :: elcod(ndime,mnode),hleng(ndime)
  !
  ! Injection through the boundary
  !
  if( kfl_boundary_injection == 1_ip) call pts_injbou()

  if( kfl_rstar_pts == 0 ) then  
     !
     ! Initial solution 
     !
     cutla_pts = 1.0e12_rp

  end if
  !
  ! Initialization of the random number generator
  !
  call random_initialization()


end subroutine pts_iniunk
