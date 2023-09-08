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
  implicit none  
  integer(ip) :: ivari,itype,ielem,pelty,pnode,ptopo,inode,ipoin
  real(rp)    :: elcod(ndime,mnode),hleng(ndime)
  !
  ! Distance to wall
  !
  call pts_waldis()
  !
  ! Injection through the boundary
  !
  if( kfl_boundary_injection == 1_ip) call pts_injbou()
  !
  ! ELement charcateristic minimum length
  !
  if( INOTMASTER ) then
     do ielem = 1,nelem
        pelty = abs(ltype(ielem))
        pnode = nnode(pelty)
        ptopo = ltopo(pelty)
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           elcod(1:ndime,inode) = coord(1:ndime,ipoin)
        end do
        call elmgeo_element_characteristic_length(& 
             ndime,pnode,elmar(pelty) % dercg,elcod,hleng)
        hleng_pts(ielem) = hleng(ndime)
     end do
  end if

  if( kfl_rstar_pts == 0 ) then  
     !
     ! Initial solution 
     !
     cutla_pts = 1.0e12_rp
  else
     !
     ! Read restart file
     !
     call pts_restar(READ_RESTART_FILE)
     !cutla_pts = cutim - dtime
  end if
  !
  ! Initialization of the random number generator
  !
  call random_initialization()

end subroutine pts_iniunk
