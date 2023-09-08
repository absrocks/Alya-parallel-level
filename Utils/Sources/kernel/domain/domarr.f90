subroutine domarr(itask)
  !-----------------------------------------------------------------------
  !****f* domain/domarr
  ! NAME
  !    domarr
  ! DESCRIPTION
  !    This routines computes some arrays that depend the mesh.
  !    If the mesh changes, all these arrays should be recomputed.
  !
  !    VMASS ... Lumped mass matrix
  !    VMASC ... Close rule mass matrix
  !    EXNOR ... Exterior normals
  !    YWALP ... Physical nodal wall distance (includes distance to wall)
  !    YWALB ... Physical boundary wall distance (includes distance to wall)
  !    WALLD ... Mesh nodal wall distance 
  !    SKCOS ... Geometrical local basis
  !    COORD ... Fringe coordinates
  !
  !    ITASK = 1 ... Allocate and compute all variables
  !          = 2 ... Recompute all variables except LPOTY
  !
  !    The order of creation of geometrical arrays is the following
  !    because some solver variables must be allocated:
  !   
  !    From domain():
  !    1. domarr(1)    computes: VMASS, VMASC, EXNOR, LPOTY, YWALP, YWALD
  !    2. opebcs()     computes: KGL_GEONO, SKCOS
  !    From Iniunk(): 
  !    3. ker_iniunk() computes: WALLD
  !    
  !    All these variables are then updated together when mesh is changing.
  !    LPOTY is not recomputed as boundary nodes are assumed not to change.
  !
  ! USED BY
  !    Turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_kermod
  use def_domain
  use mod_communications,  only : PAR_GHOST_NODE_EXCHANGE
  use mod_exterior_normal, only : extnor
  use mod_mass_matrix
  implicit none
  integer(ip), intent(in) :: itask
  !
  ! Fringe coordinates
  !
  if( ISLAVE .and. itask == 2 ) call PAR_GHOST_NODE_EXCHANGE(coord,'SUBSTITUTE','IN MY CODE')
  !
  ! VMASS: Diagonal mass matrix using lumping 
  ! 
  call mass_matrix_open_lumped()
  !
  ! VMASC: Diagonal mass matrix using close rule
  !
  call mass_matrix_close()
  !
  ! Consistent mass matrix
  !
  if( itask == 1 ) then
     call mass_matrix_consistent(CONSISTENT_MASS=.true.,CONSISTENT_WEIGHTED_MASS=.false.)
  else
     call mass_matrix_consistent(CONSISTENT_MASS=.true.,CONSISTENT_WEIGHTED_MASS=.true.)
  end if
  !
  ! Check projections
  !
  !call chkmas()  
  !
  ! EXNOR: External normal
  !
  call extnor(itask)  
  !
  ! KFL_GEONO, SKCOS: Geometrical boundary conditions
  !
  if( itask == 2 ) call geonor(itask)
  !
  ! YWALP and YWALB: Physical wall distance on boundary nodes
  !
  call waldis(itask)
  !
  ! WALLD, WALLN: Mesh wall distance and normal. Computed by Kermod: must go through
  !        main subroutine (it uses a solver)
  !
  if( itask == 2 ) call Kermod(-1_ip)
  !
  ! Element bin 
  !
  if( itask == 2 ) call elebin()

end subroutine domarr
