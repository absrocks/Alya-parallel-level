module def_neutro

  !------------------------------------------------------------------------
  !    
  ! Heading for the incompressible NASTIN routines
  !
  !------------------------------------------------------------------------

  use def_kintyp, only : ip,rp,i1p,i2p,r3p
  use def_kintyp, only : bc_nodes
  use def_kintyp, only : bc_bound
  use mod_ADR,    only : ADR_typ

  !------------------------------------------------------------------------
  !
  ! Parameters
  !
  !------------------------------------------------------------------------

  character(150)           :: &
       fil_rstar_neu
  real(rp),      parameter :: &
       zeneu = epsilon(1.0_rp)

  !--BEGIN REA GROUP
  !------------------------------------------------------------------------
  !
  ! Physical problem: read in neu_reaphy
  !
  !------------------------------------------------------------------------
 
  integer(ip)               :: &
       num_energies_neu,       &       ! Number of energy groups
       num_directions_neu,     &       ! Number of directions
       kfl_icosa_neu,          &       !
       kfl_snord_neu                   !
  real(rp)                  :: & 
       aniso_neu                       ! Linear anisotropic coefficient

  !------------------------------------------------------------------------
  !
  ! Numerical problem: read in neu_reanut
  !
  !------------------------------------------------------------------------

  integer(ip)               :: &
       miinn_neu,              &       ! Maximum inner iterations
       kfl_smobo_neu                   ! B.c. smoothing
  real(rp)                  :: & 
       cotol_neu,              &       ! Tolerance inner iterations       
       relax_neu,              &       ! Relaxation
       nitsche_neu                     ! Nitsche coefficient

  !------------------------------------------------------------------------
  !
  ! Boundary conditions: read in neu_reabcs
  !
  !------------------------------------------------------------------------

  type(bc_nodes), pointer  :: &     
       tncod_neu(:)                  ! Node code type
  type(bc_bound), pointer  :: &     
       tbcod_neu(:)                  ! Boundary code type

  !------------------------------------------------------------------------
  !
  ! Output and Postprocess: read in neu_reaous
  !
  !------------------------------------------------------------------------

  !--END REA GROUP
  !------------------------------------------------------------------------
  !
  ! Others
  !
  !------------------------------------------------------------------------
  !
  ! Dimensions, etc.
  !
  integer(ip)              :: &
       nunkn_neu,             &      ! Number of unknowns
       ncomp_neu,             &      ! Number of components
       nprev_neu,             &      ! Last time step or global iteration
       current_energy_neu,    &      ! Current energy being solved
       current_direction_neu, &      ! Current direction being solved
       kfl_goite_neu                 ! Continue inner iterations
  real(rp)                 :: &
       resid_neu                     ! Residual of outer iterations
  !
  ! Directions
  !
  real(rp),  pointer       :: &
       direc_neu(:,:),        &      ! Directions
       weigd_neu(:),          &      ! Weights of directions
       scattering_neu(:,:)           ! Scattering coefficient
  !
  ! Boundary conditions
  !
  type(i2p),   pointer     :: &
       kfl_fixno_neu(:,:)            ! Nodal fixity 
  type(i1p),   pointer     :: &
       kfl_fixbo_neu(:,:)            ! Element boundary fixity
  type(r3p),   pointer     :: &
       bvess_neu(:,:),        &      ! Essential velocity bc values
       bvnat_neu(:,:)                ! Natural bc values
  !
  ! ADR type
  !
  type(ADR_typ)            :: &
       ADR_NEU                       ! ADR type

end module def_neutro
