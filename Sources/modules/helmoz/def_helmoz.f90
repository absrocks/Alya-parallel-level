module def_helmoz

  !------------------------------------------------------------------------
  ! Sources/modules/helmoz/def_helmoz.f90
  ! NAME 
  !    def_helmoz
  ! DESCRIPTION
  !    Heading for the routines in the 'helmoz' module.
  ! USES
  ! USED BY
  !    Almost all
  !------------------------------------------------------------------------

  use def_kintyp

  !------------------------------------------------------------------------
  !
  ! Parameters and types
  !
  !------------------------------------------------------------------------
  integer(ip), parameter  :: &
       lun_pdata_hlm = 2001, lun_outpu_hlm = 2002, lun_conve_hlm = 2003,&
       lun_solve_hlm = 2004, lun_rstar_hlm = 2005, lun_setse_hlm = 2006,&
       lun_setsb_hlm = 2007, lun_setsn_hlm = 2008, lun_cvgso_hlm = 2011

  character(150)          :: fil_rstar_hlm

  !------------------------------------------------------------------------
  !
  ! Physical problem: read in 'hlm_reaphy' file
  !
  !------------------------------------------------------------------------
  integer(ip)                           :: &
       emmet_hlm,                         &    ! Controlled-source (CSEM, 1-3) or natural source (MT, 4-5) electromagnetic method
       ppcod_hlm,                         &    ! Way to read the potentials (1 - read and interpolate in the r-z plane, 2 - read directly on the nodes, 3 - generate them by itself)
       ppout_hlm,                         &    ! Output fields (1 - primary, 2 - secondary, 3 - total)
       nequs_hlm,                         &    ! Number of equations (number of unknown parameters: Asx, Asy, Asz, Psis)
       ncond_hlm,                         &    ! Number of elements in conductivity tensor
       nsite_hlm,                         &    ! Number of sites
       nmlsi_hlm,                         &    ! Number of closest mesh nodes to a site needed for the MLSI
       aniso_hlm,                         &    ! Anisotropy level: 0 - isotropy, 1 - TIV (Gx = Gy /= Gz), 2 - TIH (Gx /= Gy = Gz), 3 - Triaxial (Gx /= Gy /= Gz)
       nshot_hlm,                         &    ! Number of shots
       kfl_edges_hlm                           ! Edge element formulation
  real(rp)                              :: &
       frequ_hlm,                         &    ! f = Frequency
       anguf_hlm,                         &    ! w = Angular frequency
       xoffs_hlm,                         &    ! X-offset of the source
       yoffs_hlm,                         &    ! Y-offset of the source
       zoffs_hlm,                         &    ! Z-offset of the source
       elcur_hlm,                         &    ! I = Electric current
       length_hlm,                        &    ! dl = Dipole length
       airpl_hlm                               ! Air-plane that separates MLSI regions
  real(rp),pointer                      :: &
       xoffsv_hlm(:),                         &    ! X-offset of the source 
       yoffsv_hlm(:),                         &    ! Y-offset of the source  
       zoffsv_hlm(:),                         &    ! Z-offset of the source  
       elcurv_hlm(:),                         &    ! I = Electric current    
       lengthv_hlm(:)                              ! dl = Dipole length      
  real(rp), pointer                     :: &       
       perma_hlm(:),                      &    ! mu = Magnetic permeability
       epsil_hlm(:),                      &    ! eps = Dielectric permittivity
       sigma_hlm(:,:),                    &    ! sig = Electric conductivity tensor
       dsigma_hlm(:,:),                   &    ! dsig = Electric conductivity tensor - Background electric conductivity tensor
       bckco_hlm(:)                            ! bcksig = Background electric conductivity tensor
  complex(rp), pointer                  :: &
       pvepo_hlm(:,:)                          ! Values of primary vector potential in nodes of a reference mesh of size nr * nz, in the r-z plane
  integer(ip)                           :: &
       nr_hlm,                            &    ! Number of points in the r direction, nr
       nz_hlm                                  ! Number of points in the z direction, nz
  real(rp), pointer                     :: &    
       r_hlm(:),                          &    ! r cylindrical coordinates of nodes of a reference mesh in the r-z plane
       z_hlm(:)                                ! z cylindrical coordinates of nodes of a reference mesh in the r-z plane
  real(rp), pointer                     :: &    
       site_hlm(:,:)                           ! Cartesian coordinates of sites
  !complex(rp), pointer ::  selsp_obs(:), smgvpX_obs(:), smgvpY_obs(:), smgvpZ_obs(:)
  complex(rp), pointer :: selsp_obs(:,:), smgvpX_obs(:,:), smgvpY_obs(:,:), smgvpZ_obs(:,:)
  integer(ip), pointer :: countobs(:)
  real(rp),    pointer :: sum_obs(:)
  real(rp),    pointer :: weightv_hlm(:)
  real(rp),    pointer :: costfv_hlm(:)

  integer(ip), pointer :: incidence_obs(:,:)

  !------------------------------------------------------------------------
  !
  ! Numerical problem: read in 'hlm_reanut' file
  !
  !------------------------------------------------------------------------

  real(rp)                              :: &
       cotol_hlm                               ! Convergence tolerance

  !------------------------------------------------------------------------
  !
  ! Output and Postprocess: read in 'hlm_reaous' file
  !
  !------------------------------------------------------------------------

  !------------------------------------------------------------------------
  !
  ! Boundary conditions: read in 'hlm_reabcs' file
  !
  !------------------------------------------------------------------------

  type(bc_nodes), pointer               :: &
       tncod_hlm(:)                            ! Node code type
  type(bc_nodes), pointer               :: &
       tgcod_hlm(:)                            ! Geometrical node code type
  type(bc_nodes), pointer               :: &
       tecod_hlm(:)                            ! Edge code type
  type(bc_bound), pointer               :: &
       tbcod_hlm(:)                            ! Boundary code type

  integer(ip), pointer                  :: &
       kfl_fixno_hlm(:,:)                      ! Nodal boundary conditions

  !------------------------------------------------------------------------
  !
  ! Others
  !
  !------------------------------------------------------------------------
  integer(ip)                           :: &
       ittot_hlm,                         &    ! Total number of iterations
       ivari_hlm,                         &    ! Current variable
       cntrl_hlm                               ! Postprocess control
  real(rp)                              :: &
       resid_hlm                               ! Residual
  integer(ip), pointer                  :: &    
       clnod1_hlm(:)                           ! Global numbers of N closest mesh nodes to a test point, for all test points
  real(rp), pointer                     :: &
       clnod2_hlm(:),                     &    ! Distances between each of the N closest mesh nodes and a test point, for all test points
       clcoor_hlm(:,:)                         ! Cartesian coordinates of N closest mesh nodes to a test point, for all test points

  integer(ip), pointer                  :: &    
       clsite1_hlm(:)                          ! Global numbers of N closest mesh nodes to a site, for all sites
  real(rp), pointer                     :: &
       clsite2_hlm(:),                    &    ! Distances between each of the N closest mesh nodes and a site, for all sites
       clsite_hlm(:,:)                         ! Cartesian coordinates of N closest mesh nodes to a site, for all sites



  complex(rp), pointer                  :: &
       smgvp_hlm(:,:),                    &    ! Secondary magnetic vector potential
       pmgvp_hlm(:,:),                    &    ! Primary magnetic vector potential
       selsp_hlm(:),                      &    ! Secondary electric scalar potential
       pelsp_hlm(:),                      &    ! Primary electric scalar potential
       elefi_hlm(:,:),                    &    ! Vector of electric field intensity
       magfi_hlm(:,:)                          ! Vector of magnetic field intensity

  real(rp), pointer                  :: &
       diffj_hlm(:),                   &       ! gradient from CSEM inversion to post-processor
       design_hlm(:)                           ! design variables from CSEM inversion to post-processor
  !
  ! Edge element stuffs
  !
  real(rp), pointer                  :: &
       sign_edges_hlm(:,:),             &      ! Sign of the edges
       length_edges_hlm(:)

end module def_helmoz
