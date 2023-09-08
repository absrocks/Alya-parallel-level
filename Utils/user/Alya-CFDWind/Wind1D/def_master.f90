   module def_master
     !------------------------------------------------------------------------
     !    
     ! Heading for the master routines. This module contains the global varibles
     ! of the code.
     !
     !------------------------------------------------------------------------

!     real*8, parameter :: kar = 0.410d0
     !
     ! Logical units
     !
     integer, parameter :: &
          lun_ini = 99, lun_rst1 = 131, lun_rst2 = 132, lun_postp= 122, &
          lun_globa = 123, lun_timre = 124,  mcuts = 10, lun_lad = 141
     real, parameter :: &
          pi = 3.14159265358979323846d0,  &
!          rhocp = 1234.5    ! rho 1.229 cp =1004.5          
          rhocp = 1232.9,  & ! rho 1.225 cp =1006.43 (Allinot-Masson)
          gravi = 9.81d0,  & ! gravity acceleration
          keyam = 1.0d-4,  & ! 1.0d-4 ambient key
          lmamb  = 1.0       ! ambient mixing length

     real(8)        :: rad_time(1000), rad_heat(1000) ! time and rad heat
     real(8)        :: ug_time(1000),  u_geos(1000,2), johans(216,56) ! time and geost velocity
     real(8)        :: epsam = 7.208d-8 ! 1.4d-9 !7.208d-8   ! ambient eps   
     
     integer ::  lun_conve(5),  lun_cutre(mcuts)
     data  lun_conve  /101, 102, 103, 104, 105/
     data  lun_cutre  /111, 112, 113, 114, 115, 116,117, 118, 119, 120/ !cut results
     
     integer :: igaus, inode, ielem,  itime, ipoin, iiter, &
          kfl_goite, izdom, jnode, jpoin, istep, icuts, &
          kfl_close, & ! close integration rule
          kfl_order, & ! finite element interpolation order
          kfl_model, & ! k eps model: 0:Limited length,1: RNG,2: Realizable
          kfl_thmod, & ! thermal coupling model
          kfl_case,  & ! Transient case
          kfl_temp   ! temperature for the wall law  1:T2m   2:TSK
          
     ! REAL POINTERS

     real*8,pointer  ::  &    ! GLOBAL
          unkno(:),      &    ! unknown to the solver
          veloc(:,:,:),  &    ! Velocity  NSTINC
          tempe(:,:),  &      ! potential temperature
          keyva(:,:),    &    ! unkno: turbulent kinetic energy
          epsil(:,:),    &    ! unkno: turbulent dissipation
          veloc_ini(:,:),&    ! unkno: turbulent dissipation
          keyva_ini(:),  &    ! unkno: turbulent dissipation
          epsil_ini(:),  &    ! unkno: turbulent dissipation
          tempe_ini(:),  &    
          amatr(:),      &    ! solver matrix
          avect(:),      &
          diago(:),      &
          cvect(:),      &
          rhsid(:),      &
          shape(:,:),    & 
          deriv(:,:),    & 
          deri2(:,:),    & 
          coord(:),      &
          weigp(:),      &
          posgp(:)

     ! Tendencies 
     real*8,pointer  ::  &    !Advections and Geostrophic wind
          tempe_adv(:),  &
          u_adv(:),      &
          v_adv(:),      &
          u_geo(:),      &
          v_geo(:),      &
          u_meso(:),     &
          v_meso(:),     &
          th_meso(:)
     
     ! Model Coefficients
     real*8 ::   &    ! Model coefficients
          c1, c2, cmu, sigka, sigep, kar, cmu0, sigte, &
          !monin obukhov coefss
          mo_prand, mo_gamm1, mo_gamm2, mo_beta1, mo_beta2
     real*8  ::  &    ! GLOBAL
          dtinv, ugeos(2), vegeo, ustar, rough, densi,ustar2,  &
          dwall, toler(5), fcori, l_max, ctime, length, dz1, &
          hflx0, ztsbl, lmoni, alpha_mo, beta_mo, expon_mo, a(6,4), &
          tewal,t2m,qw,logft, &   ! wall temp, for transient temp problem, it should be a transient function
          cutpr(mcuts), &  ! cuts to be ploted
          lenmy, &         ! length mellor yamada
          tetop, &            ! temperature at top
          cdcan, LAD,  heica, LAI, ustar_Can, safet, &
          ztemin, gradbo, gradto, & 
          teref = 283.15d0    ! reference temperature

     !numerical variables
     real*8  ::  &
          damping, &
          z_damping
     

     integer*4, pointer  :: &
          lnods(:,:),       &
          ia(:),            &
          ja(:)
     
     integer*4 ::       &
          nzdom,        &     !3*(npoin-2) + 4 = 3*npoin -2  
          npoin,        &     !number of nodes
          nelem,        &     !number of element s
          ngaus,        &     !number of integration points (usually 2)       
          nnode,        &     !number of nodes per element (usually 2)
          nstep,        &     !maximum number of time steps of the run
          kfl_bouco_vel,&     !kind of veloc bottom boundary condition
          kfl_topco_vel,&     !kind of veloc top boundary condition
          maxit(4),     & 
          miite,        &
          ncuts,        &     !number od heiught cuts to postprocess
          stepr,        &     !profiles each stepr steps
          ielec(mcuts), &     !element to who icut belongs   
          kfl_canmo,    &     !canopy model
          kfl_candi           !canopy distribution
     logical ::   &
          kfl_thcou, &        !thermal coupling  
          kfl_trtem, &        !transient temperature
          kfl_canop, &        !canopy 
          kfl_logva, &        !logarithmic variables 
          kfl_abl2 , &        !ABL2 wall law
          kfl_local, &        !Local time step
          kfl_thadv, &        !Temperature advection from meso
          kfl_momadv, &       !Momentum advection from meso
          kfl_pressgr, &      !Pressure gradient from meso
          kfl_mom_nudging, &  !Momentum nudging to meso
          kfl_tem_nudging     !Potential temperature nudging to meso
     
     !
     !*** Restart variables
     !
     logical                        :: restart_in
     logical                        :: restart_out
     integer(4)                     :: freq_rst
     integer(4)                     :: nz_rst   
     real(8)                        :: fcori_rst
     real(8)                        :: ust_rst
     real(8)                        :: hflux_rst
     real(8)                        :: tstar_rst
     real(8)                        :: tewal_rst
     real(8)                        :: lmax_rst
     real(8),allocatable            :: z_rst(:)
     real(8),allocatable            :: u_rst(:)
     real(8),allocatable            :: v_rst(:)
     real(8),allocatable            :: key_rst(:)
     real(8),allocatable            :: eps_rst(:)
     real(8),allocatable            :: temp_rst(:)

     !
     !*** Netcdf input variables
     !
     character(len=120)             :: meso_file
     integer(4)                     :: nt,nz,np,lev_ref
     real(8)                        :: fc
     real(8)                        :: lat
     real(8)                        :: lon
     real(8),allocatable            :: ustini(:)  !nt
     real(8),allocatable            :: times(:)   !nt
     real(8),allocatable            :: height(:)  !nz
     real(8),allocatable            :: uini(:)    !nz
     real(8),allocatable            :: vini(:)    !nz
     real(8),allocatable            :: thini(:)   !nz
     real(8),allocatable            :: u(:,:)     !nz,nt
     real(8),allocatable            :: v(:,:)     !nz,nt    
     real(8),allocatable            :: th(:,:)    !nz,nt 
     real(8),allocatable            :: ugeo(:,:)  !nz,nt
     real(8),allocatable            :: vgeo(:,:)  !nz,nt
     real(8),allocatable            :: uadv(:,:)  !nz,nt
     real(8),allocatable            :: vadv(:,:)  !nz,nt
     real(8),allocatable            :: thadv(:,:) !nz,nt
     real(8),allocatable            :: t2(:)      !nt
     real(8),allocatable            :: tsk(:)     !nt
     real(8),allocatable            :: hflux(:)   !nt
     !Diagnosed variables
     Integer(4),allocatable         :: seconds(:) !nt   

     !
     !*** netCDF output variables
     !
     logical                        :: kfl_netCDF       !netCDF output
     logical                        :: append_netCDF    !append variables to existing file
     integer(4)                     :: freq_netcdf      !write frequency

     character(len=120)             :: rad_file   ! time dependent radiation
     character(len=120)             :: vegeo_file ! time dependent geostrophic veloc   

     
    end module def_master
