subroutine ker_allaws()

  !-----------------------------------------------------------------------
  !****f* kermod/ker_allaws
  ! NAME
  !   ker_allaws
  !   lresp = -1 ... Update property always
  ! DESCRIPTION
  !   Define laws
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_master
  use def_kermod
  use mod_ker_polynomial, only : ker_polynomial_allaws, ker_polynomial_name !< 2016JUN29   
  implicit none
  integer(ip) :: ilaws
  !
  ! Initialization
  !
  densi_ker % kfl_exist = 0
  visco_ker % kfl_exist = 0
  poros_ker % kfl_exist = 0
  condu_ker % kfl_exist = 0
  sphea_ker % kfl_exist = 0
  dummy_ker % kfl_exist = 0
  turmu_ker % kfl_exist = 0
  absor_ker % kfl_exist = 0
  scatt_ker % kfl_exist = 0
  mixin_ker % kfl_exist = 0

  densi_ker % kfl_nedsm = 0
  visco_ker % kfl_nedsm = 0
  poros_ker % kfl_nedsm = 0
  condu_ker % kfl_nedsm = 0
  sphea_ker % kfl_nedsm = 0
  dummy_ker % kfl_nedsm = 0
  turmu_ker % kfl_nedsm = 0
  absor_ker % kfl_nedsm = 0
  scatt_ker % kfl_nedsm = 0
  mixin_ker % kfl_nedsm = 0

  do ilaws = 1,mlaws_ker
     densi_ker % llaws(ilaws) % wname     = ''
     poros_ker % llaws(ilaws) % wname     = ''
     visco_ker % llaws(ilaws) % wname     = ''
     condu_ker % llaws(ilaws) % wname     = ''
     sphea_ker % llaws(ilaws) % wname     = ''
     dummy_ker % llaws(ilaws) % wname     = ''
     turmu_ker % llaws(ilaws) % wname     = ''
     absor_ker % llaws(ilaws) % wname     = ''
     scatt_ker % llaws(ilaws) % wname     = ''
     mixin_ker % llaws(ilaws) % wname     = ''

     densi_ker % llaws(ilaws) % lresp     = -2
     poros_ker % llaws(ilaws) % lresp     = -2
     visco_ker % llaws(ilaws) % lresp     = -2
     condu_ker % llaws(ilaws) % lresp     = -2
     sphea_ker % llaws(ilaws) % lresp     = -2
     dummy_ker % llaws(ilaws) % lresp     = -2
     turmu_ker % llaws(ilaws) % lresp     = -2
     absor_ker % llaws(ilaws) % lresp     = -2
     scatt_ker % llaws(ilaws) % lresp     = -2
     mixin_ker % llaws(ilaws) % lresp     = -2

     densi_ker % llaws(ilaws) % where     = ''
     poros_ker % llaws(ilaws) % where     = '' 
     visco_ker % llaws(ilaws) % where     = ''
     condu_ker % llaws(ilaws) % where     = ''
     sphea_ker % llaws(ilaws) % where     = ''
     dummy_ker % llaws(ilaws) % where     = ''
     turmu_ker % llaws(ilaws) % where     = ''
     absor_ker % llaws(ilaws) % where     = ''
     scatt_ker % llaws(ilaws) % where     = ''
     mixin_ker % llaws(ilaws) % where     = ''

     densi_ker % llaws(ilaws) % kfl_gradi =  0
     poros_ker % llaws(ilaws) % kfl_gradi =  0 
     visco_ker % llaws(ilaws) % kfl_gradi =  0
     condu_ker % llaws(ilaws) % kfl_gradi =  0
     sphea_ker % llaws(ilaws) % kfl_gradi =  0
     dummy_ker % llaws(ilaws) % kfl_gradi =  0
     turmu_ker % llaws(ilaws) % kfl_gradi =  0
     absor_ker % llaws(ilaws) % kfl_gradi =  0
     scatt_ker % llaws(ilaws) % kfl_gradi =  0
     mixin_ker % llaws(ilaws) % kfl_gradi =  0

     densi_ker % llaws(ilaws) % kfl_deriv =  0
     poros_ker % llaws(ilaws) % kfl_deriv =  0
     visco_ker % llaws(ilaws) % kfl_deriv =  0
     condu_ker % llaws(ilaws) % kfl_deriv =  0
     sphea_ker % llaws(ilaws) % kfl_deriv =  0
     dummy_ker % llaws(ilaws) % kfl_deriv =  0
     turmu_ker % llaws(ilaws) % kfl_deriv =  0
     absor_ker % llaws(ilaws) % kfl_deriv =  0
     scatt_ker % llaws(ilaws) % kfl_deriv =  0
     mixin_ker % llaws(ilaws) % kfl_deriv =  0
     
     densi_ker % llaws(ilaws) % kfl_grder =  0
     poros_ker % llaws(ilaws) % kfl_grder =  0
     visco_ker % llaws(ilaws) % kfl_grder =  0
     condu_ker % llaws(ilaws) % kfl_grder =  0
     sphea_ker % llaws(ilaws) % kfl_grder =  0
     dummy_ker % llaws(ilaws) % kfl_grder =  0
     turmu_ker % llaws(ilaws) % kfl_grder =  0
     absor_ker % llaws(ilaws) % kfl_grder =  0
     scatt_ker % llaws(ilaws) % kfl_grder =  0
     mixin_ker % llaws(ilaws) % kfl_grder =  0

     turmu_ker % llaws(ilaws) % kfl_deriv_tur =  0
     
     turmu_ker % llaws(ilaws) % kfl_deriv_vel =  0

  end do
  !
  ! Law 1 is constant for all properties
  !
  densi_ker % llaws(1) % wname = 'CONST'
  poros_ker % llaws(1) % wname = 'CONST'
  visco_ker % llaws(1) % wname = 'CONST'
  condu_ker % llaws(1) % wname = 'CONST'
  sphea_ker % llaws(1) % wname = 'CONST'
  dummy_ker % llaws(1) % wname = 'CONST'
  turmu_ker % llaws(1) % wname = 'CONST'
  absor_ker % llaws(1) % wname = 'CONST'
  scatt_ker % llaws(1) % wname = 'CONST'
  mixin_ker % llaws(1) % wname = 'CONST'

  densi_ker % llaws(1) % lresp = -2
  poros_ker % llaws(1) % lresp = -2
  visco_ker % llaws(1) % lresp = -2
  condu_ker % llaws(1) % lresp = -2
  sphea_ker % llaws(1) % lresp = -2
  dummy_ker % llaws(1) % lresp = -2
  turmu_ker % llaws(1) % lresp = -2
  absor_ker % llaws(1) % lresp = -2
  scatt_ker % llaws(1) % lresp = -2
  mixin_ker % llaws(1) % lresp = -2

  densi_ker % llaws(1) % where = 'CONST'
  poros_ker % llaws(1) % where = 'CONST'
  visco_ker % llaws(1) % where = 'CONST'
  condu_ker % llaws(1) % where = 'CONST'
  sphea_ker % llaws(1) % where = 'CONST'
  dummy_ker % llaws(1) % where = 'CONST'
  turmu_ker % llaws(1) % where = 'CONST'
  absor_ker % llaws(1) % where = 'CONST'
  scatt_ker % llaws(1) % where = 'CONST'
  mixin_ker % llaws(1) % where = 'CONST'

  !----------------------------------------------------------------------
  !
  ! Density
  !
  !----------------------------------------------------------------------

  densi_ker % llaws(2) % wname     = 'BIFLU'       ! Bifluid
  densi_ker % llaws(2) % lresp(1)  =  ID_LEVELS  
  densi_ker % llaws(2) % where     = 'IELEM'
  densi_ker % llaws(2) % kfl_gradi = 0             ! Despite there is a gradient we never use it for Biflid flow 

  densi_ker % llaws(3) % wname     = 'LOWMA'       ! Low-Mach rho = p/RT
  densi_ker % llaws(3) % lresp(1)  =  ID_TEMPER
  densi_ker % llaws(3) % lresp(2)  =  ID_NASTIN
  densi_ker % llaws(3) % where     = 'IPOIN'
  densi_ker % llaws(3) % kfl_gradi = 1

  densi_ker % llaws(4) % wname     = 'KLOWM'       ! Low Maxh with mixture rho = pW/RT
  densi_ker % llaws(4) % lresp(1)  =  ID_TEMPER
  densi_ker % llaws(4) % lresp(2)  =  ID_NASTIN
  densi_ker % llaws(4) % lresp(3)  =  ID_CHEMIC
  densi_ker % llaws(4) % where     = 'IPOIN'
  densi_ker % llaws(4) % kfl_gradi = 1
 
  densi_ker % llaws(5) % wname     = 'MIXTU'       ! Mixture of constant density fluids
  densi_ker % llaws(5) % lresp(1)  =  ID_CHEMIC
  densi_ker % llaws(5) % where     = 'IPOIN'
  densi_ker % llaws(5) % kfl_gradi = 1

  densi_ker % llaws(6) % wname     = 'BIPHA'       ! Mixture of constant density fluids in two phases
  densi_ker % llaws(6) % lresp(1)  =  ID_CHEMIC
  densi_ker % llaws(6) % lresp(2)  =  ID_LEVELS
  densi_ker % llaws(6) % where     = 'IPOIN'
  densi_ker % llaws(6) % kfl_gradi = 1

  densi_ker % llaws(7) % wname     = 'TBIPH'       ! Mixture of temperature dependent density fluids in two phases
  densi_ker % llaws(7) % lresp(1)  =  ID_CHEMIC
  densi_ker % llaws(7) % lresp(2)  =  ID_LEVELS
  densi_ker % llaws(7) % lresp(3)  =  ID_TEMPER
  densi_ker % llaws(7) % where     = 'IPOIN'
  densi_ker % llaws(7) % kfl_gradi = 1
 
  densi_ker % llaws(8) % wname     = 'DNAST'       ! Test1
  densi_ker % llaws(8) % lresp(1)  =  ID_NASTAL
  densi_ker % llaws(8) % where     = 'IPOIN'
  densi_ker % llaws(8) % kfl_gradi = 1

  densi_ker % llaws(9) % wname     = 'TEST1'       ! Test1
  densi_ker % llaws(9) % lresp(1)  =  -1
  densi_ker % llaws(9) % where     = 'IELEM'
  densi_ker % llaws(9) % kfl_gradi = 1

  densi_ker % llaws(10) % wname     = 'TEST2'       ! Test2
  densi_ker % llaws(10) % lresp(1)  =  -1
  densi_ker % llaws(10) % where     = 'IPOIN'
  densi_ker % llaws(10) % kfl_gradi = 0
 
  densi_ker % llaws(11) % wname     = 'TEST3'       ! Test3 
  densi_ker % llaws(11) % lresp(1)  =  -1
  densi_ker % llaws(11) % where     = 'IPOIN'
  densi_ker % llaws(11) % kfl_gradi = 1

  densi_ker % llaws(12) % wname     = 'LOWMG'       ! Low mach in Gauss Points
  densi_ker % llaws(12) % lresp(1)  =  ID_TEMPER
  densi_ker % llaws(12) % where     = 'IELEM'
  densi_ker % llaws(12) % kfl_gradi = 0

  densi_ker % llaws(13) % wname     = 'TLOWM'       ! Low-Mach for syncronized CFI combustion model
  densi_ker % llaws(13) % lresp(1)  =  ID_TEMPER    ! density is only updated by temper
  densi_ker % llaws(13) % where     = 'IPOIN'
  densi_ker % llaws(13) % kfl_gradi = 1

  densi_ker % llaws(14) % wname     = 'GKLOW'       ! Low Mach with mixture rho = pW/RT at Gauss Points
  densi_ker % llaws(14) % lresp(1)  =  ID_TEMPER
  densi_ker % llaws(14) % where     = 'IELEM'
  densi_ker % llaws(14) % kfl_gradi = 1
  densi_ker % llaws(14) % kfl_deriv = 1
  densi_ker % llaws(14) % kfl_grder = 1

  call ker_polynomial_allaws(densi_ker, 15_ip, ID_TEMPER, 'IELEM', 0_ip)

  !----------------------------------------------------------------------
  !
  ! Viscosity
  !
  !----------------------------------------------------------------------

  visco_ker % llaws(2) % wname     = 'BIFLU'       ! Bifluid
  visco_ker % llaws(2) % lresp(1)  =  ID_LEVELS
  visco_ker % llaws(2) % where     = 'IELEM'
  visco_ker % llaws(2) % kfl_gradi = 0             ! Despite there is a gradient we never use it for Bifulid flow - see if in nsi_elmop3 grvis

  visco_ker % llaws(3) % wname     = 'SUTHE'       ! Sutherland
  visco_ker % llaws(3) % lresp(1)  =  ID_NASTAL
  visco_ker % llaws(3) % lresp(2)  =  ID_TEMPER
  visco_ker % llaws(3) % where     = 'IPOIN'
  visco_ker % llaws(3) % kfl_gradi = 1
  
  visco_ker % llaws(4) % wname     = 'MUMIX'       ! Mixture of species
  visco_ker % llaws(4) % lresp(1)  =  ID_CHEMIC
  visco_ker % llaws(4) % where     = 'IPOIN'
  visco_ker % llaws(4) % kfl_gradi = 1

  visco_ker % llaws(5) % wname     = 'TEST4'       ! Linear viscosity
  visco_ker % llaws(5) % lresp(1)  =  -1
  visco_ker % llaws(5) % where     = 'IPOIN'
  visco_ker % llaws(5) % kfl_gradi = 1

  visco_ker % llaws(6) % wname     = 'BIFL2'       ! Bifluid
  visco_ker % llaws(6) % lresp(1)  =  ID_LEVELS
  visco_ker % llaws(6) % where     = 'IELEM'
  visco_ker % llaws(6) % kfl_gradi = 1             ! BIFL2 includes gradient

  visco_ker % llaws(7) % wname     = 'GPSUT'       ! Sutherland in Gauss points
  visco_ker % llaws(7) % lresp(1)  =  ID_NASTAL
  visco_ker % llaws(7) % lresp(2)  =  ID_TEMPER
  visco_ker % llaws(7) % where     = 'IELEM'
  visco_ker % llaws(7) % kfl_gradi = 1
  visco_ker % llaws(7) % kfl_deriv = 1
  visco_ker % llaws(7) % kfl_grder = 1

  visco_ker % llaws(8) % wname     = 'ABL  '       ! ABL at Gauss points: mu = rho*kap/Cmu*(z+z0)*U_*
  visco_ker % llaws(8) % lresp(1)  =  ID_NASTIN
  visco_ker % llaws(8) % where     = 'IELEM'
  visco_ker % llaws(8) % kfl_gradi = 1

  visco_ker % llaws(9) % wname     = 'MUTAB'       ! from table (CFI combustion model)
  visco_ker % llaws(9) % lresp(1)  =  ID_CHEMIC
  visco_ker % llaws(9) % where     = 'IPOIN'
  visco_ker % llaws(9) % kfl_gradi = 1

  visco_ker % llaws(10) % wname     = 'TMUTA'      ! from table (CFI combustion model)
  visco_ker % llaws(10) % lresp(1)  =  ID_TEMPER   ! update only in temper
  visco_ker % llaws(10) % where     = 'IPOIN'
  visco_ker % llaws(10) % kfl_gradi = 1

  call ker_polynomial_allaws(visco_ker, 11_ip, ID_TEMPER, 'IELEM', 0_ip)

  visco_ker % llaws(12) % wname     = 'GAUSS'      ! Gaussian filter for generating turbulent fluctuations 
  visco_ker % llaws(12) % lresp(1)  =  ID_NASTIN   !  
  visco_ker % llaws(12) % where     = 'IELEM'
  visco_ker % llaws(12) % kfl_gradi = 1

  !----------------------------------------------------------------------
  !
  ! Conductivity
  !
  !----------------------------------------------------------------------

  condu_ker % llaws(2) % wname    = 'KMIXT'        ! Mixture of species
  condu_ker % llaws(2) % lresp(1) =  ID_CHEMIC 
  condu_ker % llaws(2) % where     = 'IPOIN'
  condu_ker % llaws(2) % kfl_gradi = 1

  condu_ker % llaws(3) % wname    = 'SUTHE'        ! Sutherland
  condu_ker % llaws(3) % lresp(1) =  ID_TEMPER
  condu_ker % llaws(3) % where     = 'IPOIN'
  condu_ker % llaws(3) % kfl_gradi = 1

  condu_ker % llaws(4) % wname     = 'GPSUT'       ! Sutherland in Gauss points
  condu_ker % llaws(4) % lresp(1)  =  ID_NASTAL
  condu_ker % llaws(4) % lresp(2)  =  ID_TEMPER
  condu_ker % llaws(4) % where     = 'IELEM'
  condu_ker % llaws(4) % kfl_gradi =  1

  condu_ker % llaws(5) % wname    = 'KTABL'        ! from table (CFI combustion model)
  condu_ker % llaws(5) % lresp(1) =  ID_CHEMIC
  condu_ker % llaws(5) % where     = 'IPOIN'
  condu_ker % llaws(5) % kfl_gradi = 1

  condu_ker % llaws(6) % wname    = 'TKTAB'        ! from table (CFI combustion model)
  condu_ker % llaws(6) % lresp(1) =  ID_TEMPER     ! update only in temper
  condu_ker % llaws(6) % where     = 'IPOIN'
  condu_ker % llaws(6) % kfl_gradi = 1

  condu_ker % llaws(7) % wname    = 'KQUAR'        ! quartz glass, from polynomial function, for CHT
  condu_ker % llaws(7) % lresp(1) =  ID_TEMPER     ! update only in temper
  condu_ker % llaws(7) % where     = 'IPOIN'
  condu_ker % llaws(7) % kfl_gradi = 1

  condu_ker % llaws(8) % wname    = 'KSTAI'        ! stainless steel, from polynomial function, for CHT
  condu_ker % llaws(8) % lresp(1) =  ID_TEMPER     ! update only in temper
  condu_ker % llaws(8) % where     = 'IPOIN'
  condu_ker % llaws(8) % kfl_gradi = 1

  call ker_polynomial_allaws(condu_ker, 9_ip, ID_TEMPER, 'IELEM', 0_ip)

  !----------------------------------------------------------------------
  !
  ! Specific Heat
  !
  !----------------------------------------------------------------------

  sphea_ker % llaws(2) % wname    = 'CPMIX'        ! Mixture of species
  sphea_ker % llaws(2) % lresp(1) =  ID_CHEMIC 
  sphea_ker % llaws(2) % where    = 'IPOIN'

  sphea_ker % llaws(3) % wname    = 'CPTAB'        ! from table (CFI combustion model)
  sphea_ker % llaws(3) % lresp(2) =  ID_TEMPER
  sphea_ker % llaws(3) % where    = 'IPOIN'

  sphea_ker % llaws(4) % wname    = 'TCPTA'        ! from table (CFI combustion model)
  sphea_ker % llaws(4) % lresp(1) =  ID_TEMPER     ! update only in temper
  sphea_ker % llaws(4) % where    = 'IPOIN'

  sphea_ker % llaws(5) % wname    = 'CPQUA'        ! quartz glass, from polynomial function, for CHT
  sphea_ker % llaws(5) % lresp(1) =  ID_TEMPER     ! update only in temper
  sphea_ker % llaws(5) % where    = 'IPOIN'

  sphea_ker % llaws(6) % wname    = 'CPSTA'        ! stainless steel, from polynomial function, for CHT 
  sphea_ker % llaws(6) % lresp(1) =  ID_TEMPER     ! update only in temper
  sphea_ker % llaws(6) % where    = 'IPOIN'

  call ker_polynomial_allaws(sphea_ker, 7_ip, ID_TEMPER, 'IELEM', 0_ip) ! 

  !----------------------------------------------------------------------
  !
  ! Dummy variable can be used for whatever
  !
  !----------------------------------------------------------------------

  dummy_ker % llaws(2) % wname     = 'BIFL2'       ! Bifluid
  dummy_ker % llaws(2) % lresp(1)  =  ID_LEVELS
  dummy_ker % llaws(2) % where     = 'IELEM'
  dummy_ker % llaws(2) % kfl_gradi = 1  ! BIFL2 includes gradient

  !----------------------------------------------------------------------
  !
  ! Turbulent viscosity
  !
  !----------------------------------------------------------------------

  turmu_ker % llaws( 2) % wname     = 'SMAGO'       ! Turbulent viscosity mu_t from Smagorinsky
  turmu_ker % llaws( 2) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws( 2) % where     = 'IELEM'
  turmu_ker % llaws( 2) % kfl_gradi = 1

  turmu_ker % llaws( 3) % wname     = 'WALE '       ! Turbulent viscosity mu_t from WALE
  turmu_ker % llaws( 3) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws( 3) % where     = 'IELEM'
  turmu_ker % llaws( 3) % kfl_gradi = 1 

  turmu_ker % llaws( 4) % wname     = 'SSTKO'       ! Turbulent viscosity mu_t from RANS K-W SST Model
  turmu_ker % llaws( 4) % lresp(1)  =  ID_TURBUL
  turmu_ker % llaws( 4) % where     = 'IELEM'
  turmu_ker % llaws( 4) % kfl_gradi = 1
  turmu_ker % llaws( 4) % kfl_deriv_tur = 1
  turmu_ker % llaws( 4) % kfl_deriv_vel = 1

  turmu_ker % llaws( 5) % wname     = 'STDKE'       ! Turbulent viscosity mu_t from RANS K-EPS Model
  turmu_ker % llaws( 5) % lresp(1)  =  ID_TURBUL
  turmu_ker % llaws( 5) % lresp(2)  =  ID_NASTIN    ! for k e fp and realizable models
  turmu_ker % llaws( 5) % where     = 'IELEM'
  turmu_ker % llaws( 5) % kfl_gradi = 1

  turmu_ker % llaws( 6) % wname     = 'SPALA'       ! Turbulent viscosity mu_t from RANS Spalart-Almaras Model
  turmu_ker % llaws( 6) % lresp(1)  =  ID_TURBUL
  turmu_ker % llaws( 6) % where     = 'IELEM'
  turmu_ker % llaws( 6) % kfl_deriv_tur = 1

  turmu_ker % llaws( 7) % wname     = 'VRMAN'       ! Turbulent viscosity mu_t from Vreman
  turmu_ker % llaws( 7) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws( 7) % where     = 'IELEM'
  turmu_ker % llaws( 7) % kfl_gradi = 1

  turmu_ker % llaws( 8) % wname     = 'ILSA'       ! Turbulent viscosity mu_t from ILSA
  turmu_ker % llaws( 8) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws( 8) % where     = 'IELEM'
  turmu_ker % llaws( 8) % kfl_gradi = 1

  turmu_ker % llaws( 9) % wname     = 'MIXIN'       ! Turbulent viscosity mu_t from algebraic mixing-length model
  turmu_ker % llaws( 9) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws( 9) % where     = 'IELEM'
  turmu_ker % llaws( 9) % kfl_gradi = 1

  turmu_ker % llaws(10) % wname     = 'SIGMA'       ! Turbulent viscosity mu_t from Vreman
  turmu_ker % llaws(10) % lresp(1)  =  ID_NASTIN
  turmu_ker % llaws(10) % where     = 'IELEM'
  turmu_ker % llaws(10) % kfl_gradi = 1

  !----------------------------------------------------------------------
  !
  ! Canopy porosity
  !
  !----------------------------------------------------------------------

  poros_ker % llaws(2) % wname     = 'CANOP'
  poros_ker % llaws(2) % lresp(1)  =  ID_NASTIN
  poros_ker % llaws(2) % where     = 'IELEM'
  poros_ker % llaws(2) % kfl_gradi = 0

  !----------------------------------------------------------------------
  !
  ! Darcy-Forchheimer porosity
  !
  !----------------------------------------------------------------------

  poros_ker % llaws(3) % wname     = 'DAFOR'
  poros_ker % llaws(3) % lresp(1)  =  ID_NASTIN
  poros_ker % llaws(3) % where     = 'IELEM'
  poros_ker % llaws(3) % kfl_gradi = 0

    !----------------------------------------------------------------------
  !
  ! Mixing function
  !
  !----------------------------------------------------------------------

  mixin_ker % llaws(2) % wname     = 'IMMER'       ! Bifluid
  mixin_ker % llaws(2) % lresp(1)  =  ID_LEVELS
  mixin_ker % llaws(2) % where     = 'IELEM'
  mixin_ker % llaws(2) % kfl_gradi = 0             ! Despite there is a gradient we never use it for Bifulid flow - see if in nsi_elmop3 grvis

end subroutine ker_allaws
