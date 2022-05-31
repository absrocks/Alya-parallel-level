module mod_chm_sectional_soot_model

  use def_kintyp, only   : ip,rp,r3p
  use def_domain, only   : ndime,npoin,coord,mgaus,mnode,nelem
  use def_master, only   : inotmaster
  use mod_memory, only   : memory_alloca, memory_deallo
  use mod_memory, only   : memory_deallo
  use def_master, only   : mem_modul,modul,momod
  use def_chemic, only   : nclas_chm
  use def_chemic, only   : nspec_chm

  implicit none

  !
  ! Soot constants
  !
  integer(ip), parameter                 ::       &
       NCVmin = 32                                  ! Number of C atoms in Vmin

  real(rp), parameter                    ::       &
       pi = 3.1415926535897932384626433832795_rp, & ! pi
       MC = 12.0112_rp,                           & ! Atomic weight of Carbon (gm/mol)
       NAvo = 6.02214129e23_rp,                   & ! Avogadros number (1/mol)
       kBoltz = 1.38064852e-16_rp                   ! Boltzmann constsnt (cm2 g s-2 K-1)

  integer(ip)                           :: &
       indexS(5),                          & ! Index of soot sources models 
       gasCoupling_ssm,                    & ! Coupling with gas phase
       nPAH,                               & ! Number of species involved in nucleation = condensation 
       nSurf,                              & ! Number of species of the problem 
       idNucl(10),                         & ! Index species involved in nucleation 
       idCond(10),                         & ! Index species involved in condensation
       idSurf(10),                         & ! Index species involved in surface growth
       ID_NUCL,                            & ! Index nucleation process
       ID_COND,                            & ! Index condensation process
       ID_COAG,                            & ! Index coagulation process
       ID_SURF,                            & ! Index surface growth process
       nspec_ssm,                          & ! Number of species involved in surface growth
       nclas_ssm,                          & ! Total number of unknowns: Y_tot = Y_k + Y_s 
       nsect_ssm                             ! Number of sections for soot sectional model

  !
  ! Soot model variables
  !
  type(r3p),  pointer                   :: &
       Qsoot_gp_chm(:),                    & ! Soot volume fraction
       QsootDist_gp_chm(:),                & ! Soot volume fraction density
       NDsoot_gp_chm(:),                   & ! Number density
       Asoot_gp_chm(:),                    & ! Soot surface area
       Vsoot_gp_chm(:)                       ! Soot volume

  !
  ! Soot source terms
  !
  type(r3p),  pointer                   :: &
       Qnucl_gp_chm(:),                    & ! Nucleation 
       Qcoag_gp_chm(:),                    & ! Coagulation 
       Qcond_gp_chm(:),                    & ! Condensation
       Qsurf_gp_chm(:),                    & ! Surface growth
       Qtot_gp_chm(:)                        ! Total

  !
  ! Species source terms
  !
  type(r3p),  pointer                   :: &
       SSTnucl_gp_chm(:),                  & ! Nucleation
       SSTcond_gp_chm(:),                  & ! Condensation
       SSTsurf_gp_chm(:),                  & ! Surface growth and oxidation
       SSTtot_gp_chm(:)                      ! Total

  !
  ! Sectional soot model parameters
  !
  real(rp)                              :: &
       SootDensity_ssm,                    & ! Activation soot model
       RhoC_ssm,                           & ! Density of Carbon (gm/cm3)
       Vmax_ssm,                           & ! Index of soot sources models 
       RadSoot_ssm                           ! Radiation parameter

  !
  ! Sectional variables
  !
  real(rp)                              :: &
       Vsec(100),                          & ! Sectional volume discretization
       Vpmean(100),                        & ! Mean particle volume in section
       Dpmean(100)                           ! Mean particle diameter
  !
  ! Unit conversion
  !
  real(rp), parameter                   :: &
       UnC = 1000.0_rp                       ! Unit conversion factor 1: CGS InOut, 1000: MKS InOut

  integer(ip)                           :: &
       id_A4,                              & ! Index of A4
       id_H,                               & ! Index of H
       id_H2,                              & ! Index of H2
       id_OH,                              & ! Index of OH
       id_O2,                              & ! Index of O2
       id_H2O,                             & ! Index of H2O
       id_C2H2,                            & ! Index of C2H2
       id_CO                                 ! Index of CO

  private

  !
  ! Public variables
  !
  public :: SootDensity_ssm
  public :: Vmax_ssm
  public :: RadSoot_ssm
  public :: nsect_ssm
  public :: nspec_ssm
  public :: nclas_ssm
  public :: gasCoupling_ssm
  public :: RhoC_ssm
  public :: indexS,nPAH,nSurf 
  public :: idNucl,idCond,idSurf 
  public :: ID_NUCL,ID_COND,ID_COAG,ID_SURF 

  !
  ! Public subroutines 
  !
  public :: chm_assembly_soot_sources_ssm
  !!public :: chm_readingIndex_ssm
  !!public :: chm_readChem1d_ssm
  public :: chm_initialization_ssm
  !!public :: chm_calculate_soot_sources_ssm
  !!public :: chm_writeFile_ssm
  !!public :: chm_source_nucleation_ssm
  !!public :: chm_source_condensation_ssm
  !!public :: chm_source_coagulation_ssm
  !!public :: chm_source_surfgrowth_ssm
  public :: chm_validation_test_ssm
!!DMM  public :: chm_update_variables_ssm
  public :: chm_sectional_soot_memory_ssm
  public :: soot_volume_discretization_ssm

  contains

    subroutine chm_validation_test_ssm()

      use def_domain,   only : ltype,ngaus,mgaus
      use def_chemic,   only : kfl_soot_chm
      use mod_messages, only : messages_live 

      implicit none

      integer(ip), parameter :: nspec = 202   ! Input from chem1d: number of species
      integer(ip), parameter :: npts  = 200   ! Input from chem1d: number of points

      integer(ip)            :: ielem,pelty
      real(rp)               :: gpcon(mgaus,nclas_chm)
      real(rp)               :: gptem(mgaus)
      real(rp)               :: gpden(mgaus)
      real(rp)               :: gpsoot(mgaus,nsect_ssm)

      !
      ! Start soot sources calculation
      !
      call header()

      !
      ! Reading index input file
      !
      call chm_readingIndex_ssm()

      !
      ! Initialization soot source terms
      !
      call chm_initialization_ssm()

      !
      ! Reading chem1d input file
      !
      call chm_readChem1d_ssm(nspec_chm,npts,gpden,gpcon,gptem,gpsoot)

      !
      ! Compute Volume discretization
      !    Vsec,Vpmean,Dpmean
      !
      call soot_volume_discretization_ssm()

      !
      ! Converting soot mass fractions to soot sectional variables
      !
      call chm_update_variables()

      !
      ! Compute soot source terms
      !
      call chm_calculate_soot_sources_ssm(1_ip,npts,gpcon,gptem,gpden,gpsoot)

      !
      ! Writing output file 
      !
      call chm_writeFile_ssm()

      !
      ! End soot source terms calculations 
      !
      call footer()

      !
      ! Stop the calculation 
      !
      call messages_live('SECTONAL METHOD VALIDATION TEST FINISHED.')           
      call messages_live('INITIAL SOLUTION','END SECTION')           
      call runend('O.K.!')

    end subroutine chm_validation_test_ssm

    subroutine chm_assembly_soot_sources_ssm(ielem,pgaus,gpcon,gptem,gpden,gprhs,gpsoot)
      !-----------------------------------------------------------------------
      !****f* Chemic/mod_chm_element_operations/chm_turbul
      ! NAME 
      !    chm_turbul_finiteRate:w

      ! DESCRIPTION
      !    Add subgrid contribution to diffusion term
      ! USES
      ! USED BY
      !    chm_turbul_finiteRate
      !***
      !-----------------------------------------------------------------------

      use def_domain, only : ltype,ngaus
      implicit none

      integer(ip),  intent(in)     :: ielem
      integer(ip),  intent(in)     :: pgaus
      real(rp),     intent(in)     :: gpcon(pgaus,nclas_ssm,*)
      real(rp),     intent(in)     :: gptem(pgaus)
      real(rp),     intent(in)     :: gpden(pgaus)
      real(rp),     intent(in)     :: gpsoot(pgaus,nsect_ssm) ! Mass fraction soot species Ys
      real(rp),     intent(inout)  :: gprhs(pgaus,nclas_ssm)

      integer(ip)                  :: isect,iclas
      !
      ! Calculate soot source terms
      !
      call chm_calculate_soot_sources_ssm(ielem,pgaus,gpcon,gptem,gpden,gpsoot)

      !
      ! Impose RHS sectional soot terms
      !
      do isect=1,nsect_ssm
         iclas = isect + nspec_ssm
         gprhs(1:pgaus,iclas) = Qtot_gp_chm(ielem) % a (1:pgaus,isect,1)
      end do

      !
      ! Impose RHS species soot terms
      !
      do iclas=1,nspec_ssm

      end do 


    end subroutine chm_assembly_soot_sources_ssm

    subroutine chm_sectional_soot_memory_ssm()
      !-----------------------------------------------------------------------
      !****f* Chemic/mod_chm_element_operations/chm_turbul
      ! NAME 
      !    chm_turbul_finiteRate:w

      ! DESCRIPTION
      !    Add subgrid contribution to diffusion term
      ! USES
      ! USED BY
      !    chm_turbul_finiteRate
      !***
      !-----------------------------------------------------------------------

      use def_domain, only : ltype,ngaus
      implicit none

      integer(ip)     :: ielem,pgaus,pelty


         !
         ! Soot variables 
         !
         call memory_alloca(mem_modul(1:2,modul),'QSOOT_GP_CHM','chm_solmem',    Qsoot_gp_chm,nelem)
         call memory_alloca(mem_modul(1:2,modul),'QSOOTDIST_GP_CHM','chm_solmem',QsootDist_gp_chm,nelem)
         call memory_alloca(mem_modul(1:2,modul),'NDSOOT_GP_CHM','chm_solmem',   NDsoot_gp_chm,nelem)
         call memory_alloca(mem_modul(1:2,modul),'ASOOT_GP_CHM','chm_solmem',    Asoot_gp_chm,nelem)
         call memory_alloca(mem_modul(1:2,modul),'VSOOT_GP_CHM','chm_solmem',    Vsoot_gp_chm,nelem)

         !
         ! Soot source terms
         !
         call memory_alloca(mem_modul(1:2,modul),'QNUCL_GP_CHM','chm_solmem',Qnucl_gp_chm,nelem)
         call memory_alloca(mem_modul(1:2,modul),'QCOAG_GP_CHM','chm_solmem',Qcoag_gp_chm,nelem)
         call memory_alloca(mem_modul(1:2,modul),'QCOND_GP_CHM','chm_solmem',Qcond_gp_chm,nelem)
         call memory_alloca(mem_modul(1:2,modul),'QSURF_GP_CHM','chm_solmem',Qsurf_gp_chm,nelem)
         call memory_alloca(mem_modul(1:2,modul),'QTOT_GP_CHM' ,'chm_solmem', Qtot_gp_chm,nelem)

         !
         ! Species source terms
         !
         call memory_alloca(mem_modul(1:2,modul),'SSTNUCL_GP_CHM','chm_solmem',SSTnucl_gp_chm,nelem)
         call memory_alloca(mem_modul(1:2,modul),'SSTCOND_GP_CHM','chm_solmem',SSTcond_gp_chm,nelem)
         call memory_alloca(mem_modul(1:2,modul),'SSTSURF_GP_CHM','chm_solmem',SSTsurf_gp_chm,nelem)
         call memory_alloca(mem_modul(1:2,modul),'SSTTOT_GP_CHM','chm_solmem' , SSTtot_gp_chm,nelem)

         do ielem = 1,nelem 
            pelty = ltype(ielem)  
            pgaus = ngaus(pelty)      

            !
            ! Soot model variables
            !
            call memory_alloca(mem_modul(1:2,modul),'QSOOT_GP_CHM','chm_solmem',Qsoot_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
            call memory_alloca(mem_modul(1:2,modul),'QSOOTDIST_GP_CHM','chm_solmem',QsootDist_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
            call memory_alloca(mem_modul(1:2,modul),'NDSOOT_GP_CHM','chm_solmem',NDsoot_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
            call memory_alloca(mem_modul(1:2,modul),'ASOOT_GP_CHM','chm_solmem',Asoot_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
            call memory_alloca(mem_modul(1:2,modul),'VSOOT_GP_CHM','chm_solmem',Vsoot_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)

            !
            ! Soot source terms
            !
            call memory_alloca(mem_modul(1:2,modul),'QNUCL_GP_CHM','chm_solmem',Qnucl_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
            call memory_alloca(mem_modul(1:2,modul),'QCOAG_GP_CHM','chm_solmem',Qcoag_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
            call memory_alloca(mem_modul(1:2,modul),'QCOND_GP_CHM','chm_solmem',Qcond_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
            call memory_alloca(mem_modul(1:2,modul),'QSURF_GP_CHM','chm_solmem',Qsurf_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)
            call memory_alloca(mem_modul(1:2,modul),'QTOT_GP_CHM' ,'chm_solmem', Qtot_gp_chm(ielem)%a,pgaus,nsect_ssm,1_ip)

            !
            ! Species source terms
            !
            call memory_alloca(mem_modul(1:2,modul),'SSTNUCL_GP_CHM','chm_solmem',SSTnucl_gp_chm(ielem)%a,pgaus,nspec_ssm,1_ip)
            call memory_alloca(mem_modul(1:2,modul),'SSTCOND_GP_CHM','chm_solmem',SSTcond_gp_chm(ielem)%a,pgaus,nspec_ssm,1_ip)
            call memory_alloca(mem_modul(1:2,modul),'SSTSURF_GP_CHM','chm_solmem',SSTsurf_gp_chm(ielem)%a,pgaus,nspec_ssm,1_ip)
            call memory_alloca(mem_modul(1:2,modul),'SSTTOT_GP_CHM','chm_solmem' ,SSTtot_gp_chm(ielem)%a,pgaus,nspec_ssm,1_ip)

         end do

      !!DMM else
      !!DMM    !
      !!DMM    ! Soot variables 
      !!DMM    !
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QSOOT_GP_CHM','chm_solmem',    Qsoot_gp_chm,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QSOOTDIST_GP_CHM','chm_solmem',QsootDist_gp_chm,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'NDSOOT_GP_CHM','chm_solmem',   NDsoot_gp_chm,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'ASOOT_GP_CHM','chm_solmem',    Asoot_gp_chm,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'VSOOT_GP_CHM','chm_solmem',    Vsoot_gp_chm,1_ip)

      !!DMM    !
      !!DMM    ! Soot source terms
      !!DMM    !
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QNUCL_GP_CHM','chm_solmem',Qnucl_gp_chm,1_ip)  
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QCOAG_GP_CHM','chm_solmem',Qcoag_gp_chm,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QCOND_GP_CHM','chm_solmem',Qcond_gp_chm,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QSURF_GP_CHM','chm_solmem',Qsurf_gp_chm,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QTOT_GP_CHM' ,'chm_solmem', Qtot_gp_chm,1_ip)

      !!DMM    !
      !!DMM    ! Species source terms
      !!DMM    !
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'SSTNUCL_GP_CHM','chm_solmem',SSTnucl_gp_chm,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'SSTCOND_GP_CHM','chm_solmem',SSTcond_gp_chm,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'SSTSURF_GP_CHM','chm_solmem',SSTsurf_gp_chm,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'SSTTOT_GP_CHM','chm_solmem' , SSTtot_gp_chm,1_ip)

      !!DMM    !
      !!DMM    ! Soot model variables
      !!DMM    !
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QSOOT_GP_CHM','chm_solmem',Qsoot_gp_chm(ielem)%a,1_ip,1_ip,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QSOOTDIST_GP_CHM','chm_solmem',QsootDist_gp_chm(ielem)%a,1_ip,1_ip,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'NDSOOT_GP_CHM','chm_solmem',NDsoot_gp_chm(ielem)%a,1_ip,1_ip,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'ASOOT_GP_CHM','chm_solmem',Asoot_gp_chm(ielem)%a,1_ip,1_ip,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'VSOOT_GP_CHM','chm_solmem',Vsoot_gp_chm(ielem)%a,1_ip,1_ip,1_ip)

      !!DMM    !
      !!DMM    ! Soot source terms
      !!DMM    !
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QNUCL_GP_CHM','chm_solmem',Qnucl_gp_chm(ielem)%a,1_ip,1_ip,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QCOAG_GP_CHM','chm_solmem',Qcoag_gp_chm(ielem)%a,1_ip,1_ip,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QCOND_GP_CHM','chm_solmem',Qcond_gp_chm(ielem)%a,1_ip,1_ip,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QSURF_GP_CHM','chm_solmem',Qsurf_gp_chm(ielem)%a,1_ip,1_ip,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'QTOT_GP_CHM' ,'chm_solmem', Qtot_gp_chm(ielem)%a,1_ip,1_ip,1_ip)

      !!DMM    !
      !!DMM    ! Species source terms
      !!DMM    !
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'SSTNUCL_GP_CHM','chm_solmem',SSTnucl_gp_chm(ielem)%a,1_ip,1_ip,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'SSTCOND_GP_CHM','chm_solmem',SSTcond_gp_chm(ielem)%a,1_ip,1_ip,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'SSTSURF_GP_CHM','chm_solmem',SSTsurf_gp_chm(ielem)%a,1_ip,1_ip,1_ip)
      !!DMM    call memory_alloca(mem_modul(1:2,modul),'SSTTOT_GP_CHM','chm_solmem' ,SSTtot_gp_chm(ielem)%a, 1_ip,1_ip,1_ip)

      !!DMM endif


    end subroutine chm_sectional_soot_memory_ssm


    subroutine chm_update_variables()
      !-----------------------------------------------------------------------
      !****f* Chemic/mod_chm_element_operations/chm_turbul
      ! NAME 
      !    chm_turbul_finiteRate:w

      ! DESCRIPTION
      !    Add subgrid contribution to diffusion term
      ! USES
      ! USED BY
      !    chm_turbul_finiteRate
      !***
      !-----------------------------------------------------------------------

      use def_domain, only   : nelem,ltype,ngaus

      implicit none

      integer(ip)     :: ielem,igaus
      integer(ip)     :: pgaus,pelty

      !
      ! Compute sectional soot model variables from Chem1D and PSD from Sectional Method
      !
      do ielem = 1,nelem
         pelty = ltype(ielem)
         pgaus = ngaus(pelty)
         do igaus=1,pgaus
            Qsoot_gp_chm(ielem)%a(igaus,1:nsect_ssm,1_ip)     = 0.0_rp
            QsootDist_gp_chm(ielem)%a(igaus,1:nsect_ssm,1_ip) = 0.0_rp
            NDsoot_gp_chm(ielem)%a(igaus,1:nsect_ssm,1_ip)    = 0.0_rp
            Asoot_gp_chm(ielem)%a(igaus,1:nsect_ssm,1_ip)     = 0.0_rp
            Vsoot_gp_chm(ielem)%a(igaus,1:nsect_ssm,1_ip)     = 0.0_rp
         end do
      end do

    end subroutine chm_update_variables


    subroutine chm_calculate_soot_sources_ssm(ielem,pgaus,gpcon,gptem,gpden,gpsoot)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !    
      ! DESCRIPTION
      !    
      ! USES
      ! USED BY
      !    
      !***
      !-----------------------------------------------------------------------
      use def_domain, only   : nelem,ltype,ngaus
      use def_chemic, only   : kfl_soot_chm

      implicit none
      integer(ip),  intent(in)  :: ielem
      integer(ip),  intent(in)  :: pgaus
      real(rp),     intent(in)  :: gpcon(pgaus,nspec_ssm)  ! Mass fraction species Yk
      real(rp),     intent(in)  :: gptem(pgaus)            ! Temperature
      real(rp),     intent(in)  :: gpden(pgaus)            ! Density
      real(rp),     intent(in)  :: gpsoot(pgaus,nsect_ssm) ! Mass fraction soot species Ys

      !
      ! Compute sectional soot model variables
      !
      Qsoot_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip)     = soot_vol_frac( pgaus,gpden(1:pgaus), &
                                                                           gpsoot(1:pgaus,1:nsect_ssm) )

      QsootDist_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip) = soot_vol_frac_dist( pgaus, &
                                                                          Qsoot_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip) )

      NDsoot_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip)    = soot_num_den( pgaus, &
                                                                          QsootDist_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip) )

      Asoot_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip)     = soot_surf_area( pgaus, &
                                                                            QsootDist_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip) )

      Vsoot_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip)     = soot_vol( pgaus, &
                                                                      QsootDist_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip), &
                                                                      NDsoot_gp_chm(ielem)%a(1:pgaus,1:nsect_ssm,1_ip) )

      !
      ! Compute source term nucleation
      !
      call chm_source_nucleation_ssm(   pgaus,gpcon(1:pgaus,1:nspec_ssm),gptem(1:pgaus),gpden(1:pgaus), &
                                        Qnucl_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1), &
                                        SSTnucl_gp_chm(ielem) % a(1:pgaus,1:nspec_ssm,1) )

      !
      ! Compute source term coagulation
      !
      call chm_source_coagulation_ssm(  pgaus,gpcon(1:pgaus,1:nspec_ssm),gptem(1:pgaus),gpden(1:pgaus), & 
                                        NDsoot_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1) , &
                                        Qcoag_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1) ) 

      !
      ! Compute source term condensation
      !
      call chm_source_condensation_ssm( pgaus,gpcon(1:pgaus,1:nspec_ssm),gptem(1:pgaus),gpden(1:pgaus), &
                                        NDsoot_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1) , &
                                        Qcond_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1), &
                                        SSTcond_gp_chm(ielem) % a(1:pgaus,1:nspec_ssm,1) )

      !
      ! Compute source term surface growth
      !
      call chm_source_surfgrowth_ssm(   pgaus,gpcon(1:pgaus,1:nspec_ssm),gptem(1:pgaus),gpden(1:pgaus), &
                                        NDsoot_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1), &
                                        Asoot_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1), &
                                        Vsoot_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1), &
                                        Qsurf_gp_chm(ielem) % a(1:pgaus,1:nsect_ssm,1), &
                                        SSTsurf_gp_chm(ielem) % a(1:pgaus,1:nspec_ssm,1) )

      !
      ! Compute total source terms
      !
      Qtot_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) = UnC * RhoC_ssm * ( Qnucl_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) &
                                                  + Qcond_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) &
                                                  + Qsurf_gp_chm(ielem) % a (1:pgaus,1:nsect_ssm, 1) )

      SSTtot_gp_chm(ielem) % a (1:pgaus,1:nspec_ssm, 1) = UnC * ( SSTnucl_gp_chm(ielem) % a (1:pgaus,1:nspec_ssm, 1) &
                                                        + SSTcond_gp_chm(ielem) % a (1:pgaus,1:nspec_ssm, 1) &
                                                        + SSTsurf_gp_chm(ielem) % a (1:pgaus,1:nspec_ssm, 1) )

    end subroutine chm_calculate_soot_sources_ssm

    subroutine chm_initialization_ssm()
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !    
      ! DESCRIPTION
      !    
      ! USES
      ! USED BY
      !    
      !***
      !-----------------------------------------------------------------------

      implicit none
      integer(ip)     :: ielem                ! Elements of the mesh

      !
      ! Initialization soot and species source terms
      !
      do ielem = 1,nelem
         Qnucl_gp_chm(ielem) % a = 0.0_rp 
         Qcoag_gp_chm(ielem) % a = 0.0_rp
         Qcond_gp_chm(ielem) % a = 0.0_rp
         Qsurf_gp_chm(ielem) % a = 0.0_rp
         Qtot_gp_chm(ielem)  % a = 0.0_rp

         SSTnucl_gp_chm(ielem) % a = 0.0_rp
         SSTcond_gp_chm(ielem) % a = 0.0_rp 
         SSTsurf_gp_chm(ielem) % a = 0.0_rp
         SSTtot_gp_chm(ielem)  % a = 0.0_rp
      end do

    end subroutine chm_initialization_ssm

    subroutine chm_readingIndex_ssm()
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !    
      ! DESCRIPTION
      !    
      ! USES
      ! USED BY
      !    
      !***
      !-----------------------------------------------------------------------

      implicit none
      !
      ! To be read by input data and define globally
      !
      nPAH  = 2
      nSurf = 7

      idNucl(1:nPAH)  = (/138,7/)
      idCond(1:nPAH)  = (/138,7/)
      idsurf(1:nSurf) = (/3,7,5,10,8,26,22/)

      !
      ! Alternative that reads the species indexes from mechanism !!AJK
      !
      call soot_read_species()

      !
      ! Index soot physical process ==> move to reaphy 
      !
      ID_NUCL = 1
      ID_COND = 2
      ID_COAG = 3
      ID_SURF = 4

      !
      ! Soot source term model activation
      !
      indexS(ID_NUCL) = 1   ! Nucleation
      indexS(ID_COND) = 1   ! Condensation
      indexS(ID_COAG) = 1   ! Coagulation
      indexS(ID_SURF) = 1   ! Surface growth

    end subroutine chm_readingIndex_ssm

    subroutine chm_readChem1d_ssm(nspec, npts, Rho, Yspec, Temp, Ysoot)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !    
      ! DESCRIPTION
      !    
      ! USES
      ! USED BY
      !    
      !***
      !-----------------------------------------------------------------------
      use def_chemic ,    only : W_k, kfl_soot_chm
      use def_master,     only : intost
      use iso_fortran_env,only : iostat_eor,iostat_end
      implicit none

      integer(ip),  intent(in)            :: nspec
      integer(ip),  intent(in)            :: npts
      integer(ip)                         :: i,j,k
      real(rp), dimension(npts)           :: X
      real(rp), dimension(npts)           :: Temp
      real(rp), dimension(npts)           :: Rho
      real(rp), dimension(npts,nspec)     :: Yspec
      real(rp), dimension(npts,nsect_ssm) :: Ysoot
      real(rp), dimension(nspec)          :: Wk
      character(len=16)                   :: MechName
      integer(ip)                         :: ioerror

      !
      ! Read the variables from chem1d
      !
      open (111, file = 'yiend_chem1d.dat', action = 'read', status = 'old')
      do i = 1,npts
          read(111,*, iostat=ioerror) X(i),Temp(i),Rho(i),Yspec(i,:),Ysoot(i,:)
          if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
          if (ioerror /= 0 ) call runend('chm_readChem1d_ssm: cannot read chem1d solution. error code:'//trim(intost(ioerror)))
      end do
      close(111)

      print*,'    Inputs read !'
 
    end subroutine chm_readChem1d_ssm

    subroutine chm_writeFile_ssm
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !    
      ! DESCRIPTION
      !    
      ! USES
      ! USED BY
      !    
      !***
      !-----------------------------------------------------------------------
      implicit none
 
    end subroutine chm_writeFile_ssm

    subroutine chm_source_nucleation_ssm(pgaus,gpcon,gptem,gpden,Qnucl,SSTnucl)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !    
      ! DESCRIPTION
      !    
      ! USES
      ! USED BY
      !    
      !***
      !-----------------------------------------------------------------------
      use def_chemic , only: W_k, kfl_soot_chm
      implicit none
 
      integer(ip),  intent(in)  :: pgaus                    ! Vector size
      real(rp),     intent(in)  :: gpcon(pgaus,nspec_chm)   ! Mass fraction species Yk
      real(rp),     intent(in)  :: gptem(pgaus)             ! Temperature
      real(rp),     intent(in)  :: gpden(pgaus)             ! Temperature
      real(rp),     intent(out) :: Qnucl(pgaus,nsect_ssm)   ! Source term nucleation
      real(rp),     intent(out) :: SSTnucl(pgaus,nspec_chm) ! Species source term from nucleation

      integer(ip)               :: k, isect
      real(rp)                  :: gpXk(pgaus,nspec_chm)    ! Molar concentration species Ck (mol/cm3s)
      real(rp)                  :: Betaij_Nuc(pgaus)        ! Collison frequancy factor
      real(rp)                  :: NDPAH(pgaus)             ! Number density of PAH
      real(rp)                  :: VPAH                     ! Volume of PAH

      Qnucl   = 0.0_rp
      SSTnucl = 0.0_rp

      gpXk = spec_con(pgaus,gpden,gpcon)

      NDPAH = max(gpXk(:,id_A4)*NAvo, 0.0_rp)
      VPAH  = 0.5_rp * NCVmin * (MC/(NAvo*RhoC_ssm))

      ! 
      ! Compute collision frequencies for nucleation 
      ! 
      call soot_collision_freq(pgaus,VPAH,VPAH,gptem,Betaij_Nuc)

      isect = 1
      do k = 1,pgaus
         Qnucl(k,isect) = 2.0_rp*VPAH*Betaij_Nuc(k)*(NDPAH(k)**2.0_rp)
      end do

      ! 
      ! Coupling to gas phase (gm/cm3)
      ! 
      do k = 1,pgaus
         SSTnucl(k,id_A4)  = -1.0_rp*(Qnucl(k,isect)/(VPAH*NAvo))*(W_k(id_A4)*UnC)
         SSTnucl(k,id_H2)  =  5.0_rp*(Qnucl(k,isect)/(VPAH*NAvo))*(W_k(id_H2)*UnC)
      end do

      !
      ! Writing nucleation source
      !
      if (kfl_soot_chm .eq. -1) then
         call soot_write_output('Qnucl.dat',Qnucl,pgaus,nsect_ssm)
      end if

      write(momod(modul) % lun_outpu,*)'# Nucleation module finished. '

    end subroutine chm_source_nucleation_ssm

    subroutine chm_source_condensation_ssm(pgaus,gpcon,gptem,gpden,NDsoot,Qcond,SSTcond)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !    
      ! DESCRIPTION
      !    
      ! USES
      ! USED BY
      !    
      !***
      !-----------------------------------------------------------------------
      use def_chemic , only: W_k, kfl_soot_chm
      implicit none

      integer(ip),  intent(in)  :: pgaus                                ! Vector size
      real(rp),     intent(in)  :: gpcon(pgaus,nspec_chm)               ! Mass fraction species Yk
      real(rp),     intent(in)  :: gptem(pgaus)                         ! Temperature
      real(rp),     intent(in)  :: gpden(pgaus)                         ! Temperature
      real(rp),     intent(in)  :: NDsoot(pgaus,nsect_ssm)              ! NDSoot
      real(rp),     intent(out) :: Qcond(pgaus,nsect_ssm)               ! Source term condensation 
      real(rp),     intent(out) :: SSTcond(pgaus,nspec_chm)             ! Species source term from condensation

      integer(ip)               :: k, isect
      integer(ip)               :: IST_type                             ! Intersectional dynamics type
      real(rp)                  :: gpXk(pgaus,nspec_chm)                ! Molar concentration species Ck (mol/cm3s)
      real(rp)                  :: Betaij_Cond(pgaus)                   ! Collison frequancy factor
      real(rp)                  :: NDPAH(pgaus)                         ! Number density of PAH
      real(rp)                  :: VPAH                                 ! Volume of PAH
      real(rp)                  :: iVol                                 ! volume of section i
      real(rp)                  :: Deltarate(pgaus,nsect_ssm,3)         ! Intersectional matrix

      IST_type = 1

      Qcond   = 0.0_rp
      SSTcond = 0.0_rp

      gpXk = spec_con(pgaus,gpden,gpcon)

      NDPAH = max(gpXk(:,id_A4) * NAvo, 0.0_rp)
      VPAH  = 0.5_rp*NCVmin*(MC/(NAvo*RhoC_ssm))

      do  isect = 1,nsect_ssm

          iVol = Vpmean(isect)
          call soot_collision_freq(pgaus,iVol,VPAH,gptem,Betaij_Cond)

          do k = 1,pgaus
             Deltarate(k,isect,1) = VPAH*Betaij_Cond(k)*(NDPAH(k)*NDsoot(k,isect))
             Deltarate(k,isect,1) = max(Deltarate(k,isect,1),1.0e-32_rp)

             !
             ! Contribution to gas phase  (gm/cm3)
             !
             SSTcond(k,id_A4) = SSTcond(k,id_A4) &
                               -1.0_rp*(Deltarate(k,isect,1)/(VPAH*NAvo))*(W_k(id_A4)*UnC)
             SSTcond(k,id_H2) = SSTcond(k,id_H2) &
                               +5.0_rp*(Deltarate(k,isect,1)/(VPAH*NAvo))*(W_k(id_H2)*UnC)
          end do
      end do

      call soot_intersectdynamics(pgaus,IST_type,Deltarate)

      do k = 1,pgaus
           do isect = 1,nsect_ssm
                if (isect .eq. 1) then
                   Qcond(k,isect) = Deltarate(k,isect,2)
                else
                   Qcond(k,isect) = Deltarate(k,isect-1,3) + &
                                    Deltarate(k,isect,2)
                end if
           end do
      end do

      !
      ! Writing condensation source
      !
      if (kfl_soot_chm .eq. -1) then
         call soot_write_output('Qcond.dat',Qcond,pgaus,nsect_ssm)
      end if

      write(momod(modul) % lun_outpu,*)'# Condensation module finished. '

    end subroutine chm_source_condensation_ssm

    subroutine chm_source_coagulation_ssm(pgaus,gpcon,gptem,gpden,NDSoot,Qcoag)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !    
      ! DESCRIPTION
      !    
      ! USES
      ! USED BY
      !    
      !***
      !-----------------------------------------------------------------------
      use def_chemic , only: kfl_soot_chm

      implicit none

      integer(ip),  intent(in)  :: pgaus                    ! Vector size
      real(rp),     intent(in)  :: gpcon(pgaus,nspec_chm)   ! Mass fraction species Yk
      real(rp),     intent(in)  :: gptem(pgaus)             ! Temperature
      real(rp),     intent(in)  :: gpden(pgaus)             ! Temperature
      real(rp),     intent(in)  :: NDsoot(pgaus,nsect_ssm)  ! NDSoot
      real(rp),     intent(out) :: Qcoag(pgaus,nsect_ssm)   ! Source term condensation

      integer(ip)      :: i, j, k, isect, jsect, ksect, ivol
      integer(ip)      :: ConPAHCoagLoc
      integer(ip)      :: CoagCutSec,CoagCut_init,CoagCut_end
      real(rp)         :: Vol_i, Vol_j, Vol_k, VTot            ! volume of section
      real(rp)         :: ConPAHCoag(pgaus)
      real(rp)         :: rcexp,ctexp
      real(rp)         :: Epsilonij, kfm, Betaij_fm
      real(rp)         :: Beta_ij(nsect_ssm,nsect_ssm)
      real(rp)         :: VkLow, VkUp
      real(rp)         :: split_fact, kron_del, gain, loss
      integer(ip)      :: iCoagOpt

      iCoagOpt = 2           ! Optimization flag (1: Optimized, 2: Standard)

      Qcoag    = 0.0_rp
      Beta_ij  = 0.0_rp

      !
      ! Optimization of coagulation model
      !
      if (iCoagOpt .eq. 1) then
          !
          ! Cutoff section
          !
          !!DMM if (mod(nsect_ssm,2) .eq. 0) then
          CoagCutSec = nsect_ssm / 2
          !!DMM else
          !!DMM    CoagCutSec = (nsect_ssm-1)/2
          !!DMM end if

      else
          CoagCutSec   = nsect_ssm
          CoagCut_init = 1
          CoagCut_end  = pgaus
      end if

!--------------------------------------------------------------

      ctexp= 1.0_rp/3.0_rp
      rcexp= 1.0_rp/2.0_rp

      !
      ! Beta constants for free-molecular regime
      ! van der Wall collision efficiency
      !

      Epsilonij = 2.2_rp

      !
      ! Constant for Free-Molecular Regime
      !
      do k = CoagCut_init,CoagCut_end

         kfm = Epsilonij*((3.0_rp/(4.0_rp*pi))**(1.0_rp/6.0_rp)) &
                        *((6.0_rp*kBoltz*gptem(k)/RhoC_ssm)**(rcexp))

         !
         ! Precalculation of all particle collision: Beta_ij
         !
         do isect = 1, CoagCutSec
            Vol_i = Vpmean(isect)
            do jsect = 1, CoagCutSec
                Vol_j = Vpmean(jsect)
                !
                ! For general collisons in free-molecular regime!
                !
                Betaij_fm = kfm  &
                            *( Vol_i**ctexp + Vol_j**ctexp )**2.0_rp &
                            *( (1.0_rp/Vol_i + 1.0_rp/ Vol_j)**rcexp )

                Beta_ij(isect,jsect) = Betaij_fm
            end do
         end do

         !
         ! Calculation of collisions
         !
         do isect = 1,nsect_ssm
            ksect = isect
            Vol_k = Vpmean(ksect)
            !
            ! Increase Source due to particles coagulation: Vi + Vj = Vk
            ! Coag. of sections i, j <= k, volume -> k
            !
            gain = 0.0_rp

            do i = 1,ksect
               Vol_i = Vpmean(i)
               do j = i,ksect
                  Vol_j = Vpmean(j)
                  if (ksect .eq. 1) then
                     VkLow = Vpmean(ksect)
                  else
                     VkLow = Vpmean(ksect-1)
                  end if

                  if (ksect .eq. nsect_ssm) then
                     VkUp = Vpmean(ksect)
                  else
                     VkUp = Vpmean(ksect+1)  !  (i+1)th section
                  end if

                  VTot = Vol_i + Vol_j

                  !
                  ! Calculation of the spliting factor
                  !

                  if (Vol_k .le. VTot .and. VTot .le. VkUp) then
                     split_fact = (VkUp - VTot) / (VkUp - Vol_k)
                  elseif (VkLow .le. VTot .and. VTot .le. Vol_k) then
                     split_fact = (VkLow - VTot) / (VkLow - Vol_k)
                  else
                     split_fact = 0.0_rp
                  end if

                  if (j .eq. i) then
                     kron_del = 1.0_rp      ! Kronecker delta
                  else
                     kron_del = 0.0_rp
                  end if

                  gain = gain + (1.0_rp-0.5_rp*kron_del)*split_fact*Beta_ij(i,j)&
                                 *NDsoot(k,i)*NDsoot(k,j)
               end do
            end do

            !
            ! Decrease Source due to particles coagulation: Vl + Vp -> Vl
            ! Coag. of sections k, j => k, volume -< k
            !
            loss = 0.0_rp

            do i= 1,nsect_ssm
               loss = loss + Beta_ij(isect,i) &
                       * NDsoot(k,isect)*NDsoot(k,i)
            end do

            ivol = isect+1

            !
            ! Coagulation source term
            !
            Qcoag(k,isect) = (gain - loss) * &
                    ( Vsec(ivol)-Vsec(ivol-1) ) /log( Vsec(ivol)/Vsec(ivol-1) )
         end do
      end do

      Qcoag(:,nsect_ssm) = 0.0_rp

      !
      ! Writing coagulation source
      !
      if (kfl_soot_chm .eq. -1) then
         call soot_write_output('Qcoag.dat',Qcoag,pgaus,nsect_ssm)
      end if

      write(momod(modul) % lun_outpu,*)'# Coagulation module finished. '

    end subroutine chm_source_coagulation_ssm

    subroutine chm_source_surfgrowth_ssm(pgaus,gpcon,gptem,gpden,NDsoot,Asoot,Vsoot,Qsurf,SSTsurf)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !    
      ! DESCRIPTION
      !    
      ! USES
      ! USED BY
      !    
      !***
      !-----------------------------------------------------------------------
      use def_chemic , only: W_k, kfl_soot_chm

      implicit none

      integer(ip),  intent(in)   :: pgaus                       ! Vector size
      real(rp),     intent(in)   :: gpcon(pgaus,nspec_chm)      ! Mass fraction species Yk
      real(rp),     intent(in)   :: gptem(pgaus)                ! Temperature
      real(rp),     intent(in)   :: gpden(pgaus)                ! Density
      real(rp),     intent(in)   :: NDsoot(pgaus,nsect_ssm), &  ! Soot number density
                                    Asoot(pgaus,nsect_ssm),  &  ! Soot surface area density
                                    Vsoot(pgaus,nsect_ssm)      ! Soot volume density
      real(rp),     intent(out)  :: Qsurf(pgaus,nsect_ssm)      ! Source term condensation
      real(rp),     intent(out)  :: SSTsurf(pgaus,nspec_chm)    ! Species source term from surfce growth

      integer(ip)                :: k, isect
      real(rp)                   :: gpXk(pgaus,nspec_chm)       ! Molar concentration species Ck (mol/cm3s)
      real(rp), dimension(pgaus) :: conH, conH2, conOH, conO2, conH2O, conC2H2
      real(rp), dimension(pgaus) :: fR1, rR1, fR2, rR2, fR3, fR4, fR5, fR6
      real(rp), dimension(pgaus) :: RT, denomD, denomC, Ydep, Ycon, Kss
      real(rp), dimension(pgaus) :: sR1, sR2, sR3, sR4, sR5, sR6
      real(rp)                   :: VC, VolAdd_R1, VolAdd_R2, VolAdd_R3, &
                                    VolAdd_R4, VolAdd_R5, VolAdd_R6
      integer(ip)                :: IST_type_R1, IST_type_R2, IST_type_R3, &
                                    IST_type_R4, IST_type_R5, IST_type_R6   ! Intersectional dynamics type
      real(rp)                   :: kP, kR
      real(rp), dimension(pgaus) :: ProdR, ProdP,  &
                                    ProdH, ProdOH,ProdH2, ProdO2, &
                                    ProdCO, ProdH2O,ProdC2H2, ProdHR3, ProdOHR6 ! Production of gas-phase species

      real(rp), dimension(pgaus,nsect_ssm,3) :: Deltarate_sg,   &  ! Intersectional matrix Surf
                                                Deltarate_oxO2, &  ! Intersectional matrix OxO2
                                                Deltarate_oxOH, &  ! Intersectional matrix OxOH
                                                Dummyrate
      real(rp), dimension(pgaus,nsect_ssm)   :: surfSource, oxiSourceO2, oxiSourceOH

      surfSource  = 0.0_rp
      oxiSourceO2 = 0.0_rp
      oxiSourceOH = 0.0_rp
      Qsurf       = 0.0_rp
      SSTsurf     = 0.0_rp

      ProdR    = 0.0_rp    ! Dummy
      ProdP    = 0.0_rp

      ProdH    = 0.0_rp
      ProdOH   = 0.0_rp
      ProdH2   = 0.0_rp
      ProdO2   = 0.0_rp
      ProdCO   = 0.0_rp
      ProdH2O  = 0.0_rp
      ProdC2H2 = 0.0_rp
      ProdHR3  = 0.0_rp     ! Prod of H due to SR3
      ProdOHR6 = 0.0_rp     ! Prod of OH due to SR6

      RT = 1.987e-3_rp *gptem    ! R in kcal/mol

      gpXk = spec_con(pgaus,gpden,gpcon)

      conH    = max(gpXk(:,id_H),0.0_rp)
      conH2   = max(gpXk(:,id_H2),0.0_rp)
      conOH   = max(gpXk(:,id_OH),0.0_rp)
      conO2   = max(gpXk(:,id_O2),0.0_rp)
      conH2O  = max(gpXk(:,id_H2O),0.0_rp)
      conC2H2 = max(gpXk(:,id_C2H2),0.0_rp)

      ! Surface growth is currently described by the following scheme
      ! (Appel et al., 2000, Combust. Flame 121:122-136):
      !       1.  Csoot-H + H = Csoot* + H2         (fR1, rR1)
      !       2.  Csoot-H + OH = Csoot* + H2O       (fR2, rR2)
      !       3.  Csoot*  + H -> Csoot              (fR3)
      !       4.  Csoot*  + C2H2 -> Csoot-H + H     (fR4)
      !       5.  Csoot*  + O2 -> products          (fR5)
      !       6.  Csoot   + OH -> products          (fR6)
      !
      !     Surface reactions coefficients
      !     The following rate coefficients are from Appel et al. (2000)
      !     Ea in (kcal k/mol)

      !
      ! Calculating rates of surface reactions
      !
      fR1 = 4.20E+13_rp *exp(-13.0_rp/RT)*conH
      rR1 = 3.90E+12_rp *exp(-11.0_rp/RT)*conH2
      fR2 = 1.00E+10_rp *(gptem**0.734_rp)*exp(-1.43_rp/RT)*conOH
      rR2 = 3.68E+08_rp *(gptem**1.139_rp)*exp(-17.1_rp/RT)*conH2O
      fR3 = 2.00E+13_rp *conH
      fR4 = 8.00E+7_rp  *(gptem**1.56_rp)*exp(-3.8_rp/RT)*conC2H2
      fR5 = 2.20E+12_rp *exp(-7.50_rp/RT)*conO2
      fR6 = 0.13_rp * conOH                       ! gamma = 0.13 from Neoh et al.

      !
      ! Fraction of radical sites via steady-state:
      !
      denomD = rR1 + rR2 + fR3 + fR4 + fR5     ! For radical depletion
      denomC = rR1 + rR2 + fR3 + fR5           ! For radical conservation

      denomD = max(denomD,1.0e-40_rp)
      denomC = max(denomC,1.0e-40_rp)

      Ydep = (fR1 + fR2)/denomD
      Ycon = (fR1 + fR2)/denomC

      !!DMM This would never happen as we 'max' the values

!!DMM      do k = 1,pgaus
!!DMM          if (denomD(k) .eq. 0.0) then
!!DMM            Ydep(k)  = 0.0_rp
!!DMM          end if
!!DMM
!!DMM          if (denomC(k) .eq. 0.0) then
!!DMM             Ycon(k) = 0.0_rp
!!DMM          end if
!!DMM      end do
      !
      ! Steady-State ratio
      !
      Kss = RadSoot_ssm*(Ycon-Ydep)+Ydep      ! (size pgaus)
      VC  = MC/(NAvo*RhoC_ssm)

      !
      ! For surface reaction-1----------------
      !
      VolAdd_R1 = 0.0_rp
      sR1 = fR1 - rR1 * Kss
      kR = -1.0_rp
      kP = 1.0_rp
      IST_type_R1 = 0

      call soot_surf_reac_rates( pgaus,IST_type_R1,'R1',VolAdd_R1,gptem,sR1, &
                                 Dummyrate,kR,ProdH,kP,ProdH2,NDsoot,Asoot,Vsoot )

      ProdH2 = ProdH2/2.0_rp

      !
      ! For surface reaction-2----------------
      !
      VolAdd_R2 = 0.0_rp
      sR2 = fR2 - rR2 * Kss
      kR = -1.0_rp
      kP = 1.0_rp
      IST_type_R2 = 0

      call soot_surf_reac_rates( pgaus,IST_type_R2,'R2',VolAdd_R2,gptem,sR2, &
                                 Dummyrate,kR,ProdOH,kP,ProdH2O,NDsoot,Asoot,Vsoot )

      ProdH2 = ProdH2 - ProdH2O/2.0_rp

      !
      ! For surface reaction-3----------------
      !
      VolAdd_R3 = 0.0_rp
      sR3 = fR3 * Kss
      kR = -1.0_rp
      kP = 0.0_rp
      IST_type_R3 = 0

      call soot_surf_reac_rates( pgaus,IST_type_R3,'R3',VolAdd_R3,gptem,sR3, &
                                 Dummyrate,kR,ProdHR3,kP,ProdR,NDsoot,Asoot,Vsoot )

      ProdH  = ProdH  + ProdHR3
      ProdH2 = ProdH2 - ProdHR3/2.0_rp

      !
      ! For surface growth reaction ----------------
      !
      VolAdd_R4 = 2.0_rp*VC
      sR4 = fR4* Kss
      kR = -1.0_rp
      kP = 1.0_rp
      IST_type_R4 = 1

      call soot_surf_reac_rates( pgaus,IST_type_R4,'R4',VolAdd_R4,gptem,sR4, &
                                 Deltarate_sg,kR,ProdC2H2,kP,ProdH,NDsoot,Asoot,Vsoot )

      ProdH2 = ProdH2 - ProdC2H2/2.0_rp

      !
      ! For surface OxO2 reaction ----------------
      !
      VolAdd_R5 = -2.0_rp*VC
      sR5 = fR5* Kss
      kR = -1.0_rp
      kP = 2.0_rp
      IST_type_R5 = 2

      call soot_surf_reac_rates( pgaus,IST_type_R5,'R5',VolAdd_R5,gptem,sR5, &
                                 Deltarate_oxO2,kR,ProdO2,kP,ProdCO,NDsoot,Asoot,Vsoot )

      !
      ! For surface OxOH reaction ----------------
      !
      VolAdd_R6 = -1.0_rp*VC
      sR6 = fR6* NAvo
      kR = -1.0_rp
      kP = 1.0_rp
      IST_type_R6 = 2

      call soot_surf_reac_rates( pgaus,IST_type_R6,'R6',VolAdd_R6,gptem,sR6, &
                                 Deltarate_oxOH,kR,ProdOHR6,kP,ProdCO,NDsoot,Asoot,Vsoot )

      ProdOH = ProdOH + ProdOHR6
      ProdH2 = ProdH2 - ProdOHR6/2.0_rp


    do k = 1,pgaus
        do  isect = 1,nsect_ssm
            !
            ! Calculation of surface growth source term
            !
            if (isect .eq. 1) then
               surfSource(k,isect) = Deltarate_sg(k,isect,2)
            else
               surfSource(k,isect) = Deltarate_sg(k,isect-1,3) + &
                                     Deltarate_sg(k,isect,2)
            end if

            !
            ! Calculation of OH oxidation source term
            !
            if (isect .eq. nsect_ssm) then
               oxiSourceO2(k,isect) = Deltarate_oxO2(k,isect,2)
            else
               oxiSourceO2(k,isect) = Deltarate_oxO2(k,isect+1,3) + &
                                      Deltarate_oxO2(k,isect,2)
            end if

            !
            ! Calculation of OH oxidation source term
            !
            if (isect .eq. nsect_ssm) then
               oxiSourceOH(k,isect) = Deltarate_oxOH(k,isect,2)
            else
               oxiSourceOH(k,isect) = Deltarate_oxOH(k,isect+1,3) + &
                                      Deltarate_oxOH(k,isect,2)
            end if

           Qsurf(k,isect) = surfSource(k,isect) + oxiSourceO2(k,isect) + &
                            oxiSourceOH(k,isect)
        end do

        !
        ! Calculation of source terms for species due to surface growth and
        ! soot oxidation (gm/cm3)
        !
        SSTsurf(k,id_H)    = ProdH(k)   *(W_k(id_H)*UnC)
        SSTsurf(k,id_OH)   = ProdOH(k)  *(W_k(id_H)*UnC)
        SSTsurf(k,id_H2)   = ProdH2(k)  *(W_k(id_H)*UnC) 
        SSTsurf(k,id_O2)   = ProdO2(k)  *(W_k(id_H)*UnC)
        SSTsurf(k,id_CO)   = ProdCO(k)  *(W_k(id_H)*UnC)
        SSTsurf(k,id_H2O)  = ProdH2O(k) *(W_k(id_H)*UnC)
        SSTsurf(k,id_C2H2) = ProdC2H2(k)*(W_k(id_H)*UnC)
    end do

    !
    ! Writing surface growth source
    !
    if (kfl_soot_chm .eq. -1) then
       call soot_write_output('Qsurf.dat',Qsurf,pgaus,nsect_ssm)
    end if

    write(momod(modul) % lun_outpu,*)'# Surface growth module finished. ' 

    end subroutine chm_source_surfgrowth_ssm

    subroutine soot_volume_discretization_ssm()
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------
    
      ! This subroutine discretizes the PSDF into number of sections
      implicit none
      integer(ip)         :: i, j
      real(rp)            :: VC, Vmin, VPAH

      RhoC_ssm = RhoC_ssm/UnC
      Vmax_ssm = Vmax_ssm*(Unc**2.0_rp)      ! Vmax in (cm3)
    
      VC   = MC/(NAvo*RhoC_ssm)
      Vmin = NCVmin*VC
      VPAH = Vmin/2.0_rp
    
      !
      ! Volume of the smallest particle
      !
      Vsec(1) = Vmin
    
      !
      ! Mean particle size is linearily discretized in the log scale
      ! Profile of Roy 2014, PhD
      !
      do i = 1,nsect_ssm
         j = i+1
         Vsec(j) = Vmin*(Vmax_ssm / Vmin)**(real(i,rp)/real(nsect_ssm,rp))
      end do

      !
      ! Calculation of the mean particle size in the log scale
      !
      do i = 1,nsect_ssm
         j = i+1
         Vpmean(i) = exp( ( log(Vsec(j))+log(Vsec(j-1)) ) /2.0_rp)
         Dpmean(i) = (6.0_rp*Vpmean(i)/pi)**(1.0_rp/3.0_rp)
      enddo
    
      !
      ! Write volume discretization
      !
      write(momod(modul) % lun_outpu,*)'----------------------------------- '
      write(momod(modul) % lun_outpu,*)' '
      write(momod(modul) % lun_outpu,*)'# Discretization of the PSDF: '
      do i=1,nsect_ssm+1
         write(momod(modul) % lun_outpu,*)'  Section', i-1,' Vsec =',Vsec(i), ' [cm^3]'
      enddo
      write(momod(modul) % lun_outpu,*)' '
      write(momod(modul) % lun_outpu,*)'----------------------------------- '
    
    end subroutine soot_volume_discretization_ssm

    function soot_vol_frac(pgaus, mix_den, soot_massfrac)
    
      ! This function converts soot mass fraction to volume fraction
    
      implicit none
    
      integer(ip), intent(in) :: pgaus
      integer(ip)             :: i,j
      real(rp)                :: mix_den(pgaus)
      real(rp)                :: soot_massfrac(pgaus,nsect_ssm)
      real(rp)                :: soot_vol_frac(pgaus,nsect_ssm)
    
      do i=1,pgaus
          do j=1,nsect_ssm
             soot_vol_frac(i,j) = max(soot_massfrac(i,j),0.0_rp)*(mix_den(i)/(RhoC_ssm*UnC))
          end do
      end do
    
    end function soot_vol_frac

    function soot_vol_frac_dist(pgaus, sootvolfrac)
    
      !
      ! This function converts soot mass fraction to volume fraction density
      !
      integer(ip), intent(in) :: pgaus
      integer(ip)             :: i, j, ivol
      real(rp)                :: sootvolfrac(pgaus,nsect_ssm)
      real(rp)                :: soot_vol_frac_dist(pgaus,nsect_ssm)
    
      do i=1,pgaus
          do j=1,nsect_ssm
             ivol = j+1
             soot_vol_frac_dist(i,j) = sootvolfrac(i,j)/(Vsec(ivol) - Vsec(ivol-1))
          end do
      end do
    
    end function soot_vol_frac_dist

    function soot_num_den(pgaus, svfdist)

      !
      ! This function calculates soot number density in each section
      !
      integer(ip), intent(in) :: pgaus
      integer(ip)             :: i, j, ivol
      real(rp)                :: svfdist(pgaus,nsect_ssm)
      real(rp)                :: soot_num_den(pgaus,nsect_ssm)

      do i=1,pgaus
          do j=1,nsect_ssm
             ivol = j+1
             soot_num_den(i,j) = svfdist(i,j)*log(Vsec(ivol)/ Vsec(ivol-1))
          end do
      end do
    
    end function soot_num_den
    
    function soot_surf_area(pgaus, svfdist)
    
      !
      ! This function calculates soot surface area in each section
      !
      integer(ip), intent(in) :: pgaus
      integer(ip)             :: i, j, ivol
      real(rp)                :: svfdist(pgaus,nsect_ssm)
      real(rp)                :: soot_surf_area(pgaus,nsect_ssm)
      real(rp)                :: ctexp

      ctexp = 2.0_rp/3.0_rp
    
      do i=1,pgaus
          do j=1,nsect_ssm
             ivol = j+1
             soot_surf_area(i,j) = pi*(6.0_rp/pi)**ctexp &
                                * svfdist(i,j) * (1.0_rp/ctexp) &
                                * ( Vsec(ivol)**ctexp - Vsec(ivol-1)**ctexp )
          end do
      end do
    
    end function soot_surf_area
    
    function soot_vol(pgaus, svfdist, numden)
    
      !
      ! This function calculates soot volume in each section
      !
      integer(ip), intent(in) :: pgaus
      integer(ip)             :: i, j, ivol
      real(rp)                :: svfdist(pgaus,nsect_ssm)
      real(rp)                :: numden(pgaus,nsect_ssm)
      real(rp)                :: soot_vol(pgaus,nsect_ssm)
    
      do i=1,pgaus
          do j=1,nsect_ssm
             ivol = j+1
             soot_vol(i,j) = svfdist(i,j) *( Vsec(ivol) - Vsec(ivol-1) ) &
                             / (numden(i,j)+1.0_rp)
          end do
      end do
    
    end function soot_vol

    function spec_con(pgaus, mix_den, spec_massfrac )

      !
      ! This function converts mass fraction to molar concentration
      !
      use def_chemic,  only   : W_k, kfl_soot_chm
      implicit none

      integer(ip)             :: i
      integer(ip), intent(in) :: pgaus
      real(rp),    intent(in) :: mix_den(pgaus)
      real(rp)                :: spec_massfrac(pgaus,nspec_chm)
      real(rp)                :: spec_con(pgaus,nspec_chm)

      do i = 1,nspec_chm
         spec_con(:,i) = (mix_den/UnC)*(spec_massfrac(:,i)/(W_k(i)*UnC)) !to convert from (kg/m3) to (g/cm3)
      end do

    end function spec_con

    subroutine soot_collision_freq(pgaus,Vi,Vj,T,Betaij)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------
    
      ! This subroutine calculates collision coefficients
    
      implicit none
    
      integer(ip), intent(in)  :: pgaus
      real(rp),    intent(in)  :: Vi,Vj
      real(rp),    intent(in)  :: T(pgaus)
      real(rp),    intent(out) :: Betaij(pgaus)
    
      real(rp) :: Epsilonij
      real(rp) :: ctexp,rcexp
      real(rp) :: kfm(pgaus)
    
      ctexp = 1.0_rp/3.0_rp
      rcexp = 1.0_rp/2.0_rp
    
      !
      ! Beta constants for free-molecular regime
      ! van der Wall collision efficiency
      !
      Epsilonij = 2.2_rp
    
      !
      ! Constant for Free-Molecular Regime
      !
      kfm = Epsilonij *((3.0_rp/(4.0_rp*pi))**(1.0_rp/6.0_rp)) &
                      *( (6.0_rp*kBoltz*T/RhoC_ssm)**(rcexp) )
    
      !
      ! For general collisons in free-molecular regime
      !
    
      Betaij = kfm * ( Vi**ctexp + Vj**ctexp )**2.0_rp &
                   * ( (1.0_rp/Vi + 1.0_rp/ Vj)**rcexp )
    
    end subroutine soot_collision_freq
    
    subroutine soot_intersectdynamics(pgaus,ST_type,DeltaQ)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------
    
      ! This subroutine distributes soot mass over sections
    
      implicit none
    
      integer(ip), intent(in) :: pgaus
      integer(ip), intent(in) :: ST_type
      integer(ip)             :: isec, ivol, imeshp
      real(rp)                :: IBout, IBin
      real(rp)                :: DeltaQ(pgaus,nsect_ssm,3)
    
      do imeshp=1,pgaus
          do isec=1,nsect_ssm
             ivol=isec+1

             !
             ! Inter-sectional balance for Surface-Growth
             !
             if (ST_type==1) then

                 !
                 ! Sections i < nsect_ssm-1
                 !
                 if (isec .lt. nsect_ssm) then
                    !
                    ! Inter-sectional Balance In : Dq_in
                    !
                    IBin = ((Vsec(ivol+1)-Vsec(ivol))/(Vsec(ivol)-Vsec(ivol-1)))* &
                             (log(Vsec(ivol)/Vsec(ivol-1))/log(Vsec(ivol+1)/Vsec(ivol)))
    
                    DeltaQ(imeshp,isec,2) = DeltaQ(imeshp,isec,1)/(1.0_rp-IBin)
    
                    !
                    ! Inter-sectional Balance Out : Dq_out
                    !
                    IBout = ((Vsec(ivol)-Vsec(ivol-1))/(Vsec(ivol+1)-Vsec(ivol)))*&
                            (log(Vsec(ivol+1)/Vsec(ivol))/log(Vsec(ivol)/Vsec(ivol-1)))
    
                    DeltaQ(imeshp,isec,3) = DeltaQ(imeshp,isec,1)/(1.0_rp-IBout)
    
                 !
                 ! Last section
                 !
                 elseif  (isec .eq. nsect_ssm) then
                    DeltaQ(imeshp,isec,2) = DeltaQ(imeshp,isec,1)
                    DeltaQ(imeshp,isec,3) = 0.0_rp
                 endif
    
             !
             ! Inter-sectional balance for Oxidation
             !
             elseif (ST_type==2) then
    
                 if (isec .gt. 1) then
                    !
                    ! Intersectional Balance In
                    !
                    IBIn  = ((Vsec(ivol-1)-Vsec(ivol-2))/(Vsec(ivol)-Vsec(ivol-1)))*&
                            (log(Vsec(ivol)/Vsec(ivol-1))/log(Vsec(ivol-1)/Vsec(ivol-2)))
    
                    DeltaQ(imeshp,isec,2) = DeltaQ(imeshp,isec,1)/(1.0_rp-IBIn)
    
                    !
                    ! Intersectional Balance Out
                    !
                    IBOut = ((Vsec(ivol)-Vsec(ivol-1))/(Vsec(ivol-1)-Vsec(ivol-2)))*&
                            (log(Vsec(ivol-1)/Vsec(ivol-2))/log(Vsec(ivol)/Vsec(ivol-1)))
    
                    DeltaQ(imeshp,isec,3) = DeltaQ(imeshp,isec,1)/(1.0_rp-IBOut)
    
                 elseif (isec .eq. 1) then
                    DeltaQ(imeshp,isec,2) = DeltaQ(imeshp,isec,1)
                    DeltaQ(imeshp,isec,3) = 0.0_rp
                 end if
             end if
          end do
      end do
    
    end subroutine soot_intersectdynamics

    subroutine soot_surf_reac_rates(pgaus,ST_type,RName,VolAdd,Temp,Rcoef,DeltaQ,stoiMR,ProdR,stoiMP,ProdP,NDsoot,Asoot,Vsoot)
      !-----------------------------------------------------------------------
      !
      ! DESCRIPTION : This subroutine calculates the contribution of surface reaction
      !               in production rates of species and soot source
      !
      !
      ! INPUTS:
      !         pgaus   = Vector size
      !         ST_type = Type of intersectional dynamics
      !         RName   = Name of surface reaction
      !         VolAdd  = Volume of C added due to surface reaction
      !         Temp    = Gas temperature (K)
      !         Rcoef   = Surface rate of reaction
      !         DeltaQ  = Surface rate with inter-sectional distribution (cm3/cm3)
      !         stoiMR  = stoichimetric coefficient for reactant
      !         ProdR   = Production rate of reactant species (mol/cm3 s)
      !         stoiMR  = stoichimetric coefficient for product
      !         ProdP   = Production rate of product species (mol/cm3 s)
      !         NDsoot  = Number density of soot (# / cm3)
      !         Asoot   = Surface area density of soot (cm2 / cm3)
      !         NDsoot  = Volume of soot (cm3)
      !
      ! OUTPUT:
      !         DeltaQ  = Surface rate with inter-sectional distribution (cm3/cm3)
      !         ProdR   = Production rate of reactant species (mol/cm3 s)
      !         ProdP   = Production rate of product species (mol/cm3 s)
      !
      !-----------------------------------------------------------------------
    
      implicit none
    
      integer(ip), intent(in)        :: pgaus
      integer(ip), intent(in)        :: ST_Type
      character(len = 2), intent(in) :: RName
      real(rp),    intent(in)        :: VolAdd
      real(rp),    intent(in)        :: Temp(pgaus),Rcoef(pgaus)
      real(rp)                       :: DeltaQ(pgaus,nsect_ssm,3)
      real(rp),    intent(in)        :: stoiMR,stoiMP
      real(rp),    intent(in)        :: NDsoot(pgaus,nsect_ssm), &
                                        Asoot(pgaus,nsect_ssm),  &
                                        Vsoot(pgaus,nsect_ssm)
      real(rp)                       :: ProdR(pgaus),ProdP(pgaus)
    
      integer(ip)                    :: isect, k, iStericFac
      real(rp)                       :: Xsoot, VC, VOH, Vol_i
      real(rp), dimension(pgaus)     :: TotalRate,par_a,par_b,nCarbon, &
                                        StericFact,Csoot,RateQ,AS,ModRate, &
                                        Betaij_OH

      iStericFac = 2       ! 1 = unity ; 2 = ABF2000
    
      !
      ! Nominal number of sites per cm^2
      !
      Xsoot = 2.3e+15_rp
    
      !
      ! Initilizing the rate summation for species production
      !
      TotalRate = 0.0_rp
    
      do isect = 1,nsect_ssm
          !
          ! Reaction name
          !
          if (RName .eq. 'R6') then
             VOH = 1.2_rp/NAvo ! Molecular volume of OH
             Vol_i = Vpmean(isect)
             call soot_collision_freq(pgaus,Vol_i,VOH,Temp,Betaij_OH)
             RateQ = Betaij_OH*NDsoot(1:pgaus,isect) / NAvo

          !
          ! Steric factor  definition  (fraction of actie sites that will react)
          !
          else
             if (iStericFac .eq. 1) then
                StericFact = 1.0_rp
             else if (iStericFac .eq. 2) then
                VC = MC/(NAvo*RhoC_ssm)
                par_a = 12.65_rp - 5.63e-03_rp * Temp
                par_b = -1.38_rp + 6.80e-04_rp * Temp
                nCarbon = max(Vsoot(1:pgaus,isect),1.0e-40_rp)/VC            ! Average number of C atoms in a particles of section isec
                StericFact =  tanh((par_a/ log10(nCarbon)) + par_b)
                do k = 1,pgaus
                   if (StericFact(k) .lt. 0.0_rp) then
                       StericFact(k) = 0.0_rp
                   end if
                end do
             end if
             AS = Asoot(1:pgaus,isect)                    ! total soot surface area in section i
             Csoot = Xsoot * AS /NAvo                     ! concentration of sites available
             RateQ = Csoot*StericFact                     ! concentration of reactive sites
          end if
    
          TotalRate = TotalRate + RateQ                   ! summation of [S] over sections
          ModRate   = max(NAvo*(Rcoef*RateQ), 1.0e-32_rp)
          DeltaQ(1:pgaus,isect,1) = VolAdd*ModRate        ! soot surface rate in section i
      end do
    
      !
      ! Inter-sectional distribution
      !
      if (ST_Type .gt. 0) &
         call soot_intersectdynamics(pgaus,ST_Type,DeltaQ)
    
      !
      ! Total specie rate over all sections for Reation
      !
      ProdR = ProdR + stoiMR*(Rcoef*TotalRate)
      ProdP = ProdP + stoiMP*(Rcoef*TotalRate)
    
    end subroutine soot_surf_reac_rates

    subroutine soot_write_output(OutFname,OutVar,nrows,ncols)
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------
      implicit none

      integer(ip)   :: i,nrows,ncols
      character(9)  :: OutFname
      real(rp)      :: OutVar(nrows,ncols)

      ! Writing nucleation source and contribution to gas phase

      open(1, file = OutFname, status='replace')

      do i=1,nrows
         write(1,*) OutVar(i,1:ncols)
      end do

      close(1)

      print*, '    Outputfile ', OutFname,' written!'

    end subroutine soot_write_output

    subroutine soot_read_species()
      !-----------------------------------------------------------------------
      !****f* Chemic/
      ! NAME :w
      !
      !
      ! DESCRIPTION
      !
      ! USES
      ! USED BY
      !
      !***
      !-----------------------------------------------------------------------

#ifdef CANTERA
      use cantera
      use def_chemic,    only : gas_chm
#endif

      implicit none
      integer                                     :: i,j
      character(Len=100)                          :: soot_species(8)
      character(len=:), allocatable               :: soot_species_trimmed
      integer                                     :: index_sootspec(8)

      soot_species(1) = 'A4'
      soot_species(2) = 'H'
      soot_species(3) = 'H2'
      soot_species(4) = 'OH'
      soot_species(5) = 'O2'
      soot_species(6) = 'H2O'
      soot_species(7) = 'C2H2'
      soot_species(8) = 'CO'

#ifdef CANTERA
      do i = 1,8
         soot_species_trimmed = trim(adjustl(soot_species(i)))
         index_sootspec(i) = speciesIndex(gas_chm, soot_species_trimmed)
      end do
#endif
      id_A4   = index_sootspec(1)
      id_H    = index_sootspec(2)
      id_H2   = index_sootspec(3)
      id_OH   = index_sootspec(4)
      id_O2   = index_sootspec(5)
      id_H2O  = index_sootspec(6)
      id_C2H2 = index_sootspec(7)
      id_CO   = index_sootspec(8)

    end subroutine soot_read_species

    subroutine header()
    
        implicit none
    
        print*,' '
        print*,'---------------------------------------------------------------------'
        print*,'---------------------------------------------------------------------'
        print*,' '
        print*,'    SOOT SOURCE TERMS COMPUTATION MODULE '
        print*,' '
        print*,'---------------------------------------------------------------------'
        print*,'---------------------------------------------------------------------'
        print*,' '
        print*,'    Starting computations..'
        print*,' '
    
    end subroutine header
    
    subroutine footer()
    
        implicit none
    
        print*,' '
        print*,'    Computations done!'
        print*,' '
        print*,'---------------------------------------------------------------------'
        print*,'---------------------------------------------------------------------'
        print*,' '
    
    end subroutine footer

end module mod_chm_sectional_soot_model

