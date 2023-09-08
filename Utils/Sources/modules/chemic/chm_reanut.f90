subroutine chm_reanut
  !-----------------------------------------------------------------------
  !****f* partis/chm_reanut
  ! NAME 
  !    chm_reanut
  ! DESCRIPTION
  !    This routine reads the numerical treatment
  ! USES
  !    ecoute
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_solver
  use def_chemic
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none

  if( INOTSLAVE ) then
     ! 
     !  Initializations (defaults)
     !
     kfl_dttyp_chm = 1                                ! Module time step strategy

     kfl_ellen_chm = 0                                ! Minimum element length
     kfl_taust_chm = 1                                ! Tau strategy
     kfl_shock_chm = 0                                ! Shock capturing off
     kfl_stabi_chm = -2                               ! Galerkin
     kfl_limit_chm = 0                                ! No limiter
     kfl_wallc_chm = 0                                ! CFI combustion model correction source term at walls (=0 OFF, =1 ON)
     kfl_tibub_chm = 0                                ! Time integration of Bubble function
     kfl_tiacc_chm = 1                                ! First order time integ.
     kfl_stagg_chm = 0                                ! 1 is for staggered step in combustion
     kfl_assem_chm = 1                                ! Assembly strategy
     miinn_chm     = 1                                ! Max inner iterations
     neule_chm     = 0                                ! Number of Euler time steps
     kfl_tisch_chm = 1                                ! Trapezoidal rule
     kfl_normc_chm = 2                                ! L2 norm for convergence
     kfl_coupl_chm = 0                                ! Decoupled equations
     kfl_dtcri_chm = 1                                ! Time step criteria
     kfl_dttar_chm = -1                               ! Time step target
     kfl_sgsti_chm = 0                                ! No SGS tracking
     kfl_negat_chm = 0                                ! Strategy for negative concentrations
     kfl_posit_chm = 0                                ! Strategy for too positive concentrations
     kfl_warni_chm = 1                                ! DEfault warn about points where mass sums to zero
     staco_chm(1)  = 1.0_rp                           ! Diffusive term
     staco_chm(2)  = 1.0_rp                           ! Convective term
     staco_chm(3)  = 1.0_rp                           ! Reactive term
     shock_chm     = 0.0_rp                           ! SC parameter
     bemol_chm     = 0.0_rp                           ! Bemol (convectice term)
     temli_chm     = 0.0_rp                           ! T limiter to compute reaction rates 
     cotol_chm     = 1.0e-3_rp                        ! Convergence tolerance
     safet_chm     = 1.0_rp                           ! Safety factor
     chemical_time_factor = 1.0_rp                    ! Safety factor exclusively for the source term
     cutof_chm     = 1.0e-8                           ! Concentration cutoff for critical time computation
     sstol_chm     = -1.0e-5_rp                       ! Steady-state tolerance
     strec_chm     = 2.0_rp                           ! Adaptive dt: Streatching factor
     dampi_chm     = 2.0_rp                           ! Adaptive dt: damping
     epsht_chm     = 0.025_rp                         ! Adaptive dt: eps_R
     epstr_chm     = 0.025_rp                         ! Adaptive dt: eps_A
     dtmin_chm     = 1.0e-12_rp                       ! Minimum time step
     dtmax_chm     = 1.0e12_rp                        ! Maximum time step
     kfl_spite_chm = 1                                ! Intra species iterations to improve shock capturing
     kfl_meshi_chm = 0                                ! Mesh interpolator activation (OFF=0,ON=1)
     kfl_temli_chm = 0                                ! Flag to activate a T limiter to compute reaction rates 
     relax_chm = 1.0_rp                               ! Relaxation

     kfl_gauss_chm = 1_ip                             ! Level of Gauss-Seidel updating 

     ! Staggered time step: Solve for ADR of concentrations, and then solve for chemical reactions
     ! The variable  kfl_stagg_chm is set in chm_reaphy
     initial_fraction_step_chm=4                    ! Initial try time step will be Delta t divided by this number
     max_fixed_point_iterations_chm=10000             ! Maximum iterations in the fixed point search of the ODE integrator
     fixed_point_tolerance_chm=1.0D-6                 ! Tolerance for fixed point iterations
     odeint_tolerance_chm=1.0D-4                      ! Tolerance for adaptive time stepping
     timestep_min_chm=1.0D-5                          ! Minimum time step allowed in the adaptive time stepping


     solve_sol     => solve                           ! Solver type
     !
     ! Reach the section
     !
     call ecoute('chm_reanut')
     do while(words(1)/='NUMER')
        call ecoute('chm_reanut')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDNU')
        call ecoute('chm_reanut')

        !----------------------------------------------------------------
        !
        ! Stabilization strategy
        !
        !----------------------------------------------------------------

        if( words(1) == 'TAUST' ) then
           call reatau(kfl_taust_chm)

        else if( words(1) == 'STRAT' .or. words(1) == 'STABI' ) then
           if( words(2) == 'GALER' ) then
              kfl_stabi_chm = -2
           else if( words(2) == 'SU   ' .or. words(2) == 'FIRST' ) then
              kfl_stabi_chm = -1
           else if( words(2) == 'ASGS ' ) then
              kfl_stabi_chm =  0
           else if( words(2) == 'FULLO' ) then
              kfl_stabi_chm =  1
           else if( words(2) == 'OSS  ' .or. words(2) == 'AOSS ') then
              kfl_stabi_chm =  2
              if( exists('NOLIM') ) kfl_limit_chm =  0  ! No limiter
              if( exists('SOTO ') ) kfl_limit_chm =  1  ! Soto limiter
              if( exists('DIFFU') ) kfl_limit_chm =  2  ! Very diffusive limiter
              if( exists('FIRST') ) kfl_limit_chm = -1  ! First order
           else if( words(2) == 'AROSS' ) then
              kfl_stabi_chm =  3
              if( exists('NOLIM') ) kfl_limit_chm =  0  ! No limiter
              if( exists('SOTO ') ) kfl_limit_chm =  1  ! Soto limiter
              if( exists('DIFFU') ) kfl_limit_chm =  2  ! Very diffusive limiter
              if( exists('FIRST') ) kfl_limit_chm = -1  ! First order
           else if( words(2) == 'MARGA' ) then
              kfl_stabi_chm =  4
           else if( words(2) == 'SUPG ' ) then
              kfl_stabi_chm =  -3
           else if( words(2) == 'BUBBL' ) then
              kfl_stabi_chm =  -4
              if( words(3) == 'TIMET' ) kfl_tibub_chm = 1
           end if

        else if( words(1) == 'ELEME' ) then
           call realen(kfl_ellen_chm)

        else if( words(1) == 'SHOCK' ) then
           if(exists('ISOTR').or.exists('ON   ')) then
              kfl_shock_chm = 1
              shock_chm  = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           else if(exists('ANISO')) then
              kfl_shock_chm = 2
              shock_chm  = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           end if

        else if( words(1) == 'TIMES' ) then
           if(words(2)=='LOCAL' ) then
              kfl_dttyp_chm = 1
           else if(words(2)=='ADAPT' ) then
              kfl_dttyp_chm = 2
              strec_chm = getrea('STRET',2.0_rp,     '#Adaptive dt: Stretching factor')
              dampi_chm = getrea('DAMPI',2.0_rp,     '#Adaptive dt: damping')
              epsht_chm = getrea('EPSR ',0.025_rp,   '#Adaptive dt: eps_R')
              epstr_chm = getrea('EPSA ',0.025_rp,   '#Adaptive dt: eps_A')
              dtmin_chm = getrea('DTMIN',1.0e-12_rp, '#Adaptive dt: minimum dt')
              dtmax_chm = getrea('DTMAX',1.0e12_rp,  '#Adaptive dt: maximum dt')
              if( exists('BOTHC') ) then
                 kfl_dtcri_chm = 2
              else
                 kfl_dtcri_chm = 1
              end if
           else if(words(2)=='PID' ) then
              kfl_dttyp_chm = 3
           end if

        else if( words(1) == 'NEGAT' ) then
           if(words(2)=='ON   ' ) then
              if (exists('LAGRA')) then
                 kfl_negat_chm = 0_ip  ! no negative strategy
                 kfl_norma_chm = -2_ip  ! Lagrangian multiplier normalization
              else
                 kfl_negat_chm = 1_ip  ! take previous steps
              endif
           else
              kfl_negat_chm = 0_ip     ! Leave as is, no strategy
           endif

        else if( words(1) == 'POSIT' ) then
           if(words(2)=='ON   ' ) then
              kfl_posit_chm = 1
           end if

        else if( words(1) == 'WARNI' ) then
           if(words(2)=='ON   ' ) then
              kfl_warni_chm = 1
           else if(words(2)=='OFF  ' ) then
              kfl_warni_chm = 0
           end if

        else if( words(1) == 'GAUSS' ) then
           kfl_gauss_chm = getint('LEVEL',0_ip,  'Level of Gauss-Seidel update')

        else if( words(1) == 'STAGG' ) then
           if (exists('TEMPER')) then
              kfl_stagg_chm = 2                                    ! 2 is for staggered step in chemic+temper
           else
              kfl_stagg_chm = 1                                    ! 1 is for staggered step in combustion
           endif
           initial_fraction_step_chm      = getint('INITI',4_ip,'#Staggered Combustion: Initial time step fraction')
           max_fixed_point_iterations_chm = getint('MAXIM',10000_ip,&
                '#Staggered Combustion: Max number of fixed point iterations')
           fixed_point_tolerance_chm      = getrea('FIXED ',1.0e-6_rp,'#Staggered Combustion: Fixed point tolerance')
           odeint_tolerance_chm           = getrea('TOLER ',1.0e-4_rp,'#Staggered Combustion: Adaptive Time Step tolerance')
           timestep_min_chm               = getrea('DTMIN',1.0e-5_rp,'#Staggered Combustion: minimum dt')

        else if( words(1) == 'RELAX' ) then
           if(words(2)=='ON   ' ) then
              relax_chm = getrea('PARAM',0.5_rp,     '#Relaxation factor for update')
           else
              relax_chm = 1.0_rp
           endif

        else if( words(1) == 'INTRA' ) then
           kfl_spite_chm = int(param(1))

        else if( words(1) == 'TARGE' ) then
           if( words(2) == 'ALLSP' ) then
              kfl_dttar_chm = -1
           else if( words(2) == 'MAXIM' ) then
              kfl_dttar_chm = getint('MAXIM',10_ip,  '#Adaptive dt: number of targets')
           end if

        else if( words(1) == 'MAXIM' ) then
           miinn_chm = int(param(1))

        else if( words(1) == 'ELEME' ) then
           call realen(kfl_ellen_chm)
           if(words(2)=='NEW  ') kfl_ellen_chm=-1

        else if( words(1) == 'TIMEI' ) then
           if(exists('TRAPE')) kfl_tisch_chm=1
           if(exists('BDF  ')) kfl_tisch_chm=2
           kfl_tiacc_chm = getint('ORDER',1_ip,'#Time integration order')
           neule_chm     = getint('EULER',0_ip,'#EULER TIME STEPS')

        else if( words(1) == 'ASSEM' ) then
           if( words(2) == 'SPLIT' ) kfl_assem_chm = 2
           if( exists('SAVE ') ) kfl_assem_chm = 3

        else if( words(1) == 'TIMEA' ) then
           kfl_tiacc_chm = int(param(1))
           if(kfl_timei_chm==0) kfl_tiacc_chm = 1
           neule_chm = getint('EULER',0_ip,'#EULER TIME STEPS')

        else if( words(1) == 'CONVE' ) then
           cotol_chm = getrea('CONVE',1.0e-3_rp,'#CONVERGENCE TOLERANCE')

        else if( words(1) == 'SAFET' ) then
           safet_chm = param(1)
           if (exists('SOURC')) chemical_time_factor = getrea('SOURC',1.0_rp,'#Source term safety factor')

        else if( words(1) == 'CONCE' ) then
           cutof_chm = param(1)

        else if( words(1) == 'STEAD' ) then
           sstol_chm = param(1)

        else if (words(1) == 'WALLC' ) then
            if(words(2)=='ON   ' )  kfl_wallc_chm = 1_ip 

        else if( words(1) == 'NORMO' ) then
           if(exists('L1   ')) then
              kfl_normc_chm = 1
           else if(exists('L-INF')) then
              kfl_normc_chm = 0
           else if(exists('LINF ')) then
              kfl_normc_chm = 0
           else if(exists('L2   ')) then
              kfl_normc_chm = 2
           else if(exists('ALGEB')) then
              kfl_normc_chm = 3
           end if

        else if( words(1) == 'BEMOL' ) then
           bemol_chm = getrea('BEMOL',0.0_rp,'#Bemol of convective term')

        else if( words(1) == 'COUPL' ) then
           if(words(2)=='ON   ' ) then
              kfl_coupl_chm = 1
           else if(words(2)=='MONOL' ) then
              kfl_coupl_chm = 2
           else
              kfl_coupl_chm = 0
           end if

        else if( words(1) == 'ALGEB' ) then
           call reasol(1_ip)

        else if( words(1) == 'MESHI' ) then
           if(exists('ON   ') ) kfl_meshi_chm = 1

        else if( words(1) == 'LIMIT' ) then
           kfl_temli_chm = 1
           temli_chm = getrea('TEMPE',-1000.0_rp,'#Temperature limiter to compute reaction rates')

        else if( words(1) == 'PRECO' ) then 
           call reasol(2_ip)

        else if( words(1) == 'TRACK' ) then
           if( exists('TIME ') ) kfl_sgsti_chm = 1

        end if
     end do

     if (kfl_coupl_chm == 2) then !Monolithic turns off other options
        kfl_spite_chm = 1  ! Intra species iterations to improve shock capturing
        kfl_gauss_chm = 0_ip !Gauss-Seidel is zero
        kfl_stagg_chm=0_ip
     endif

  end if

end subroutine chm_reanut
