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
  integer(ip) :: ivari
  
  if( INOTSLAVE ) then
     ! 
     !  Initializations (defaults)
     !
     kfl_ellen_chm = 0                                ! Minimum element length
     kfl_taust_chm = 1                                ! Tau strategy
     kfl_shock_chm = 0                                ! Shock capturing off
     kfl_stabi_chm = 0                                ! Galerkin
     kfl_limit_chm = 0                                ! No limiter
     kfl_wallc_chm = 0                                ! Flamelet combustion model correction source term at walls (=0 OFF, =1 ON)
     kfl_tibub_chm = 0                                ! Time integration of Bubble function
     kfl_tiacc_chm = 1                                ! First order time integ.
     kfl_tisch_chm = 1                                ! Trapezoidal rule
     kfl_normc_chm = 2                                ! L2 norm for convergence
     kfl_dtcri_chm = 1                                ! Time step criteria
     kfl_negat_chm = 0                                ! Strategy for negative concentrations
     kfl_posit_chm = 0                                ! Strategy for too positive concentrations
     kfl_warni_chm = 1                                ! DEfault warn about points where mass sums to zero
     kfl_split_chm = -1                               ! Splitting algrithm, ORDER = 1,2 
     staco_chm(1)  = 1.0_rp                           ! Diffusive term
     staco_chm(2)  = 1.0_rp                           ! Convective term
     staco_chm(3)  = 1.0_rp                           ! Reactive term
     shock_chm     = 0.0_rp                           ! SC parameter
     bemol_chm     = 0.0_rp                           ! Bemol (convectice term)
     temli_chm     = 0.0_rp                           ! T limiter to compute reaction rates 
     cotol_chm     = 1.0e-3_rp                        ! Convergence tolerance
     safet_chm     = 1.0_rp                           ! Safety factor
     chemical_time_factor = 1.0_rp                    ! Safety factor exclusively for the source term
     cutof_chm     = 1.0e-8_rp                        ! Concentration cutoff for critical time computation
     sstol_chm     = -1.0e-5_rp                       ! Steady-state tolerance
     strec_chm     = 2.0_rp                           ! Adaptive dt: Streatching factor
     dampi_chm     = 2.0_rp                           ! Adaptive dt: damping
     epsht_chm     = 0.025_rp                         ! Adaptive dt: eps_R
     epstr_chm     = 0.025_rp                         ! Adaptive dt: eps_A
     dtmin_chm     = 1.0e-12_rp                       ! Minimum time step
     dtmax_chm     = 1.0e12_rp                        ! Maximum time step
     kfl_temli_chm = 0                                ! Flag to activate a T limiter to compute reaction rates 
     relax_chm = 1.0_rp                               ! Relaxation

     kfl_entropy_chm = 0                              ! entropy viscosity stablization method

     solve_sol     => solve                           ! Solver type
     ivari         =  1
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
           if( words(2) == 'GALER' .or. words(2) == 'OFF  ' ) then
              kfl_stabi_chm = -2
           else if ( words(2) == 'ENTRO') then
              kfl_stabi_chm = -2
              kfl_entropy_chm = 1
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
           else if(words(2)=='CONST' ) then ! reads stab constants
              staco_chm(1) = param(2)
              staco_chm(2) = param(3)
              staco_chm(3) = param(4)
           end if

        else if( words(1) == 'ELEME' ) then
           call realen(kfl_ellen_chm)

        else if( words(1) == 'SPLIT' ) then
              kfl_split_chm = getint('ORDER',kfl_split_chm,'#Splitting algorithm order')
              write(momod(modul) % lun_outpu,*)''
              write(momod(modul) % lun_outpu,*)'SPLITTING ORDER = ',kfl_split_chm
              write(momod(modul) % lun_outpu,*)''

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

        else if( words(1) == 'RELAX' ) then
           if(words(2)=='ON   ' ) then
              relax_chm = getrea('PARAM',0.5_rp,     '#Relaxation factor for update')
           else
              relax_chm = 1.0_rp
           endif

        else if( words(1) == 'ELEME' ) then
           call realen(kfl_ellen_chm)
           if(words(2)=='NEW  ') kfl_ellen_chm=-1

        else if( words(1) == 'TIMEI' ) then
           if(exists('TRAPE').or. exists('BDF  ') ) call runend('In chm_reanut: implicit schemes not available anymore')

           if(exists('ADAMS')) then
              kfl_tisch_chm    = 3
           end if

           if(exists('RUNGE')) then
              kfl_tisch_chm    = 4
           end if

           kfl_tiacc_chm = getint('ORDER',1_ip,'#Time integration order')

        else if( words(1) == 'TIMEA' ) then
           kfl_tiacc_chm = int(param(1),ip)
           if(kfl_timei_chm==0) kfl_tiacc_chm = 1

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

        else if( words(1) == 'CONSI' ) then
           !
           ! Consistent matrix
           !
           ivari = 2

        else if( words(1) == 'ALGEB' ) then
           solve_sol => solve(ivari:)
           call reasol(1_ip)

        else if( words(1) == 'LIMIT' ) then
           kfl_temli_chm = 1
           temli_chm = getrea('TEMPE',-1000.0_rp,'#Temperature limiter to compute reaction rates')

        else if( words(1) == 'PRECO' ) then 
           call reasol(2_ip)

        end if
     end do

  end if

end subroutine chm_reanut
