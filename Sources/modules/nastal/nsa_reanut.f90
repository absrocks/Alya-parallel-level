subroutine nsa_reanut
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_reanut
  ! NAME 
  !    nsa_reanut
  ! DESCRIPTION
  !    This routine reads the numerical treatment 
  ! USES
  !    ecoute
  ! USED BY
  !    nsa_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk
  use      def_solver

  use      def_nastal
  use mod_ecoute, only :  ecoute
  use mod_messages, only : livinf

  implicit none
  integer(ip) :: irang,ivari,iauxi
  character(300)           :: messa

  irang= 0


  if(INOTSLAVE) then
     !
     !  Initializations (defaults)
     !
     kfl_timei_nsa = 1                                    ! Transient-like scheme
     kfl_timet_nsa = 1                                    ! Explicit 
     kfl_matri_nsa = 0                                    ! RHS straigthly computed, no matrix 
     kfl_algor_nsa = 1                                    ! Monolithic
     kfl_diagi_nsa = 0                                    ! No diagonal terms appart from closed integ. mass matrix
     kfl_lotim_nsa = 0                                    ! No diagonal terms from local time step
     kfl_tisch_nsa = 1                                    ! Forward-Euler time scheme
     kfl_tiacc_nsa = 1                                    ! Default in nastal for the explicit case
!!!     kfl_tiacc_nsa = 2                                    ! Default in nastal for the explicit case
     kfl_stabi_nsa = 0                                    ! No stabilization at all
     kfl_stafl_nsa = 0                                    ! Stabilization term formulation using the Jacobian matrixes
     kfl_galer_nsa = 0                                    ! Galerkin terms with jacobian change of variables 
     kfl_taufa_nsa = 1                                    ! All tau factors active: convection, reaction, acoustics, viscosity
     kfl_repro_nsa = 0                                    ! Res. projection not used 
     kfl_hconv_nsa = 1                                    ! Correct hconv according to speed vector
     kfl_shock_nsa = 0                                    ! Shock capturing off
     shock_nsa     = 0.0_rp                               ! SC parameter
     shtol_nsa     = 0.0_rp                               ! SC residual tolerance
     kfl_weigh_nsa = 1                                    ! dT/dt is in the residual
     sstol_nsa     = - 1.0_rp                             ! Steady-satate tolerance. Default: NEVER converge
     cotol_nsa     = - 1.0_rp                             ! Subiterations tolerance. Default: NEVER reach convergence
     corat_nsa     = - 1.0_rp                             ! Subiterations tolerance ratio.
     kfl_adres_nsa = 0                                    ! Subiterations adaptive tolerance. DEFAULT: no adaptive
     kfl_normc_nsa = 2                                    ! L2 norm for convergence
     miinn_nsa     = 1                                    ! One internal iteration
     miinn_pseud_nsa = 1                                  ! One pseudo-time iteration
     kfl_algso_nsa = 8                                    ! Algebraic solver is GMRES
     kfl_penal_nsa = 0                                    ! No penalization
     penal_nsa     = 0.0_rp                               ! Penalization factor
     kfl_linea_nsa = 1                                    ! Linearization (RHS=0, Picard=1, Newton=2)
     kfl_ximpl_nsa = 0                                    ! No explicit terms in the implicit
     kfl_delun_nsa = 0                                    ! Delta form is off
     minew_nsa     = 1                                    ! Initial Newton iteration
     npica_nsa     = 1                                    ! Number of Picard iteration (Newton's lin.)   
     kfl_modfi_nsa = 0                                    ! Do not modify kfl_fixno_nsi
     kfl_dttyp_nsa = 0                                    ! Time increment type set to global for all eqs...
     kfl_dtadj_nsa = 0                                    ! ... and for all times.
     kfl_taudi_nsa = 0                                    ! Diagonal or non-diagonal tau
     kfl_track_nsa = 0                                    ! No subscales tracking

     kfl_higha_nsa = 0                                    ! No high aspect ratio elements
     kfl_resmo_nsa = 0                                    ! No residual smoothing
     kfl_skews_nsa = 0

     kranr_nsa     = 0

     kfl_unkse_nsa = 0                                    ! Conservative unknowns set
     kfl_lopre_nsa = 0                                    ! Local preconditioning inactive
     kfl_pseud_nsa = 0                                    ! No pseudo time step

     kfl_rayle_nsa = 0                                    ! No Rayleigh dumping

     kfl_tiext_nsa = 0                                    ! Externally fixed time step flag
     dtext_nsa     = 0.0_rp                               ! Externally fixed time step

     safet_nsa     = 1.0_rp                               ! Safety factor
     safet_pseud_nsa     = 1.0_rp                        ! Safety factor for the pseudo time step
     safrk_nsa     = 1.0_rp                               ! Delta time shrink for RK schemes

     safex_nsa     = 1.0_rp                               ! Time function parameter for safety factor
     safma_nsa     = 1.0e9_rp                             ! Maximum safety factor
     nisaf_nsa     = 0_ip
     safet_table_nsa = 0.0_rp
     kfl_safet_table_nsa = 0

     theta_nsa     = 0.0_rp                               ! Explicitness weigthing factors

     mfreq_nsa     = 1                                   ! Base frequency mode for subscales
     frmax_nsa     = 1                                    ! Maximum frequency mode for subscales

     kfl_fasts_nsa = 0
     nfrap_nsa=1                                          ! no fractional rk 
     kfl_nacdr_nsa = 0                                    ! no nastal cdr
     kfl_reate_nsa = 0                                    ! no reate term neither in tau nor dt

     dtlim_nsa = 1.0_rp
     xfree_nsa = 100000000.0_rp

     kfl_cfllo_nsa = 0
     mcfll_nsa = 0
     fcfll_nsa = 1.0_rp

     neule_nsa = 1                                        ! number of initial Euler timesteps for the BDF or CN

     kfl_mod_elmop_nsa = 0                                ! No elmoperations (debug)

       
     !
     ! Local variables
     !
     ivari = 1
     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute('reanut')
     do while(words(1)/='NUMER')
        call ecoute('reanut')
     end do
     !
     ! Begin to read data
     !
     !-----------------------------------------------------------------------
     ! ADOC[0]> $-----------------------------------
     ! ADOC[0]> $-- Numerical treatment
     ! ADOC[0]> $-----------------------------------
     ! ADOC[0]> NUMERICAL_TREATMENT
     !-----------------------------------------------------------------------
     messa = &
          '        READING NUMERICS...'
     call livinf(0_ip,messa,one)

     do while(words(1)/='ENDNU')
        call ecoute('reanut')
        

        if(words(1)=='GALER') then
           if(words(2)=='JACOB') then
              kfl_galer_nsa = 0
           else if(words(2)=='FLUXE') then
              kfl_galer_nsa = 1
           end if
        else if(words(1)=='HIGHA') then
           kfl_higha_nsa= 1
        else if(words(1)=='RESID') then
           kfl_resmo_nsa= 1
        else if(words(1)=='ELMOP') then
           kfl_mod_elmop_nsa= 1
        else if(words(1)=='SKEWS') then  ! SKEW_SYMMETRIC 
           kfl_skews_nsa= 1
        else if(words(1)=='STABI') then

           if(exists('CHARA')) kfl_stabi_nsa = 1    !Characteristic galerkin
           if(exists('LOCTA')) kfl_stabi_nsa = 2    !Characteristic galerkin, with local pressure tau

           if(exists('MULTI')) then
              kfl_stabi_nsa = 5    !Variational multiscale
              kfl_taudi_nsa = 1                       !Diagonal tau
              if(exists('1DNON')) then
                 kfl_taudi_nsa = 5   !Non-diagonal tau, forced isentropic, forced 1D
                 if (kfl_isent_nsa == 0) then
                    kfl_taudi_nsa = 6   !Non-diagonal tau, non-forced, forced 1D
                 end if
              end if
              if(exists('NONDI')) then
                 kfl_taudi_nsa = 7   !Non-diagonal tau, forced isentropic, multidimensional
                 if (kfl_isent_nsa == 0) then
                    kfl_taudi_nsa = 8   !Non-diagonal tau, non-forced, multidimensional
                 end if
              end if
              if(exists('FLUXE')) kfl_stafl_nsa = 1   !Stabilization convective term with fluxes 
              if(exists('STATI')) kfl_taudi_nsa = 10  !Non-diagonal tau, stationary
              if(exists('TRACK')) kfl_track_nsa = 1   !Subscales tracking
           else if(exists('GENSU')) then
              kfl_stabi_nsa = 6    !Generalized SUPG
              kfl_taudi_nsa = 1                       !Diagonal tau              
              
           end if

        else if(words(1)=='TAUFA') then

           if(exists('NOCON')) kfl_taufa_nsa(1,2) = 0    !Tau factors no-convection
           if(exists('NOREA')) kfl_taufa_nsa(2,2) = 0    !Tau factors no-reaction
           if(exists('NOACO')) kfl_taufa_nsa(3,2) = 0    !Tau factors no-acoustics
           if(exists('NOVIS')) kfl_taufa_nsa(4,2) = 0    !Tau factors no-viscosity

        else if(words(1)=='TIMEF') then

           if(exists('NOCON')) kfl_taufa_nsa(1,1) = 0    !Time factors idem
           if(exists('NOREA')) kfl_taufa_nsa(2,1) = 0    !Time factors
           if(exists('NOACO')) kfl_taufa_nsa(3,1) = 0    !Time factors
           if(exists('NOVIS')) kfl_taufa_nsa(4,1) = 0    !Time factors

        else if(words(1)=='HCONV') then

           if(exists('INACT')) kfl_hconv_nsa = 0    !No correction for hconv 

        else if(words(1)=='SHOCK') then
           shock_nsa= 0.7_rp
           if(exists('ANISO')) then   ! default value
              kfl_shock_nsa = 1
              shock_nsa     = getrea('VALUE',0.7_rp,'Shock capturing parameter')
           end if
           if(exists('ISOTR').or.exists('ON   ')) then
              kfl_shock_nsa = 2
              shock_nsa     = getrea('VALUE',0.7_rp,'Shock capturing parameter')
           end if
           if (exists('TOLER')) then
              shtol_nsa     = getrea('TOLER',0.0_rp,'Shock capturing residual tolerance')
           end if
           if (exists('STATI')) kfl_shock_nsa = 11  ! stationary residual, no transient terms
           if (exists('TRANS')) kfl_shock_nsa = 1   ! transient terms (DEFAULT VALUE)
           
           if (exists('NOMOM')) kfl_shock_nsa(1:ndime) = 0 
           if (exists('NOCON')) kfl_shock_nsa(1+ndime) = 0 
           if (exists('NOENE')) kfl_shock_nsa(2+ndime) = 0 

        else if(words(1)=='UNKNO') then                    ! Unknowns 
           if(words(2)  =='CONSE') kfl_unkse_nsa=0         !   0. conservative set rho-U-E
           if(words(2)  =='PRIMI') kfl_unkse_nsa=1         !   1. primitive set p-u-T
           if(words(2)  =='HEATC') kfl_unkse_nsa=10        !  10. heat and conservative set rho-U-theta
           
        else if(words(1)=='LOCPR') then  ! Local preconditioning          
           do while(words(1)/='ENDLO')
              call ecoute('reanut')
              if(words(1)=='PRECO') then  ! Local preconditioner type
                 if(words(2)  =='IDENT') then     !   1: No preconditioner (Identity preconditioner)
                    kfl_lopre_nsa=1
                 else if(words(2)  =='VLR  ') then        !   VAN LEER-LEE-ROE - steady inviscid problems
                    kfl_lopre_nsa=2
                 else if(words(2)  =='CM   ') then        !   CHOI & MERKLE - steady/unsteady viscous/inviscid problems 
                    kfl_lopre_nsa=3
                 end if
              end if
           end do

        else if(words(1)=='MARCH') then
           if(exists('FASTS').or.exists('ON   ').or.exists('ACTIV')) then
              kfl_fasts_nsa = getint('LAYER',1_ip,'Number of neighbors layers')
           end if           

        else if(words(1)=='TIMES') then        
           if(words(2)  =='ALLGL') kfl_dttyp_nsa    =  0            ! All eqs. global
           if(words(2)  =='ALLLO') kfl_dttyp_nsa    =  1            ! All eqs. local
           if(words(2)  =='GLOBA') kfl_dttyp_nsa(1) =  0            ! Momentum eq.
           if(words(2)  =='LOCAL') kfl_dttyp_nsa(1) =  1
           if(words(3)  =='GLOBA') kfl_dttyp_nsa(2) =  0            ! Continuity eq.
           if(words(3)  =='LOCAL') kfl_dttyp_nsa(2) =  1
           if(words(4)  =='GLOBA') kfl_dttyp_nsa(3) =  0            ! Energy eq.
           if(words(4)  =='LOCAL') kfl_dttyp_nsa(3) =  1
           if(words(5)  =='GLOBA') kfl_dttyp_nsa(4) =  0            ! ... eq.
           if(words(5)  =='LOCAL') kfl_dttyp_nsa(4) =  1
           if(exists('GRADUA'))    kfl_dtadj_nsa    =  1            ! Gradual adjustment  
           dtlim_nsa= getrea('LIMIT',1.0_rp,'Local time step limiter')
           kfl_lotim_nsa = kfl_dttyp_nsa(1) + kfl_dttyp_nsa(2) + kfl_dttyp_nsa(3) + kfl_dttyp_nsa(4)

        else if(words(1)=='CFLLO') then        
           kfl_cfllo_nsa= 1
           iauxi= 0
           call ecoute('nsa_reanut')
           do while(words(1)/='ENDCF')
              if(words(1)=='TOTAL') then  ! Total number of mach values
                 mcfll_nsa= getint('TOTAL',0_ip,'Total number of CFL values')                 
              else
                 iauxi= iauxi + 1
                 if (iauxi > 500) call runend('NSA_REANUT: LOCAL CFL READS LARGER 500. REDUCE NUMBER!')
                 fcfll_nsa(iauxi,1)= param(1)
                 fcfll_nsa(iauxi,2)= param(2)
!! print*, "DEBUG: ", iauxi, mcfll_nsa, fcfll_nsa(iauxi,1), fcfll_nsa(iauxi,2)
          !!    if (iauxi==10) then
          !!       write(789,*) param(1),param(2)                    
          !!       iauxi= 0
          !!    end if
                 if (iauxi > mcfll_nsa) call runend('NSA_REANUT: LOCAL CFL READS LARGER THAN DECLARED')
              end if
              call ecoute('nsa_reanut')
           end do

        else if( words(1) == 'MODIF' ) then
           !
           ! Users defined boundary conditions
           !
           if( words(2) == 'FREE ' ) then
              kfl_modfi_nsa = 1
              xfree_nsa = getrea('XFREE',0.0_rp,'#X Coordinate of the plane where to free')
           end if
           if( words(2) == 'PRNAC') kfl_modfi_nsa=2
           if( words(2) == 'PRNA3') kfl_modfi_nsa=3
           if( words(2) == 'PRNA4') kfl_modfi_nsa=4
           if( words(2) == 'PRNA5') kfl_modfi_nsa=5
           if( words(2) == 'PRNA6') kfl_modfi_nsa=6

        else if(words(1)=='SAFET') then
           safet_nsa = param(1)
           if (exists('EXPON')) &
                safex_nsa  = getrea('EXPON',1.0_rp,'#Safety Factor Exponential time function')
           if (exists('MAXIM')) &
                safma_nsa  = getrea('MAXIM',1.0e9_rp,'#Maximum safety factor function')
           if (exists('INITI')) &
                nisaf_nsa  = getint('INITI',0_ip    ,'#Initial step for variable CFL')
           if (exists('TABLE')) then
              iauxi= 0
              call ecoute('nsa_reanut')
              do while(words(1)/='ENDTA')
                 iauxi= iauxi+1
                 safet_table_nsa(1,iauxi)  = getrea('VALUE',1.0_rp,'#Safety Factor')
                 safet_table_nsa(2,iauxi)  = getrea('RESID',0.0_rp,'#Residual')
                 if (iauxi == 10) &
                      call runend('NSA_REANUT: SAFET TABLE MUST BE SMALLER THAN 10.')
                 call ecoute('nsa_reanut')
              end do
              safet_nsa= safet_table_nsa(1,1)  ! starting safety factor value
              kfl_safet_table_nsa = iauxi
           end if
        else if(words(1)=='NSCDR') then                     ! nastal cdr
           kfl_nacdr_nsa= 1
        else if(words(1)=='TAUAN') then                     ! nastal cdr
           kfl_reate_nsa= 1
        else if(words(1)=='STEAD') then
           sstol_nsa = param(1)
        else if(words(1)=='NORMO') then
           if(exists('L1   ')) then
              kfl_normc_nsa = 1
           else if(exists('L-inf')) then
              kfl_normc_nsa = 0
           end if
        else if(words(1)=='MAXIM') then         ! Linearization iterations  TO BE COMPATIBLE WITH OTHER MODULES...
           miinn_nsa  = int(param(1))

        else if(words(1)=='LINEA') then         !linearization or outer iterations strategy, OLD WAY
           if(exists('RHS  ')) then
              kfl_linea_nsa = 0
           else if(exists('PICAR')) then
              kfl_linea_nsa = 1
              miinn_nsa     = getint('PICAR',   1_ip,'Number of Picard Linearization iterations')
              cotol_nsa     = getrea('TOLER',-1.0_rp,'Tolerance of Picard Linearization iterations _rp') 
           end if

        else if(words(1)=='SUBIT' .or. words(1)=='OUTER') then         !Sub- or outer iterations strategy, NEW WAY!!
           iauxi = 1
!           do while(words(1)/='ENDNE')
           do while (iauxi==1)
              call ecoute('nsa_reanut')
              if (words(1)=='ENDSU' .or. words(1)=='ENDOU') iauxi = 0    ! temporary, for back-compatibility
              if (exists('CONVE')) then
                 miinn_nsa= getint('ITERA',   1_ip ,'Maximum number of sub-iterations')
                 cotol_nsa= getrea('TOLER', -1.0_rp,'Tolerance for the sub-iterations')
                 corat_nsa= getrea('RATIO',  0.1_rp,'Tolerance for the sub-iterations')
                 if (exists('ADAPT')) kfl_adres_nsa = 1
              end if
              if (words(1)=='STRAT') then     
                 if (words(2)=='NEWRA') then     ! NEWTON-RAPHSON ()
                    kfl_linea_nsa = 2
                    kfl_delun_nsa = 1
                    minew_nsa= getint('START',   1_ip ,'Starting iteration for the Newton scheme')

!!                    if (exists('NODEL')) kfl_delun_nsa = 0   !only for debugging
!!                    if (words(3)=='INEXA') then  ! INEXACT NEWTON
!!                       kfl_linea_nsa = 2
!!                       if (words(4)=='RHSCO') then   ! INEXACT NEWTON WITH RHS CORRECTION
!!                          kfl_linea_nsa = 3
!!                          call runend('NSA_REANUT: RHS CORRECTION NOT PROGRAMMED YET')
!!                       end if
!!                    end if

                    if (exists('XVISC')) kfl_ximpl_nsa(1) = 1   
                    if (exists('XCONV')) kfl_ximpl_nsa(2) = 1                       
                    if (exists('XSTAB')) kfl_ximpl_nsa(3) = 1   
                    if (exists('XSHOC')) kfl_ximpl_nsa(4) = 1   
                 else if (words(2)=='JACOB') then  ! JACOBI or fix point 
                    kfl_linea_nsa = 1
                    kfl_delun_nsa = 1
                    if (exists('XVISC')) kfl_ximpl_nsa(1) = 1   
                    if (exists('XCONV')) kfl_ximpl_nsa(2) = 1                       
                    if (exists('XSTAB')) kfl_ximpl_nsa(3) = 1   
                    if (exists('XSHOC')) kfl_ximpl_nsa(4) = 1   
                 else if (words(2)=='FUMPR') then  ! full matrix point relaxation 
                    kfl_linea_nsa = 3
!                    kfl_diagi_nsa = 1
                    kfl_delun_nsa = 1
!                    kfl_timet_nsa = 1                                    ! Explicit 
                 end if
              else if(words(1)  =='PSEUD') then                        ! Tau-subiterations (i.e. dual or pseudo time stepping)
                 kfl_pseud_nsa = 1         
                 safet_pseud_nsa = getrea('SAFET',1.0e9_rp,'Safety factor for the pseudo-timestep')                                     
              end if

           end do
        else if(words(1)=='TIMET') then              ! Explicit (1) or Implicit (2) time treatment

           if(words(2)=='EXPLI') then
              ! Default for the explicit is forward euler
              kfl_timet_nsa=1
           else if(words(2)=='IMPLI') then
              ! Default for the implicit is backwards euler
              kfl_timet_nsa=2
              if(exists('BDF  ')) then
                 kfl_tisch_nsa=2
                 kfl_tiacc_nsa=2   ! BDF will be always 2nd order
                 call runend('NSA_REANUT: DO NOT USE BDF. USE CRANK-NICOLSON INSTEAD!')
              end if
              if(exists('CRANK')) then
                 kfl_tisch_nsa=3
                 kfl_tiacc_nsa=1   ! CN is 2nd order, but uses only one previous time step
              end if
              neule_nsa     = getint('EULER',0_ip,'#EULER TIME STEPS')
           end if

        else if( words(1) == 'PHYSI' ) then
           ivari = 1

        else if( words(1) == 'MATRI' ) then
           kfl_matri_nsa = 1

        else if(words(1)=='ALGEB') then
           !
           ! Solver
           !
           if (kfl_timet_nsa == 1) then              
              kfl_timet_nsa = 2   ! this is only meaningful for IMPLICIT, so implicit is forced when algeb is read
              kfl_tiacc_nsa = 1   ! the scheme is the default, i.e. backwards euler
           end if

           solve_sol => solve(ivari:)
           call reasol(1_ip)

        else if(words(1)=='ORTPA') then
           !
           ! Solver for the v_ort
           !
           solve_sol => solve(1:)
           call reasol(1_ip)

        else if(words(1)=='PRECO') then 
           !
           ! Preconditioner
           !
           solve_sol => solve(ivari:)
           call reasol(2_ip)
           if (solve(1)%kfl_preco == 14) kfl_diagi_nsa = 1

        else if(words(1)=='DAMPI') then
           kfl_rayle_nsa = 1
           kranr_nsa = 0
           vranr_nsa = 0.0_rp
           do while(words(1)/='ENDDA')
              call ecoute('nsa_reanut')
              if(words(1)=='RANGE') then
                 do while(words(1)/='ENDRA')
                    call ecoute('nsa_reanut')
                    if (words(1)=='COORD') then
                       if (words(2)=='XCOOR') irang= 1
                       if (words(2)=='YCOOR') irang= 2
                       if (words(2)=='ZCOOR') irang= 3
                    else if(words(1)=='LOWER') then                             
                       kranr_nsa(irang,1)   =  1
                       vranr_nsa(irang,1)   = param(1)    ! damping starts below this...
                       vranr_nsa(irang,3)   = param(2)    ! ... down to this
                    else if(words(1)=='GREAT') then
                       kranr_nsa(irang,2)   =  1       
                       vranr_nsa(irang,2)   =  param(1)   ! damping starts above this...
                       vranr_nsa(irang,3)   =  param(2)   ! ... up to this
                    else if(words(1)=='LIMIT') then
                       kranr_nsa(irang,1)   =  2
                       kranr_nsa(irang,2)   =  2
                       vranr_nsa(irang,1)   =  param(1)   ! lower limit
                       vranr_nsa(irang,2)   =  param(2)   ! upper limit
                       vranr_nsa(irang,3)   =  (param(1) + param(2)) / 2.0_rp   ! mean value
                    else if(words(1)=='ACOUS') then
                       kranr_nsa(irang,3)   =  1
                       frayl_nsa(irang,1)   =  param(1)   ! alpha / 2
                       frayl_nsa(irang,2)   =  param(2)   ! d
                    else if(words(1)=='GRAVI') then
                       kranr_nsa(irang,4)   =  1
                       frayl_nsa(irang,3)   =  param(1)
                       frayl_nsa(irang,4)   =  param(2)
                    end if
                 end do
              end if
           end do
    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------! 
        else if(words(1)=='NSCBC') then
           lodi_nsa       = 1
           prefe_lodi_nsa = getrea('PRESS',  -103400.0_rp, 'LODI, pressure parameter')
           sigma_lodi_nsa = getrea('SIGMA',     -0.25_rp, 'LODI, sigma parameter')
           
           euler_nsa = .false. 
           if(words(2)=='EULER') euler_nsa = .true. 
        end if
    !-----------------------------------------------------------------------!
 
     end do
     !-----------------------------------------------------------------------
     ! ADOC[0]> END_NUMERICAL_TREATMENT
     !-----------------------------------------------------------------------

     ! Check 
!     if (kfl_linea_nsa == 3) kfl_timet_nsa = 1 ! explicit

     !
     ! External timestep is fixed
     !
     if (dtime > 0._rp ) then
        dtext_nsa = dtime
        kfl_tiext_nsa = 1 
     end if

     iauxi= kfl_dttyp_nsa(1)+kfl_dttyp_nsa(2)+kfl_dttyp_nsa(3)
     if (iauxi > 0 .and. kfl_timet_nsa == 1) then
        ! explicit cases
        solve_sol(1) % kfl_preco  = SOL_LOCAL_DIAGONAL
     end if

     if (kfl_mod_elmop_nsa > 0 .and. kfl_timet_nsa == 1)  then
        ! for explicit elmoperations, always do this
        solve_sol(1) % kfl_preco  = SOL_LOCAL_DIAGONAL
     end if


  end if

end subroutine nsa_reanut
