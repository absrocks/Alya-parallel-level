!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_reanut.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Read numerical treatment
!> @details This routine reads the numerical treatment for solidz module
!> @}
!-----------------------------------------------------------------------

subroutine sld_reanut()

  use def_kintyp, only: ip, rp
  use def_inpout, only: words, param, exists
  use def_inpout, only: getint, getrea
  use def_master, only: solve, solve_sol, INOTSLAVE
  use def_domain, only: ndime
  use mod_ecoute, only: ecoute
  use def_solidz

  implicit none

  real(rp)     :: vauxi
  integer(ip)  :: icavi,i,iauxi

  if( INOTSLAVE ) then
     !
     !  Initializations (defaults)
     !
     kfl_stabi_sld = 0                          ! Stabilization
     kfl_xfeme_sld = 0                          ! XFEM family of enrichment strategies
     kfl_resid_sld = 0                          ! Residual by increment of displacement
     kfl_timet_sld = 1                          ! Time treatment (explicit=1, implicit=2)
     kfl_ninex_sld = 0                          ! Inexact newton (only for implicit cases)
     kfl_tisch_sld = 1                          ! Time integration scheme
     kfl_serei_sld = 0                          ! Selective reduced integration is off
     kfl_limit_sld = 0
     kfl_pseud_sld = 0
     kfl_penal_sld = 0
     kfl_volca_sld = 0
     kfl_prest_sld = 0                          ! Do not compute pre-stress
     kfl_gdepo     = 0                          ! Compute and store inverse of the deformation gradient
     kfl_safet_table_sld = 0
     minex_sld     = 2                          ! Newton inexact counter
     miinn_sld     = 1                          ! Max # of N-R iterations
     cotol_sld     = 1.0e-5_rp                  ! Convergence tolerance
     safet_sld     = 1.0_rp                     ! Safety factor for time step
     safex_sld     = 1.0_rp                     ! Safety factor for time step
     safma_sld     = 1.0_rp                     ! Safety factor for time step
     safet_pseud_sld     = 1.0_rp               ! Safety factor for pseudotime step (this value is changed when pseudo is used)
     safet_table_sld = 0.0_rp
     nisaf_sld     = 0
     factor_penal_sld  = 1.0_rp                 ! Penalty factor (this value is changed when penalization is used)
     dafac_sld     = 0.0_rp                     ! Damping factor for automatic stabilization of static problems
     sstol_sld     = 1.0e-5_rp                  ! Steady state tolerance
     masss_sld     = 1.0_rp                     ! Mass scaling factor
     mcavi_sld     = 1                          ! number of cavities (ventricles) - only to compute ejection rate
     kcavi_sld     = 0
     ocavi_sld     = 0.0_rp
     epsex_sld     = 1.0e-8_rp
     !
     ! Reach the section
     !
     call ecoute('sld_reanut')
     do while(words(1)/='NUMER')
        call ecoute('sld_reanut')
     end do
     !
     ! Begin to read data
     !
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> $ Numerical treatment definition
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> NUMERICAL_TREATMENT
     do while(words(1)/='ENDNU')
        call ecoute('sld_reanut')

        if(words(1)=='TIMET') then              ! Time treatment
           !
           ! ADOC[1]> TIME_TREATMENT:    IMPLICIT | EXPLICIT      $ Temporal integration scheme
           ! ADOC[d]> TIME_TREATMENT:
           ! ADOC[d]> Set the temporal integration method to IMPLICIT or EXPLICIT methods.
           !
           if(words(2)=='EXPLI') then
              kfl_timet_sld=1
           else if(words(2)=='IMPLI') then
              kfl_timet_sld=2
              if(exists('INEXA')) then
                 kfl_ninex_sld= 1
                 epsex_sld= getrea('EPSIL',1.0e-8_rp,'Perturbation to compute the secant stiffness matrix.')
              end if
           end if

        else if(words(1)=='TIMEI') then         ! Time integration scheme
           !
           ! ADOC[1]> TIME_INTEGRATION:  NEWMARK, DAMPED | UNDAMPED | CDIFF $ Integration scheme method
           ! ADOC[d]> TIME_INTEGRATION:
           ! ADOC[d]> Select the time integration scheme. The Newmark Beta-method, NEWMARK, is used by default.
           ! ADOC[d]> The option UNDAMPED uses the undamped trapezoidal rule (NEWBE=0.25, NEWGA=0.5, NEWAL=0.0).
           ! ADOC[d]> The option DAMPED uses the numerically damped integrator with damping proportional to NEWGA-0.5
           ! ADOC[d]> (NEWBE=0.65, NEWGA=0.9, NEWAL=0.0).
           ! ADOC[d]> The option CDIFF uses the explicit central differences method (NEWBE=0, NEWGA=0.5, NEWAL=0.0).
           ! ADOC[d]> The CDIFF option can only be used for EXPLICIT analysis. Moreover, the user can set manually
           ! ADOC[d]> NEWBE=<tt>real</tt>, NEWGA=<tt>real</tt> and NEWAL=<tt>real</tt>, when one of the previous
           ! ADOC[d]> options are missing.
           !
           if(exists('NEWMA')) then
              !
              ! Newmark-Beta Scheme
              !
              kfl_tisch_sld = 1_ip

              ! Defaults
              if (kfl_timet_sld == SLD_IMPLICIT_SCHEME .and. kfl_timei_sld == SLD_STATIC_PROBLEM ) then      ! Implicit
                 ! This parameters are not used in equilibrium (static) problems
                 tifac_sld(1)  = getrea('NEWBE',0.00_rp,'Beta factor')
                 tifac_sld(2)  = getrea('NEWGA',0.00_rp,'Gamma factor')
                 tifac_sld(3)  = getrea('NEWAL',0.00_rp,'Alpha factor')
              else if (kfl_timet_sld == SLD_IMPLICIT_SCHEME .and. kfl_timei_sld == SLD_DYNAMIC_PROBLEM ) then ! Implicit
                 ! UNDAM
                 tifac_sld(1)  = getrea('NEWBE',0.25_rp,'Beta factor')
                 tifac_sld(2)  = getrea('NEWGA',0.50_rp,'Gamma factor')
                 tifac_sld(3)  = getrea('NEWAL',0.00_rp,'Alpha factor')
                 ! Compatibility check (only when NEWBE is defined by the user)
                 if (tifac_sld(1) == 0.00_rp) call runend("SLD_REANUT: IMPLICIT SCHEMES ARE MEANINGLESS WHEN BETA=0")
              else                                                              ! Explicit
                 ! CDIFF
                 tifac_sld(1)  = getrea('NEWBE',0.00_rp,'Beta factor')
                 tifac_sld(2)  = getrea('NEWGA',0.50_rp,'Gamma factor')
                 tifac_sld(3)  = getrea('NEWAL',0.00_rp,'Alpha factor')
              end if

              if (exists('DAMPE')) then
                 tifac_sld(1)  = 0.65_rp
                 tifac_sld(2)  = 0.90_rp
                 tifac_sld(3)  = 0.00_rp
              else if (exists('UNDAM')) then
                 tifac_sld(1)  = 0.25_rp
                 tifac_sld(2)  = 0.50_rp
                 tifac_sld(3)  = 0.00_rp
              else if (exists('CDIFF')) then
                 tifac_sld(1)  = 0.00_rp
                 tifac_sld(2)  = 0.50_rp
                 tifac_sld(3)  = 0.00_rp
              end if

           else if( exists('TCHAM') ) then
              !
              ! Dissipative Tchamwa - Wielgosz scheme
              !
              kfl_tisch_sld = 2_ip
              tifac_sld(4)  = getrea('PHITW',1.033_rp,'Phi factor')

           else if( exists('RUNGE') ) then
              !
              ! Runge-Kutta Scheme
              !
              kfl_tisch_sld = 3_ip

           else

              call runend('SLD_REANUT: TIME INTEGRATION SCHEME NOT DEFINED')

           end if

        else if(words(1)=='SAFET') then         ! Safety factor
           !
           ! ADOC[1]> SAFETY_FACTOR:     real                     $ Temporal Safety factor
           ! ADOC[d]> SAFETY_FACTOR: Set this parameter equal to the safety factor. This parameter is
           ! ADOC[d]> applied to the critical time step. In Explicit analysis is recommended a value smaller
           ! ADOC[d]> than one, while in Implicit analysis greater than one. By default, this parameter is set to 1.
           !
           safet_sld = param(1)
           !
           safex_sld  = getrea('EXPON',1.0_rp,'#Safety Factor Exponential time function')
           safma_sld  = getrea('MAXIM',1.0e9_rp,'#Maximum safety factor function')
           nisaf_sld  = getint('INITI',0_ip    ,'#Initial step for variable CFL')
           if (exists('TABLE')) then
              iauxi= 0
              call ecoute('sld_reanut')
              do while(words(1)/='ENDTA')
                 iauxi= iauxi+1
                 safet_table_sld(1,iauxi)  = getrea('VALUE',1.0_rp,'#Safety Factor')
                 safet_table_sld(2,iauxi)  = getrea('RESID',0.0_rp,'#Residual')
                 if (iauxi == 10) &
                      call runend('SLD_REANUT: SAFET TABLE MUST BE SMALLER THAN 10.')
                 call ecoute('sld_reanut')
              end do
              safet_sld= safet_table_sld(1,1)  ! starting safety factor value
              kfl_safet_table_sld = iauxi
           end if

        else if ( words(1) == 'STABI' ) then
           !
           ! ADOC[1]> STABILIZATION: real                         $ Damping factor
           ! ADOC[d]> STABILIZATION: Automatic stabilization with constant damping for static problems.
           !
           kfl_stabi_sld = 1
           dafac_sld = param(1)

        else if(words(1)=='SUBIT') then         !Sub-iterations strategy, NEW WAY!!

           do while (words(1)/='ENDSU')
              call ecoute('sld_reanut')
              if(words(1)  =='PSEUD') then   ! Tau-subiterations (i.e. dual or pseudo time stepping)
                 kfl_pseud_sld = 1
!                 safet_pseud_sld = getrea('SAFET',0.1_rp,'Safety factor for the pseudo-timestep')
                 vauxi = getrea('FACTOR',100.0_rp,'Factor for the pseudo-timestep')
                 safet_pseud_sld = sqrt(1.0_rp / vauxi / vauxi)
              else if(words(1)  =='PENAL') then                        ! Penalization
                 kfl_penal_sld = 1
                 factor_penal_sld = getrea('FACTOR',10.0_rp,'Penalty factor')
              else if(words(1)=='MAXIM') then         ! Linearization iterations
                 miinn_sld  = int(param(1))
              else if(words(1)=='CONVE') then         ! Convergence tolerance
                 cotol_sld = param(1)
              else if(words(1)=='RESID') then         ! Convergence criterion
                 if(exists('DISPL')) then
                    kfl_resid_sld = 0
                 else if(exists('FORCE')) then
                    kfl_resid_sld = 1
                 else if(exists('TOTAL')) then
                    kfl_resid_sld = 2
                 end if
              end if
           end do

        else if(words(1)=='SELEC') then         ! Selective reduced integration enabled
           kfl_serei_sld = 1

        else if(words(1)=='STEAD') then         ! Steady state tolerance
           !
           ! ADOC[1]> STEADY_STATE_TOL:  real                     $ Tolerance for steady-state
           ! ADOC[d]> STEADY_STATE_TOL: Tolerance for detection of steady state in transient problems.
           ! ADOC[d]> When the time residual calculated between time steps achieves the
           ! ADOC[d]> steady state tolerance the analysis stop.
           !
           sstol_sld = param(1)

        else if(words(1)=='NOSTE') then         ! Force no steady state check
           !
           ! ADOC[1]> NO_STEADY_STATE                             $ Steady state deactivated
           ! ADOC[d]> NO_STEADY_STATE: Steady state deactivated
           !
           sstol_sld = -1.0_rp

        else if(words(1)=='MASSS') then         ! Mass scaling
           !
           ! ADOC[1]> MASS_SCALING:   real                        $ Mass scaling
           ! ADOC[d]> MASS_SCALING: Mass scaling factor for solving quasi-static simulations in Explicit analysis.
           ! ADOC[d]> Artificially increaseing the material density, <tt>rho</tt>, by a factor <tt>f^2</tt>
           ! ADOC[d]> reduces <tt>n</tt> to <tt>n/f</tt>, just like decreasing <tt>T</tt> to <tt>T/f</tt>.
           ! ADOC[d]> This concept reduces the ratio of the event time to time for wave propagation across an
           ! ADOC[d]> element while leaving the event time fixed, which allows rate-dependent behavior to be
           ! ADOC[d]> included in the analysis. Mass scaling has exactly the same effect on inertia forces as
           ! ADOC[d]> speeding up the time of the simulations.
           !
           masss_sld = param(1)
           densi_sld(1,1:nmate_sld) = masss_sld*densi_sld(1,1:nmate_sld)

        else if(words(1)=='RESID') then         ! Convergence criterion
           !
           ! ADOC[1]> RESIDUAL:          DISPL | FORCE | ENERG | TOTAL    $ Convergence criteria (Implicit analysis)
           ! ADOC[d]> RESIDUAL: This option is used to choose the convergence criteria in the Newton-Raphson iterations.
           ! ADOC[d]> When DISPL is used the convergence criteria to terminate the iterations is a criterion based on
           ! ADOC[d]> the magnitude of the displacement increments. Set to FORCE in order to use a criterion based on
           ! ADOC[d]> the magnitude of the residual of forces. Set to ENERGY in order to use an energy error criterion.
           ! ADOC[d]> Set TOTAL, the three criteria have to be fulfilled with the same tolerance. This option is only
           ! ADOC[d]> applicable to Implicit analysis. By default, the displacement increment error, DISPL, is used.
           !
           if(exists('DISPL')) then
              kfl_resid_sld = 0_ip
           else if(exists('FORCE')) then
              kfl_resid_sld = 1_ip
           else if(exists('ENERG')) then
              kfl_resid_sld = 2_ip
           else if(exists('TOTAL')) then
              kfl_resid_sld = 3_ip
           end if

        else if(words(1)=='MAXIM') then         ! Maximum N-R iterations
           !
           ! ADOC[1]> MAXIMUM_ITERATION: int                      $ Max no. of N-R iterations (Implicit analysis)
           ! ADOC[d]> MAXIMUM_ITERAITON: Maximum number of N-R iterations. This option is only applicable to
           ! ADOC[d]> Implicit analysis. By default is set to 1.
           !
           miinn_sld  = int(param(1))

        else if(words(1)=='CONVE') then         ! Convergence tolerance
           !
           ! ADOC[1]> CONVERGENCE_TOL:   real                     $ Convergence tolerance (Implicit analysis)
           ! ADOC[d]> CONVERGENCE_TOL: Convergence tolerance of the Newton-Raphson iterations. This tolerance is
           ! ADOC[d]> the same for all the convergence criteria when TOTAL is used. By default is set to 1e-5.
           !
           cotol_sld = param(1)

        else if(words(1)=='ALGEB') then
           !
           ! ADOC[1]> ALGEBRAIC_SOLVER
           ! ADOC[2]>   SOLVER:         DIRECT | GMRES | CG
           ! ADOC[2]>   CONVERGENCE:    ITERA=int, TOL=real
           ! ADOC[2]>   PRECONDITIONER: DIAGONAL | RAS
           ! ADOC[2]>   RESIDUAL:       RHS
           ! ADOC[2]>   OPTIONS:        ZERO_FIXITY
           ! ADOC[2]>   OUTPUT:         CONVERGENCE
           ! ADOC[1]> END_ALGEBRAIC_SOLVER
           ! ADOC[d]> ALGEBRAIC_SOLVER: Definition of the algebraic solver. For more details, see Alya
           ! ADOC[d]> documentation for Solvers.
           call reasol(1_ip)

        else if(words(1)=='PRECO') then
           call reasol(2_ip)

        else if(words(1)=='ROTAT')then
           !
           ! ADOC[1]> ROTATION:          On | Off                 $ Rotation to ref. configuration (push forward)
           ! ADOC[d]> ROTATION: Rotation to reference configuration (push forward)
           !
           if(words(2)=='ON') kfl_gdepo = 1_ip

        else if(words(1)=='ENRIC') then         ! Enrichement strategy: XFEM/GFEM
           kfl_xfeme_sld= 1
           kfl_xfcra_sld= 1                     ! Local crack definition (as opposed to level set definition)
           if(exists('X-FEM')) kfl_xfeme_sld= 1     ! Default value
           if(exists('LOCAL')) kfl_xfcra_sld= 1     ! Local crack definition
           if(exists('LEVEL')) kfl_xfcra_sld= 2     ! Level set crack definition

        else if(words(1)=='SHELL') then
           do while( words(1) /= 'ENDSH' )
              if(words(1)=='ALGEB') then
                 solve_sol => solve(2:)
                 call reasol(1_ip)
              end if
              call ecoute('SLD_REANUT')
           end do
           solve_sol => solve(1:)

        else if(words(1)=='LIMIT') then         ! Limiter
           kfl_limit_sld = 1_ip

        else if(words(1)=='PREST') then         ! Compute pre-stress
           !
           ! ADO[1]> PRESTRESS                                   $ Prestress initial condition
           ! ADO[d]> PRESTRESS                                   $ Prestress initial condition
           ! ADO[d]> Compute the pre-stress field
           !
           kfl_prest_sld = 1

        else if(words(1)=='CAVIT') then         ! TO BE DEPRECATED

           call runend('SLD_REANUT: ESTO CAMBIO!!! PREGUNTARLE A MARIANO')

           kfl_volca_sld = 1_ip
           mcavi_sld = getint('TOTAL',1_ip,'Number of cavities to compute the volume')
           call ecoute('sld_reanut')
           do while(words(1)/='ENDCA')
              !              if (words(1)=='VENTR') then
              !                 mcavi_sld= int(param(1))               ! Number of cavities
              if (words(1)=='ORIGI') then
                 do icavi=1,mcavi_sld
                    i= (icavi-1) * ndime + 1
                    ocavi_sld(1:ndime,icavi)= param(i:i+ndime-1)
                 end do
              else if (words(1)=='BOUND') then          ! Boundary set that defines de inner surface
                 do icavi=1,mcavi_sld
                    iocav_sld(icavi)= int(param(icavi))
                 end do
              else if (words(1)=='MESHN') then          ! MESH_NODE
                 call runend('SLD_REANUT: MESH NODE FOR CAVITIES IS DEPRECATED. USE ORIGIN.')
              end if

              call ecoute('sld_reanut')
           end do

        end if

     end do
     ! ADOC[0]> END_NUMERICAL_TREATMENT

     !
     ! Final compatibility corrections
     !
     if (kfl_timet_sld==1 .and. kfl_pseud_sld == 0) then
        miinn_sld     = 1
     end if

  end if

end subroutine sld_reanut
