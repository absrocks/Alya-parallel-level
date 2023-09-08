!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_reanut.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Read numerical treatment
!> @details Read numerical treatment
!> @} 
!------------------------------------------------------------------------
subroutine por_reanut
  use def_parame
  use def_inpout
  use def_master
  use def_solver
  use def_porous
  use def_domain
  implicit none

  integer(ip) :: ipres,isatu


  if(kfl_paral<=0) then
     ! 
     !  Initializations (defaults)
     !
     kfl_ellen_por = 0                                ! Minimum element length - for the moment unique
     do kprsa_por=1,nprsa_por
        kfl_sgsti_por(kprsa_por) = 0                      ! Subgrid scale time tracking
        kfl_sgsno_por(kprsa_por) = 0                      ! Subgrid scale non-lniear tracking
        kfl_taust_por(kprsa_por) = 1                      ! Tau calculation option
        kfl_ortho_por(kprsa_por) = 0                      ! Orthogonal SGS
        kfl_limit_por(kprsa_por) = 0                      ! Lmiter
        kfl_shock_por(kprsa_por) = 0                      ! Shock capturing off
        kfl_tiacc_por(kprsa_por) = 1                      ! First order time integ.
     end do
     kfl_difwe_por = 0                                    ! No diffusion to be added at the wells

     miinn_por     = 1                                ! One internal iteration

     solve_sol     => solve                           ! Solver type
     ! 
     !  Initializations (defaults)
     !
     kfl_ellen_por = 0                                ! Minimum element length
     kfl_sgsti_por = 0                                ! Subgrid scale time tracking
     kfl_sgsno_por = 0                                ! Subgrid scale non-lniear tracking
     kfl_taust_por = 1                                ! Tau calculation option
     kfl_ortho_por = 0                                ! Orthogonal SGS
     kfl_limit_por = 0                                ! Lmiter
     neule_por     = 0                                ! Number of Euler time steps
     kfl_tisch_por = 1                                ! Trapezoidal rule
     kfl_normc_por = 2                                ! L2 norm for convergence
     miinn_por     = 1                                ! One internal iteration
!     misgs_por     = 1                                ! Max # of SGS iterations

     staco_por(1)  = 1.0_rp                           ! Diffusive term
     staco_por(2)  = 1.0_rp                           ! Convective term
     staco_por(3)  = 1.0_rp                           ! Reactive term
     shock_por     = 0.0_rp                           ! SC parameter
     difwe_por     = 0.0_rp                           ! No diffusion to be added at the wells
     safet_por     = 1.0e10_rp                        ! Safety factor
     sstol_por     = 1.0e-5_rp                        ! Steady-state tolerance
     cotol_por     = 0.1_rp                           ! Internal tolerance
     relax_por     = 1.0_rp                           ! Relaxation factor
     bemol_por     = 0.0_rp                           ! Bemol (convectice term)
     relsg_por     = 1.0_rp                           ! Relaxation parameter of subgrid scale
     tosgs_por     = 0.01_rp                          ! Subgrid scale tolerance
     !
     ! Reach the section
     !
     call ecoute('por_reanut')
     do while(words(1)/='NUMER')
        call ecoute('por_reanut')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDNU')
        call ecoute('por_reanut')

        if(words(1)=='ELEME') then
           call realen(kfl_ellen_por)

        else if(words(1)=='MAXIM') then
           miinn_por = int(param(1))

        else if(words(1)=='CONVE') then
           cotol_por = param(1)

!        else if(words(1)=='RELAX') then
!           relax_por = param(1)

        else if( words(1) == 'DIFWE' ) then
           !
           ! Add diffusion at the wells
           !
           kfl_difwe_por     = 1_ip
           if(exists('VALUE')) &
                difwe_por     = getrea('VALUE',0.0_rp,'#Diffusion to be added at the wells')

        else if( words(1) == 'PRESS' ) then
           !
           ! Current problem is pressure equation (used to define solver)
           !
           kprsa_por = 1
           ipres     = 1

        else if( words(1) == 'ENDPR' ) then
           !
           ! Current problem is no longer pressure equation (used to define solver)
           !
           ipres     = 0

        else if( words(1) == 'SATUR' ) then
           !
           ! Current problem is saturation equation (used to define solver)
           !
           kprsa_por = 2
           isatu     = 1

        else if( words(1) == 'ENDSA' ) then
           !
           ! Current problem is no longer saturation equation (used to define solver)
           !
           isatu     = 0

        else if(words(1)=='SHOCK') then
           if(exists('ISOTR').or.exists('ON   ')) then
              kfl_shock_por(kprsa_por) = 1
              if(exists('VALUE')) &
                   shock_por     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')  ! beware for the moment it is the same for pr & sat
           else if(exists('ANISO')) then
              kfl_shock_por(kprsa_por) = 2
              if(exists('VALUE')) &
                   shock_por     = getrea('VALUE',0.0_rp,'#Shock capturing parameter')
           end if

        else if( words(1) == 'ALGEB' ) then
           !
           ! ADOC[1]> XXX_EQUATION                                                                              $ Define solver for equation XXX
           ! ADOC[2]>   ALGEBRAIC_SOLVER
           ! ADOC[2]>     INCLUDE ./sample-solver.dat
           ! ADOC[2]>   END_ALGEBRAIC_SOLVER
           ! ADOC[1]> END_XXX_EQUATION
           ! ADOC[d]> XXX_EQUATION:
           ! ADOC[d]> Define the algebraic solver options for equation XXX.
           ! ADOC[d]> XXX: PRESSURE | SATURATION
           !  
           solve_sol => solve(kprsa_por:)
           call reasol(1_ip)

        else if( words(1) == 'PRECO' ) then  
           !
           ! See previous option
           !
           solve_sol => solve(kprsa_por:)
           call reasol(2_ip)

        end if
     end do
  end if

end subroutine por_reanut
