!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_solsch.f90
!> @author  Guillaume Houzeaux
!> @date    20/09/2016
!> @brief   Fractional step
!> @details Fractional step
!> @} 
!-----------------------------------------------------------------------
module mod_nsi_multi_step_fs

  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_nastin
  use def_kermod,                      only : kfl_noslw_ker
  use def_kermod,                      only : gasco
  use mod_ker_proper,                  only : ker_proper
  use mod_memory,                      only : memory_alloca
  use mod_memory,                      only : memory_deallo
  use mod_memory,                      only : memory_size
  use mod_nsi_schur_operations
  use mod_matrices,                    only : matrices_gradient_divergence
  use mod_matrices,                    only : matrices_laplacian
  use mod_nsi_dirichlet_global_system, only : nsi_dirichlet_matrix_rhs
  use mod_solver,                      only : solver_lumped_mass_system
  use mod_messages,                    only : livinf
  use mod_solver,                      only : solver_solve
  use mod_solver,                      only : solver_exchange
  use mod_solver,                      only : solver_explicit
  use mod_local_basis,                 only : local_basis_global_to_local
  use mod_local_basis,                 only : local_basis_local_to_global
  implicit none
  private

  real(rp)               :: gamma1

  real(rp),   pointer    :: vv(:)       ! auxiliar vector
  real(rp),   pointer    :: xx(:)       ! auxiliar vector
  real(rp),   pointer    :: aux(:)      ! auxiliar vector
  real(rp),   pointer    :: uu0(:)      ! auxiliar vector
  real(rp),   pointer    :: rm(:,:)     ! residual of momentum at time 
  real(rp),   pointer    :: rc(:)       ! residual of continuity 
  real(rp),   pointer    :: divup(:)    ! div of up 
  real(rp),   pointer    :: rt(:)       ! Total residual of momentum 
  real(rp),   pointer    :: rt_tmp(:)   ! Total residual of momentum
  
  real(rp),   pointer    :: bu(:)         
  real(rp),   pointer    :: bp(:)
  real(rp),   pointer    :: uu(:)
  real(rp),   pointer    :: pp(:)
  
  real(rp),   pointer    :: Q(:)
  real(rp),   pointer    :: Dp(:)
  real(rp),   pointer    :: Aupp(:)
  
  real(rp),   pointer    :: Tim(:)      ! Nodal projection of dt/rho
  real(rp),   pointer    :: Mass(:)     ! Mass matrix
  real(rp),   pointer    :: Grad(:,:)   ! Gradient matrix   G_ij=-\int (div Ni) Nj
  real(rp),   pointer    :: Div(:,:)    ! Divergence matrix D_ij= \int Ni (div Nj) => D=-G^t
  real(rp),   pointer    :: Lapl(:)     ! Laplacian matrix

  real(rp)               :: c3 
  real(rp)               :: a(5)
  real(rp)               :: b(5,5)
  real(rp)               :: sigma(5)
  integer(ip), parameter :: conservative = 0
  integer(ip), parameter :: n_moin_steps = 4   ! Number of steps as Moin would do it. That is, without the enhanced extrapolation proposed by Capuano for the case kfl_fscon_nsi==0 

  real(rp),   pointer    :: Nua(:)      ! Nodal projection of mu/rho
  real(rp),   pointer    :: Nul(:)      ! Nodal projection of M rho not updated
  real(rp),   pointer    :: Nu0l(:)     ! Nodal projection of M rho not updated

  real(rp)               :: rho(1)
  real(rp)               :: dt_inv0
  real(rp)               :: dt_inv00
  real(rp)               :: dt_inv000
  integer(ip)            :: time_iter
  integer(ip)            :: dummi
 
  type(soltyp), pointer  :: solve_mome(:)
  type(soltyp), pointer  :: solve_mass(:)
  
  public :: nsi_multi_step_fs_solution
  public :: nsi_multi_step_fs_initialization
  public :: nsi_multi_step_fs_memory
  public :: nsi_multi_step_fs_matrices

  public :: Grad
  
contains 

  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Initializaation 
  !> @details Initializaation
  !
  !-----------------------------------------------------------------------

  subroutine nsi_multi_step_fs_initialization()

    integer(ip)       :: idofn,ipoin,idime
    real(rp)          :: auxra, epsle = 0.05_rp
    
    nullify(vv)
    nullify(xx)
    nullify(aux)
    nullify(uu0)
    nullify(rm)
    nullify(rc)
    nullify(divup)
    nullify(rt)
    nullify(rt_tmp)

    nullify(bu)
    nullify(bp)
    nullify(uu)
    nullify(pp)

    nullify(Q)
    nullify(Dp)
    nullify(Aupp)
        
    nullify(Tim)
    nullify(Mass)
    nullify(Grad)
    nullify(Div)
    nullify(Lapl)

    nullify(Nua)    
    nullify(Nul)    
    nullify(Nu0l)   

    time_iter = 1

  end subroutine nsi_multi_step_fs_initialization
  
  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Compute some matrices
  !> @details Compute some matrices
  !
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_matrices()
    
    if( kfl_grad_div_nsi /= 0 ) then
       solve_sol => solve(2:)
       call matrices_gradient_divergence(Grad,Div,DEALLOCATE_MATRICES=.true.)
       call matrices_laplacian          (Lapl,DEALLOCATE_MATRIX=.true.)
       call nsi_dirichlet_matrix_rhs    (Q=Lapl,ROTATION=.false.,DIRICHLET=.true.)
       call nsi_dirichlet_matrix_rhs    (Aup=Grad,Apu=Div,ROTATION=.true.,DIRICHLET=.false.)
       Q => Lapl
    else
       Q => lapla_nsi
    end if
    
    bu     => rhsid
    bp     => rhsid(ndbgs_nsi+1:)
    uu     => unkno
    pp     => unkno(ndbgs_nsi+1:)   
    Tim    => dt_rho_nsi
    if( INOTEMPTY ) then
       Nul    => mass_rho_nsi(:,1)
       Nu0l   => mass_rho_nsi(:,2)
    else
       nullify(Nul)
       nullify(Nu0l)
    end if
    Mass   => vmass
  
  end subroutine nsi_multi_step_fs_matrices

  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Allocate memory
  !> @details Allocate memory and define module constants
  !
  !-----------------------------------------------------------------------

  subroutine nsi_multi_step_fs_memory()

    integer(ip) :: ipoin
    
    gamma1 =  1.0_rp - gamma_nsi

    call memory_alloca(mem_modul(1:2,modul),'VV'    ,'mod_nsi_multi_step_fs',vv,    max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'XX'    ,'mod_nsi_multi_step_fs',xx,    max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'UU0'   ,'mod_nsi_multi_step_fs',uu0,   max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'RM'    ,'mod_nsi_multi_step_fs',rm,    max(1_ip,ndime*npoin),5_ip) 
    call memory_alloca(mem_modul(1:2,modul),'RT'    ,'mod_nsi_multi_step_fs',rt,    max(1_ip,ndime*npoin)) 
    call memory_alloca(mem_modul(1:2,modul),'RTTMP' ,'mod_nsi_multi_step_fs',rt_tmp,max(1_ip,ndime*npoin)) 
    call memory_alloca(mem_modul(1:2,modul),'RC'    ,'mod_nsi_multi_step_fs',rc,    max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'AUX'   ,'mod_nsi_multi_step_fs',aux,   max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'DIVUP' ,'mod_nsi_multi_step_fs',divup, max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'DP'    ,'mod_nsi_multi_step_fs',Dp,    max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'AUPP'  ,'mod_nsi_multi_step_fs',Aupp,  max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'NUA'   ,'mod_nsi_multi_step_fs',Nua,   max(1_ip,npoin))

    solve_mome => solve(NSI_SOLVER_MOMENTUM:)
    solve_mass => solve(NSI_SOLVER_CONSISTENT_MASS:)
    
    do ipoin = 1,npoin
       Nua(ipoin) = 1.0_rp
    end do
    
    a     = 0.0_rp
    b     = 0.0_rp
    sigma = 0.0_rp

    if(kfl_tiacc_nsi == 2) then !Heun’s method or RK2 just for validation

       a(4)     = 1.0_rp
       
       b(4,3)   = 1.0_rp

       sigma(3) = 0.5_rp
       sigma(4) = 0.5_rp

    else if(kfl_tiacc_nsi == 4) then ! standard 4 order RK

       a(2)     = 0.5_rp
       a(3)     = 0.5_rp
       a(4)     = 1.0_rp

       b(2,1)   = 0.5_rp 
       b(3,2)   = 0.5_rp
       b(4,3)   = 1.0_rp

       sigma(1) = 1.0_rp/6.0_rp 
       sigma(2) = 1.0_rp/3.0_rp 
       sigma(3) = 1.0_rp/3.0_rp 
       sigma(4) = 1.0_rp/6.0_rp 

    else ! energy preserving 3 order RK 

       if(conservative == 1) then 

          c3       = 1.0_rp/4.0_rp

          a(2)     = (c3 - 1.0_rp)/(4.0_rp*c3 - 3.0_rp)
          a(3)     = c3
          a(4)     = 1.0_rp

          b(2,1)   = (c3-1.0_rp)/(4.0_rp*c3-3.0_rp) 
          b(3,1)   = c3 - ((2_rp*c3-1.0_rp)*(4.0_rp*c3-3.0_rp))/(2.0_rp*(c3-1.0_rp)) 
          b(4,1)   = -1.0_rp * ((2.0_rp*c3-1.0_rp)**2/(2.0_rp*(c3-1.0_rp)*(4.0_rp*c3-3.0_rp)))  

          b(3,2)   = ((2_rp*c3-1.0_rp)*(4.0_rp*c3-3.0_rp))/(2.0_rp*(c3-1.0_rp)) 
          b(4,2)   = (6.0_rp*c3**2-8.0_rp*c3+3.0_rp)/(2_rp*(c3-1.0_rp)*(2_rp*c3-1.0_rp))

          b(4,3)   = (c3-1.0_rp)/((2*c3-1.0_rp)*(4.0_rp*c3-3.0_rp))

          sigma(1) = -(1.0_rp)/(12.0_rp*(c3-1.0_rp)) 
          sigma(2) = ((4.0_rp*c3-3.0_rp)**2)/(12.0_rp*(c3-1.0_rp)*(2.0_rp*c3-1.0_rp)) 
          sigma(3) = (-1.0_rp)/(12.0_rp*(c3-1.0_rp)*(2.0_rp*c3-1.0_rp))
          sigma(4) = (4.0_rp*c3-3.0_rp)/(12.0_rp*(c3-1.0_rp))

       else

          !  a(3) = 0.5_rp
          !  a(4) = 1.0_rp  

          !  b(3,2) = 0.5_rp
          !  b(4,2) = -1.0_rp
          !  b(4,3) = 2.0_rp

          !  sigma(2) = 1.0_rp/6.0_rp 
          !  sigma(3) = 2.0_rp/3.0_rp
          !  sigma(4) = 1.0_rp/6.0_rp

          a(3)     = 1.0_rp
          a(4)     = 0.5_rp  

          b(3,2)   = 1_rp
          b(4,2)   = 1.0_rp/4.0_rp
          b(4,3)   = 1.0_rp/4.0_rp

          sigma(2) = 1.0_rp/6.0_rp 
          sigma(3) = 1.0_rp/6.0_rp
          sigma(4) = 2.0_rp/3.0_rp

          ! a(3) = 8.0_rp/15.0_rp
          ! a(4) = 2.0_rp/3.0_rp  

          ! b(3,2) = 8.0_rp/15.0_rp
          ! b(4,2) = 1.0_rp/4.0_rp
          ! b(4,3) = 5.0_rp/12.0_rp

          ! sigma(2) = 1.0_rp/4.0_rp 
          ! sigma(3) = 0.0_rp
          ! sigma(4) = 3.0_rp/4.0_rp

       end if

    end if

  end subroutine nsi_multi_step_fs_memory

  !-----------------------------------------------------------------------
  !> 
  !> @author  lehmkuhl and houzeaux
  !> @date    2020-03-21
  !> @brief   Multi-step FS
  !> @details Multistep fractional step method from:
  !>          Explicit Runge–Kutta schemes for incompressible
  !>          flow with improved energy-conservation properties
  !>          F.Capuano, G.Coppola, L.Rández, L.deLuca, Journal of
  !>          ComputationalPhysics 328: 86–94 (2017)
  !>
  !>          Schemme 3p5q(4) from the paper 3o time, 5o energy and 4 steps
  !
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_solution()
    
    integer(ip) :: idofn,ipoin,idime
        
     do ipoin = 1,npoin
        do idime = 1,ndime
           idofn = (ipoin-1)*ndime + idime
           uu0(idofn)   = uu(idofn)
        end do
     end do
     if( kfl_regim_nsi == 3 ) then
       do ipoin = 1,npoin
          Nu0l(ipoin) = Nul(ipoin)
        end do
     end if

     if(kfl_tiacc_nsi == 2) then
        !
        ! Order 2
        !
        call nsi_multi_step_fs_eval(3_ip)
        call nsi_multi_step_fs_eval(4_ip)
        
     else if(kfl_tiacc_nsi == 3) then  
        !
        ! Order 3
        !
        if (conservative == 1) then 
           call nsi_multi_step_fs_eval(1_ip)
           call nsi_multi_step_fs_eval(2_ip)
           call nsi_multi_step_fs_eval(3_ip)
           call nsi_multi_step_fs_eval(4_ip)
        else 
           call nsi_multi_step_fs_eval(2_ip)
           call nsi_multi_step_fs_eval(3_ip)
           call nsi_multi_step_fs_eval(4_ip)           
        end if
        
     else 
        !
        ! Order 4
        !
        call nsi_multi_step_fs_eval(1_ip)
        call nsi_multi_step_fs_eval(2_ip)
        call nsi_multi_step_fs_eval(3_ip)
        call nsi_multi_step_fs_eval(4_ip)
        
     end if

     dt_inv000 = dt_inv00 
     dt_inv00  = dt_inv0 
     dt_inv0   = 1.0_rp/dtinv_nsi
     time_iter = time_iter + 1

  end subroutine nsi_multi_step_fs_solution

  !-----------------------------------------------------------------------
  !> 
  !> @author  lemhkuhl and houzeaux
  !> @date    2020-03-20
  !> @brief   Assemble equations
  !> @details Assemble equations
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_eval(istep)

    integer(ip) , intent(in) :: istep
    integer(ip)              :: idofn,ipoin,idime
    real(rp)                 :: timea,timeb

    !-----------------------------------------------------------------
    !
    ! Go back to local system
    !
    !-----------------------------------------------------------------
    
    call nsi_solini(3_ip)                                     ! Initialize solver

    call local_basis_global_to_local(kfl_fixrs_nsi,uu)        ! Global to local
    call local_basis_global_to_local(kfl_fixrs_nsi,uu0)       ! Global to local

    !-----------------------------------------------------------------
    !
    ! Assemble equations
    !
    !-----------------------------------------------------------------

    call nsi_matrix()                                           

    !-----------------------------------------------------------------
    !
    ! Compute density
    !
    !-----------------------------------------------------------------

    if( kfl_regim_nsi == 3 .and. kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0 ) then
       !
       ! Low Mach
       !
       do ipoin = 1,npoin
          if (kfl_lookg_nsi > 0) then
             call ker_proper('DENSI','IPOIN',ipoin,dummi,rho)
             Nua(ipoin)  = rho(1)
          else
             Nua(ipoin)  = (prthe(1)/gasco) * (wmean(ipoin,1)/tempe(ipoin,1))
          endif
       end do
       !Nua = 0.5*(3.0_rp*Nu-Nu0)

    else if( kfl_regim_nsi == 3 .and. (time_iter /= 1) .and. kfl_coupl(ID_NASTIN,ID_CHEMIC) == 0 ) then
       
       do ipoin = 1,npoin
          Nua(ipoin)  = (prthe(1)/gasco) * (1.0_rp/tempe(ipoin,1))
       end do

    else if( kfl_surte_nsi == 2_ip ) then

       Nua = 1.0_rp

    else
       !
       ! Incompressible
       !
       do ipoin = 1,npoin
          Nua(ipoin) = mass_rho_nsi(ipoin,1)
       end do
       call solver_lumped_mass_system(1_ip,Nua,EXCHANGE=.false.)
    
    end if

    call livinf(56_ip,' ',modul)    
    call livinf(160_ip,' ',1_ip)  
    call cputim(timea)

    if(istep == 4_ip) then
       call nsi_multi_step_fs_solution_final()
    else
       call nsi_multi_step_fs_solution_sj(istep)
    endif

    !-----------------------------------------------------------------
    !
    ! Go back to global system and actualize velocity and pressure,
    ! required for next assembly
    !
    !-----------------------------------------------------------------

    call livinf(165_ip,'C',0_ip)

    call local_basis_local_to_global(kfl_fixrs_nsi,uu)
    call local_basis_local_to_global(kfl_fixrs_nsi,uu0)

    call nsi_updunk(ITASK_INNITE)
 
    call livinf(164_ip,' ',1_ip) 
    call cputim(timeb)
    cpu_ass_sol_nsi(2) = timeb - timea

  end subroutine nsi_multi_step_fs_eval  

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Momentum
  !> @details Intermediate solution of the momentum equations
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_solution_sj(istep)

    integer(ip) , intent(in) :: istep
    integer(ip)              :: idofn,ipoin,idime,jstep
    real(rp)                 :: aux2 
    real(rp)                 :: time1,time2

    call cputim(time1)
    !
    ! Momentum residual rm^k
    !
    if( INOTEMPTY ) then
       if( kfl_fscon_nsi == 0 ) then 
          if( ittim > n_moin_steps ) then 
             do ipoin = 1,npoin
                do idime = 1,ndime
                   idofn = (ipoin-1)*ndime + idime
                   rm(idofn,istep) = bu(idofn)
                end do
                aux(ipoin) = 0.5_rp*(3.0_rp*press(ipoin,3)-press(ipoin,4)) + 0.5_rp*a(istep+1)*(press(ipoin,3)-press(ipoin,4)) !Capuano et al
             end do
          else
             do ipoin = 1,npoin
                do idime = 1,ndime
                   idofn = (ipoin-1)*ndime + idime
                   rm(idofn,istep) = bu(idofn)
                end do
                aux(ipoin) = pp(ipoin)  ! pressure solved from Laplacian
             end do
          endif
       else
          do ipoin = 1,npoin
             do idime = 1,ndime
                idofn = (ipoin-1)*ndime + idime
                rm(idofn,istep) = bu(idofn)
             end do
             aux(ipoin) = 0.0_rp
          end do
       end if
       !
       ! vv = GRAD*press
       !
       if( kfl_grad_div_nsi /= 0 )  then
          call nsi_aupvec(1_ip,Grad,aux,vv)           
       else
          call nsi_aupvec(1_ip,amatr(poaup_nsi:),aux,vv) 
       end if
       
       call rhsmod(ndime,rm(:,istep))
       call rhsmod(ndime,vv)

       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn = (ipoin-1)*ndime + idime
             aux2 = 0.0_rp
             do jstep = 1,istep
                aux2 = aux2 + rm(idofn,jstep)*b(istep+1,jstep)
             end do
             rt_tmp(idofn) = (aux2 - a(istep+1)*vv(idofn))
          end do
       end do
    end if
    !
    ! Solve system
    !
    call solver_explicit(solve_mome,rt_tmp,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)    
    !
    ! Update solution
    !
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu0(idofn) + (rt_tmp(idofn)/dtinv_nsi)/Nua(ipoin)       ! u = u0 + r * dt/rho
       end do
    end do
    call nsi_multi_step_fs_dirichlet()                                        ! Prescribe Drichlet bc on  u'^{k+1}

    call livinf(165_ip,'M',0_ip)
    call cputim(time2)
    cputi_nsi(1) = cputi_nsi(1) + time2 - time1

    if( kfl_fscon_nsi == 1 ) then 
       call nsi_multi_step_fs_pressure_and_correction(a(istep+1),0_ip)
    end if

  end subroutine nsi_multi_step_fs_solution_sj

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Momentum
  !> @details Fimal solution of the momentum equations
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_solution_final()

    integer(ip) :: idofn,ipoin,idime
    real(rp)    :: aux
    real(rp)    :: time1,time2

    call cputim(time1)

    if( INOTEMPTY ) then
       !
       ! Momentum residual rm^k
       !
       do idofn = 1,ndime*npoin                                        ! rm = bu 
          rm(idofn,4) = bu(idofn)  
       end do
       call rhsmod(ndime,rm(:,4))

       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn         = (ipoin-1)*ndime + idime
             aux           = rm(idofn,1)*sigma(1) + rm(idofn,2)*sigma(2) + rm(idofn,3)*sigma(3) + rm(idofn,4)*sigma(4)
             rt(idofn)     = aux
             rt_tmp(idofn) = (aux)
          end do
       end do
    end if
    !
    ! Solve
    !
    call solver_explicit(solve_mome,rt_tmp,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu0(idofn) + (rt_tmp(idofn)/dtinv_nsi)/Nua(ipoin) ! ORIOL
       end do
    end do
    call nsi_multi_step_fs_dirichlet()                            ! Prescribe Dirichlet bc on u'^{k+1}

    call livinf(165_ip,'M',0_ip)
    call cputim(time2)
    cputi_nsi(1) = cputi_nsi(1) + time2 - time1

    call nsi_multi_step_fs_pressure_and_correction(1.0_rp,1_ip)

  end subroutine nsi_multi_step_fs_solution_final

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Impose Dirichlet
  !> @details Impose Dirichlet conditions
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_dirichlet(xarray)

    real(rp),   intent(inout), pointer, optional :: xarray(:)
    integer(ip)                                  :: idofn,ipoin,idime

    if( present(xarray) ) then
       idofn = 0
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn = idofn + 1
             if( kfl_fixno_nsi(idime,ipoin) > 0 .and. kfl_noslw_ker /= 0) &
                  vafor_nsi(idime,ipoin) = xarray(idofn) ! actualy in the multistep case we could think of averaging all the substeps
             xarray(idofn) = 0.0_rp
          end do
       end do
    else
       idofn = 0
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn = idofn + 1
             if( kfl_fixno_nsi(idime,ipoin) > 0 ) &
                  uu(idofn) = bvess_nsi(idime,ipoin,1)
          end do
       end do
    end if

  end subroutine nsi_multi_step_fs_dirichlet

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Pressure and correction
  !> @details Solve pressure equation and correct velocity
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_pressure_and_correction(fact_substep,eval_rc)

    real(rp),    intent(in) :: fact_substep
    integer(ip), intent(in) :: eval_rc
    integer(ip)             :: idofn,ipoin,idime,jj
    real(rp)                :: aux2
    real(rp)                :: time1,time2

    !-----------------------------------------------------------------
    !
    ! Solve pressure
    !
    !-----------------------------------------------------------------

    call cputim(time1)
    if( INOTEMPTY ) then
       !
       ! Scale solution
       !
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn = (ipoin-1)*ndime + idime
             uu(idofn) = uu(idofn)*Nua(ipoin)
          end do
       end do
       !
       ! Continuity residual rc'^{k+1}
       !        
       if( kfl_grad_div_nsi /= 0 ) then
          call nsi_apuvec(1_ip,Div,uu,divup)
       else
          call nsi_apuvec(1_ip,amatr(poapu_nsi:),uu,divup)             ! rc' =  Apu u*
       end if
       !
       ! Residual
       !
       do ipoin = 1,npoin
          rc(ipoin) = (bp(ipoin) - divup(ipoin) )*(dtinv_nsi/fact_substep)
       end do
       call nsi_appvec(0_ip,Q,pp,rc,-1.0_rp)                           ! rc <= rc - Q p^k
       call rhsmod(1_ip,rc)
       !
       ! Add mass source from spray
       ! mass_sink is already exchanged, thus it's after rhsmod
       !
       if (associated(mass_sink)) then
          do ipoin = 1,npoin
             rc(ipoin) = rc(ipoin) +  mass_sink(ipoin)*(dtinv_nsi/fact_substep)
          end do
       endif
       !
       ! Low mach: drho/dt from old values.
       ! DRHODT is assumed to be already assembled
       !
       if( kfl_regim_nsi == 3 ) then
          if (ittim /= 1 .or. kfl_rstar /= 0 ) then
             call nsi_multi_step_fs_drho_dt(0_ip, fact_substep)
          else
             do ipoin = 1,npoin
                drhodt_nsi(ipoin) = 0.0_rp
             end do
          endif
          do ipoin = 1,npoin
             rc(ipoin) = rc(ipoin) - drhodt_nsi(ipoin)*(dtinv_nsi/fact_substep)
          end do
       end if

       do ipoin = 1,npoin
          Dp(ipoin) = 0.0_rp
       end do

    end if
    !
    ! Solve: Q z = rc'^k
    !
    call nsi_solini(2_ip)                                                         ! Initialize solver

    call solver_solve(solve_sol,Q,rc,Dp,EXCHANGE_RHS=.false.,ONLY_SOLVE=.true.)   ! Solve system Q Dp = rc^k

    if( INOTMASTER ) then
       !
       ! Dp      = Dp' + p^k
       ! p^{k+1} = Dp  +  p^k
       !
       do ipoin = 1,npoin
          Dp(ipoin) = Dp(ipoin) + pp(ipoin)  
          pp(ipoin) = Dp(ipoin) !- fsrot_nsi*Nu(ipoin)*divup(ipoin)
       end do
    end if

    if(eval_rc == 1) call norm2x(1_ip,rc,resin_nsi(2))                                  ! || rc ||

    call cputim(time2)
    cputi_nsi(2) = cputi_nsi(2) + time2 - time1
    call livinf(165_ip,'S',0_ip)

    !-----------------------------------------------------------------
    !
    ! Correct velocity
    !
    !-----------------------------------------------------------------
    
    call cputim(time1)
    if( INOTMASTER ) then

       if( kfl_grad_div_nsi /= 0 )  then
          call nsi_aupvec(1_ip,Grad,Dp,vv)                              ! rm =  G Dpp^{k+1}
       else
          call nsi_aupvec(1_ip,amatr(poaup_nsi:),Dp,vv)                 ! rm =  Aup Dpp^{k+1}
       end if
       call rhsmod(ndime,vv)
       !
       ! Obtain vafor without the contribution from the pressure
       !
#ifdef VAFOR_NO_PR
       rt_tmp = rt
       call nsi_multi_step_fs_dirichlet(rt_tmp)    ! This obtains vafor - done here the master will not enter but it does not matter
#endif
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn = (ipoin-1)*ndime + idime
             rt_tmp(idofn) = fact_substep*vv(idofn)
             rt(idofn) = rt(idofn) - vv(idofn) 
          end do
       end do
    end if
    !
    ! Solve
    !
    call solver_explicit(solve_mome,rt_tmp,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = (uu(idofn) - rt_tmp(idofn)/dtinv_nsi ) /Nua(ipoin)
       end do
    end do
    call nsi_multi_step_fs_dirichlet()                            ! Prescribe bc 

    if( kfl_matdi_nsi == NSI_DIRICHLET_ALGORITHM ) then
       do idofn = 1,npoin*ndime
          rt_tmp(idofn) = rt(idofn)
       end do
#ifndef VAFOR_NO_PR
       call nsi_multi_step_fs_dirichlet(rt_tmp)
#endif
       call norm2x(ndime,rt_tmp,resin_nsi(1))                        ! || rm ||
    else
       call norm2x(ndime,rt,resin_nsi(1))                            ! || rm ||       
    end if
    call cputim(time2)
    cputi_nsi(4) = cputi_nsi(4) + time2 - time1

  end subroutine nsi_multi_step_fs_pressure_and_correction

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Density
  !> @details Time derivative of density
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_drho_dt(i_ord, fact_substep)
    integer(ip), intent(in)     :: i_ord
    real(rp),    intent(in)     :: fact_substep
    real(rp)                    :: drho_dt, drho_dt0, drho_dt00, drho_dt000
    real(rp)                    :: aa, bb, cc, del_t_1, del_t_2, del_t_step
    integer(ip)                 :: ipoin

    select case (i_ord)
       
    case (-1_ip)
       !
       ! DO NOT CALCULATE
       return
       
    case (0_ip)
       !
       ! Time derivative is computed from two old values
       !
       do ipoin = 1,npoin
          drhodt_nsi(ipoin) = ( Nul(ipoin) - Nu0l(ipoin) ) * dtinv_nsi ! time derivative at n-0.5
       end do
       
    case (1_ip)
       ! Time derivative linearly extrapolated to current timestep: 
       ! drho/dt(s) = drho/dt(n-0.5) + (t(s)-t(n-0.5)) * d^2rho/dt^2  
       ! on paper this is equivalent to Nicoud's formulation in the final timestep 
       ! drho_dt = (((dt_inv0+dt_inv00)**2-dt_inv0**2)*Nul(ipoin) - ((dt_inv0+dt_inv00)**2)*Nu0l(ipoin) + (dt_inv0**2)*Nu00l(ipoin))/(dt_inv0*dt_inv00*(dt_inv0+dt_inv00))
       ! (F. Nicoud, Conservative high-order finite-difference schemes for low-Mach number flows,
       ! Journal of Computational Physics 158 (2000) 71-97)
       !do ipoin = 1,npoin
       !   drho_dt0   = (Nul(ipoin) -Nu0l(ipoin) ) * dtinv_nsi    ! time derivative at n-0.5
       !   drho_dt00  = (Nu0l(ipoin)-Nu00l(ipoin)) / dt_inv00     ! time derivative at n-1.5
       !   drhodt_nsi(ipoin) = drho_dt0 + (fact_substep/dtinv_nsi - 0.5_rp/dtinv_nsi)  *  (drho_dt0 - drho_dt00) / ( 0.5_rp * (dt_inv00 + 1.0_rp/dtinv_nsi) )
       !enddo

    case (2_ip)
       ! Parabola through the 3 drho/dt with t_loc=0 at n-0.5
       ! drho/dt (t) = a t^2 + b t + c
       ! 
       !           dt_inv000        dt_inv00      1/dtinv_nsi
       !       |================|==============|===============|
       !      n-3              n-2            n-1              n
       !             n-2.5           n-1.5           n-0.5    
       !          drho_dt000       drho_dt00       drho_dt0
       !         t = -del_t_2     t = -del_t_1       t = 0
       ! 
       !do ipoin = 1,npoin
       !   drho_dt0   = (Nul(ipoin) -Nu0l(ipoin) ) * dtinv_nsi    ! time derivative at n-0.5
       !   drho_dt00  = (Nu0l(ipoin)-Nu00l(ipoin)) / dt_inv00     ! time derivative at n-1.5
       !   drho_dt000 = (Nu00l(ipoin)-Nu000l(ipoin)) / dt_inv000  ! time derivative at n-2.5

       !   aa = 0.0_rp
       !   bb = 0.0_rp
       !   cc = 0.0_rp
       !   del_t_1 = 0.5_rp*(dt_inv00 + 1.0_rp/dtinv_nsi )
       !   del_t_2 = 0.5_rp*(dt_inv000 + 1.0_rp/dtinv_nsi ) + dt_inv00
       !   cc = drho_dt0
       !   aa = ( (drho_dt00-cc) - del_t_1/del_t_2 * (drho_dt000-cc) ) / ( del_t_1**2 - del_t_1*del_t_2 )
       !   bb = ( (drho_dt000-cc) - del_t_2**2 * aa  ) / (-1.0_rp * del_t_2)
       !   del_t_step = fact_substep/dtinv_nsi - 0.5_rp/dtinv_nsi
       !   drhodt_nsi(ipoin) = aa * del_t_step**2 + bb * del_t_step + cc
       !enddo
    end select

  end subroutine nsi_multi_step_fs_drho_dt

end module mod_nsi_multi_step_fs
