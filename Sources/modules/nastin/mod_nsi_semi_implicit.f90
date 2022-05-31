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
module mod_nsi_semi_implicit

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
  use mod_matrix,                      only : matrix_copy_matrix_to_block
  use mod_matrices,                    only : matrices_gradient_divergence
  use mod_matrices,                    only : matrices_laplacian
  use mod_nsi_dirichlet_global_system, only : nsi_dirichlet_matrix_rhs
  use mod_solver,                      only : solver_lumped_mass_system
  use mod_messages,                    only : livinf
  use mod_messages,                    only : messages_live
  use mod_solver,                      only : solver_solve
  use mod_solver,                      only : solver_exchange
  use mod_solver,                      only : solver_explicit
  use mod_solver,                      only : solver_copy_matrix, solver_add_matrix
  use mod_solver,                      only : solver_add_diagonal, solver_SpMV
  use mod_local_basis,                 only : local_basis_global_to_local
  use mod_local_basis,                 only : local_basis_local_to_global
  use mod_mass_matrix,                 only : mass_matrix_consistent
  implicit none
  private

  real(rp)               :: gamma1
  real(rp)               :: gamma

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
  real(rp),   pointer    :: Visc(:)

  real(rp),   pointer    :: Tim(:)      ! Nodal projection of dt/rho
  real(rp),   pointer    :: Mass(:)     ! Mass matrix
  real(rp),   pointer    :: Grad(:,:)   ! Gradient matrix   G_ij=-\int (div Ni) Nj
  real(rp),   pointer    :: Div(:,:)    ! Divergence matrix D_ij= \int Ni (div Nj) => D=-G^t
  real(rp),   pointer    :: Lapl(:)     ! Laplacian matrix
  real(rp),   pointer    :: ViscM(:)    ! Modified viscous term matrix

  real(rp)               :: c3 

  real(rp),   pointer    :: Nua(:)      ! Nodal projection of mu/rho
  real(rp),   pointer    :: Nul(:)      ! Nodal projection of M rho not updated
  real(rp),   pointer    :: Nu0l(:)     ! Nodal projection of M rho not updated

  real(rp)               :: rho(1)
  real(rp)               :: dt_inv0
  real(rp)               :: dt_inv00
  real(rp)               :: dt_inv000
  integer(ip)            :: time_iter
  integer(ip)            :: dummi

  real(rp), pointer :: rhs(:)
  real(rp), pointer :: bun(:)
  real(rp), pointer :: bu1(:)
  real(rp), pointer :: bu2(:)
  real(rp), pointer :: bu3(:)
  real(rp), pointer :: Mun(:)
  real(rp), pointer :: Auuu0(:)
  real(rp), pointer :: Auuu1(:)
  real(rp), pointer :: Auuu2(:)
  real(rp), pointer :: Auuu3(:)
  real(rp), pointer :: Auuu4(:)
  real(rp), pointer :: cm(:)


 
  type(soltyp), pointer  :: solve_mome(:) ! Momentum equations (explicit system) 
  type(soltyp), pointer  :: solve_mass(:) ! Consistent mass sytem
  type(soltyp), pointer  :: solve_visc(:) ! Viscous term solver

  public :: nsi_semi_implicit_solution
  public :: nsi_semi_implicit_initialization
  public :: nsi_semi_implicit_memory
  public :: nsi_semi_implicit_matrices

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

  subroutine nsi_semi_implicit_initialization()

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
    nullify(ViscM)

    nullify(Nua)    
    nullify(Nul)    
    nullify(Nu0l)   

    nullify(rhs,bun,bu1,bu2,bu3,mun,Auuu1,Auuu2,Auuu3,Auuu4,Auuu0,cm)

    time_iter = 1

  end subroutine nsi_semi_implicit_initialization
  
  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Compute some matrices
  !> @details Compute some matrices
  !
  !-----------------------------------------------------------------------
  
  subroutine nsi_semi_implicit_matrices()
    
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
    Visc   => visco_nsi 
    
    if( INOTEMPTY ) then
       Nul    => mass_rho_nsi(:,1)
       Nu0l   => mass_rho_nsi(:,2)
    else
       nullify(Nul)
       nullify(Nu0l)
    end if
    Mass   => vmass
  
  end subroutine nsi_semi_implicit_matrices

  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Allocate memory
  !> @details Allocate memory and define module constants
  !
  !-----------------------------------------------------------------------

  subroutine nsi_semi_implicit_memory()

    integer(ip) :: ipoin,nz
    integer(ip)       :: nn
    
    nn = 0
    gamma1 = 1.0_rp - gamma_nsi
    nz     = solve(1) % nzmat


    call memory_alloca(mem_modul(1:2,modul),'VV'    ,'mod_nsi_semi_implicit',vv,    max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'XX'    ,'mod_nsi_semi_implicit',xx,    max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'UU0'   ,'mod_nsi_semi_implicit',uu0,   max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'RM'    ,'mod_nsi_semi_implicit',rm,    max(1_ip,ndime*npoin),5_ip) 
    call memory_alloca(mem_modul(1:2,modul),'RT'    ,'mod_nsi_semi_implicit',rt,    max(1_ip,ndime*npoin)) 
    call memory_alloca(mem_modul(1:2,modul),'RTTMP' ,'mod_nsi_semi_implicit',rt_tmp,max(1_ip,ndime*npoin)) 
    call memory_alloca(mem_modul(1:2,modul),'RC'    ,'mod_nsi_semi_implicit',rc,    max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'AUX'   ,'mod_nsi_semi_implicit',aux,   max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'DIVUP' ,'mod_nsi_semi_implicit',divup, max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'DP'    ,'mod_nsi_semi_implicit',Dp,    max(1_ip,npoin))
    call memory_alloca(mem_modul(1:2,modul),'AUPP'  ,'mod_nsi_semi_implicit',Aupp,  max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'NUA'   ,'mod_nsi_semi_implicit',Nua,   max(1_ip,npoin))


    call memory_alloca(mem_modul(1:2,modul),'RHS'  ,'mod_nsi_semi_implicit',rhs,  max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'BUN'  ,'mod_nsi_semi_implicit',bun,  max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'BU1'  ,'mod_nsi_semi_implicit',bu1,  max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'BU2'  ,'mod_nsi_semi_implicit',bu2,  max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'BU3'  ,'mod_nsi_semi_implicit',bu3,  max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'MUN'  ,'mod_nsi_semi_implicit',Mun,  max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'AUUU0','mod_nsi_semi_implicit',Auuu0,max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'AUUU1','mod_nsi_semi_implicit',Auuu1,max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'AUUU2','mod_nsi_semi_implicit',Auuu2,max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'AUUU3','mod_nsi_semi_implicit',Auuu3,max(1_ip,ndime*npoin))
    call memory_alloca(mem_modul(1:2,modul),'AUUU4','mod_nsi_semi_implicit',Auuu4,max(1_ip,ndime*npoin))


    solve_mome => solve(NSI_SOLVER_MOMENTUM:)
    solve_mass => solve(NSI_SOLVER_CONSISTENT_MASS:)
    solve_visc => solve(NSI_SOLVER_VISCOUS_TERM:)

    if( memory_size(cmass) > 0 ) then
       nn = memory_size(cmass)* solve_visc(1) % ndofn * solve_visc(1) % ndofn
       call memory_alloca(memit,'CM','mod_nsi_semi_implicit',cm,nn)
       call matrix_copy_matrix_to_block(npoin,solve_visc(1) % ndofn,r_dom,c_dom,cmass,cm)
    else
       call memory_alloca(memit,'CM','mod_nsi_semi_implicit',cm,1_ip)             
    end if

    
    do ipoin = 1,npoin
       Nua(ipoin) = 1.0_rp
    end do
    
  end subroutine nsi_semi_implicit_memory

  !-----------------------------------------------------------------------
  !> 
  !> @author  lehmkuhl and houzeaux
  !> @date    2020-03-21
  !> @brief   Multi-step FS
  !> @details Multistep fractional step method from:
  !
  !-----------------------------------------------------------------------
  
  subroutine nsi_semi_implicit_solution()

    integer(ip) :: idofn,ipoin,idime
    real(rp)    :: time1,time2

    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn = (ipoin-1)*ndime + idime
          uu0(idofn)   = uu(idofn)
       end do
    end do
    !  if( kfl_regim_nsi == 3 ) then
    !    do ipoin = 1,npoin
    !       Nu0l(ipoin) = Nul(ipoin)
    !     end do
    !  end if

    !  if( kfl_regim_nsi == 3 .and. kfl_coupl(ID_TEMPER,ID_CHEMIC) /= 0 ) then
    !     !
    !     ! Low Mach
    !     !
    !     do ipoin = 1,npoin
    !        if (kfl_lookg_nsi > 0) then
    !           call ker_proper('DENSI','IPOIN',ipoin,dummi,rho)
    !           Nua(ipoin)  = rho(1)
    !        else
    !           Nua(ipoin)  = (prthe(1)/gasco) * (wmean(ipoin,1)/tempe(ipoin,1))
    !        endif
    !     end do
    !     !Nua = 0.5*(3.0_rp*Nu-Nu0)
    !
    !  else if( kfl_regim_nsi == 3 .and. (time_iter /= 1) .and. kfl_coupl(ID_TEMPER,ID_CHEMIC) == 0 ) then
    !
    !     do ipoin = 1,npoin
    !        Nua(ipoin)  = (prthe(1)/gasco) * (1.0_rp/tempe(ipoin,1))
    !     end do
    ! 
    !  else if( kfl_surte_nsi == 2_ip ) then
    !
    !     Nua = 1.0_rp
    !
    ! else
    !    !
    !    ! Incompressible
    !    !
    !    do ipoin = 1,npoin
    !       Nua(ipoin) = mass_rho_nsi(ipoin,1)
    !    end do
    !    call solver_lumped_mass_system(1_ip,Nua,EXCHANGE=.false.)
    ! 
    ! end if
    call cputim(time1)
    call livinf(56_ip,' ',modul)    
    call messages_live(' (',  ADVANCE='no')

    if(kfl_tiacc_nsi == 2) then 
       call nsi_semi_implicit_viscous_matrix2()
    else if (kfl_tiacc_nsi == 3) then 
       call nsi_semi_implicit_viscous_matrix5()
    else 
       call runend ('IMEX is only avaliable with 2 and 3 order')
    end if
    
    call messages_live(')')

    dt_inv000 = dt_inv00 
    dt_inv00  = dt_inv0 
    dt_inv0   = 1.0_rp/dtinv_nsi
    time_iter = time_iter + 1

    call cputim(time2)
    cpu_ass_sol_nsi(2) = time2 - time1

  end subroutine nsi_semi_implicit_solution

  ! CN+RK3
  subroutine nsi_semi_implicit_viscous_matrix2()

    integer(ip)       :: ipoin,idime,itotv,idofn,nn
    real(rp)          :: xfact

    call solver_SpMV(solve_visc(1),cm,uu0,mun)

    !
    ! u1
    !
    call messages_live('M',ADVANCE='no')
    
    call nsi_solini(3_ip)                                     ! Initialize solver

    call local_basis_global_to_local(kfl_fixrs_nsi,uu)        ! Global to local
    call local_basis_global_to_local(kfl_fixrs_nsi,uu0)       ! Global to local

    solve(NSI_SOLVER_VISCOUS_TERM) % kfl_do_assembly = SOL_YES
    call nsi_matrix()
    
    if( ittim > 4 ) then 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bun(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = 0.5_rp*(3.0_rp*press(ipoin,3)-press(ipoin,4)) + 0.5_rp*(8.0_rp/15.0_rp)*(press(ipoin,3)-press(ipoin,4))
       end do
    else 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bun(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = pp(ipoin)
       end do
    end if
    call rhsmod(ndime,bun)

    call solver_SpMV(solve_visc(1),Visc,uu0,Auuu0)

    call solver_copy_matrix(solve_visc(1),Visc,ViscM)
    xfact = (4.0_rp/15.0_rp)* dtime
    call solver_add_matrix(solve_visc(1),cm,ViscM,MATRIX_FACTOR=xfact,ADD_FACTOR=1.0_rp)

    do itotv = 1,npoin*ndime
       rhs(itotv) = mun(itotv) - dtime*(4.0_rp/15.0_rp)*Auuu0(itotv) + dtime*(8.0_rp/15.0_rp)*bun(itotv)
    end do
    call solver_solve(solve_visc,ViscM,rhs,uu,EXCHANGE_RHS=.false.,ONLY_SOLVE=.true.)
    !call nsi_semi_implicit_pressure_and_correction(8.0_rp/15.0_rp,0_ip) 

    call nsi_aupvec(1_ip,Grad,aux,vv)  
    call rhsmod(ndime,vv)
    call solver_explicit(solve_mome,vv,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu(idofn) - (8.0_rp/15.0_rp)*vv(idofn)*dtime /Nua(ipoin)
       end do
    end do
    call nsi_semi_implicit_dirichlet() 

    call local_basis_local_to_global(kfl_fixrs_nsi,uu)
    call local_basis_local_to_global(kfl_fixrs_nsi,uu0)

    call nsi_updunk(ITASK_INNITE)


    !
    ! u2
    !
    call messages_live('M',ADVANCE='no')
    
    call nsi_solini(3_ip)                                     ! Initialize solver

    call local_basis_global_to_local(kfl_fixrs_nsi,uu)        ! Global to local
    call local_basis_global_to_local(kfl_fixrs_nsi,uu0)       ! Global to local

    solve(NSI_SOLVER_VISCOUS_TERM) % kfl_do_assembly = SOL_NO
    call nsi_matrix()

    if( ittim > 4 ) then 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bu1(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = 0.5_rp*(3.0_rp*press(ipoin,3)-press(ipoin,4)) + 0.5_rp*(2.0_rp/3.0_rp)*(press(ipoin,3)-press(ipoin,4))
       end do
    else 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bu1(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = pp(ipoin)
       end do
    end if
    call rhsmod(ndime,bu1)

    call solver_SpMV(solve_visc(1),Visc,uu,Auuu1)

    call solver_copy_matrix(solve_visc(1),Visc,ViscM)
    xfact = (1.0_rp/15.0_rp)* dtime
    call solver_add_matrix(solve_visc(1),cm,ViscM,MATRIX_FACTOR=xfact,ADD_FACTOR=1.0_rp)

    do itotv = 1,npoin*ndime
       rhs(itotv) = mun(itotv) + (1.0_rp/4.0_rp)*dtime * bun(itotv) + (5.0_rp/12.0_rp)*dtime * bu1(itotv)& 
                  - dtime*(4.0_rp/15.0_rp) * Auuu0(itotv) - dtime*(1.0_rp/3.0_rp) * Auuu1(itotv)
    end do
    call solver_solve(solve_visc,ViscM,rhs,uu,EXCHANGE_RHS=.false.,ONLY_SOLVE=.true.)

    !call nsi_semi_implicit_pressure_and_correction(2.0_rp/3.0_rp,0_ip) 
    
    call nsi_aupvec(1_ip,Grad,aux,vv)  
    call rhsmod(ndime,vv)
    call solver_explicit(solve_mome,vv,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu(idofn) - (2.0_rp/3.0_rp)*vv(idofn)*dtime /Nua(ipoin)
       end do
    end do
    call nsi_semi_implicit_dirichlet() 

    call local_basis_local_to_global(kfl_fixrs_nsi,uu)
    call local_basis_local_to_global(kfl_fixrs_nsi,uu0)

    call nsi_updunk(ITASK_INNITE)

    !
    ! u3
    !
    call messages_live('M',ADVANCE='no')
    call nsi_solini(3_ip)                                     ! Initialize solver

    call local_basis_global_to_local(kfl_fixrs_nsi,uu)        ! Global to local
    call local_basis_global_to_local(kfl_fixrs_nsi,uu0)       ! Global to local

    solve(NSI_SOLVER_VISCOUS_TERM) % kfl_do_assembly = SOL_NO
    call nsi_matrix()

    if( ittim > 4 ) then 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bu2(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = 0.5_rp*(3.0_rp*press(ipoin,3)-press(ipoin,4)) + 0.5_rp*(press(ipoin,3)-press(ipoin,4))
       end do
    else 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bu2(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = pp(ipoin)
       end do
    end if
    call rhsmod(ndime,bu2)

    call solver_SpMV(solve_visc(1),Visc,uu,Auuu2)

    call solver_copy_matrix(solve_visc(1),Visc,ViscM)
    xfact = (1.0_rp/6.0_rp)* dtime
    call solver_add_matrix(solve_visc(1),cm,ViscM,MATRIX_FACTOR=xfact,ADD_FACTOR=1.0_rp)

    do itotv = 1,npoin*ndime
       rhs(itotv) = mun(itotv) + (1.0_rp/4.0_rp)*dtime * bun(itotv) + (3.0_rp/4.0_rp)*dtime * bu2(itotv)& 
                  - dtime*(4.0_rp/15.0_rp) * Auuu0(itotv) - dtime*(1.0_rp/3.0_rp) * Auuu1(itotv) - dtime*(7.0_rp/30.0_rp) * Auuu2(itotv)
    end do
    call solver_solve(solve_visc,ViscM,rhs,uu,EXCHANGE_RHS=.false.,ONLY_SOLVE=.true.)

    call nsi_aupvec(1_ip,Grad,aux,vv)  
    call rhsmod(ndime,vv)
    call solver_explicit(solve_mome,vv,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu(idofn) - vv(idofn)*dtime /Nua(ipoin)
       end do
    end do
    call nsi_semi_implicit_dirichlet() 

    call local_basis_local_to_global(kfl_fixrs_nsi,uu)
    call local_basis_local_to_global(kfl_fixrs_nsi,uu0)

    call nsi_updunk(ITASK_INNITE)

    !
    ! u+1
    !
    call messages_live('M',ADVANCE='no')
    call nsi_solini(3_ip)                                     ! Initialize solver

    call local_basis_global_to_local(kfl_fixrs_nsi,uu)        ! Global to local
    call local_basis_global_to_local(kfl_fixrs_nsi,uu0)       ! Global to local

    call solver_SpMV(solve_visc(1),Visc,uu,Auuu3)

    do itotv = 1,npoin*ndime
       rhs(itotv) = (1.0_rp/4.0_rp)*dtime * bun(itotv) + (3.0_rp/4.0_rp)*dtime * bu2(itotv)& 
                  - dtime*(4.0_rp/15.0_rp) * Auuu0(itotv) - dtime*(1.0_rp/3.0_rp) * Auuu1(itotv) - dtime*(7.0_rp/30.0_rp) * Auuu2(itotv) &
                  - dtime*(1.0_rp/6.0_rp) * Auuu3(itotv)
    end do
    call solver_explicit(solve_mome,rhs,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu0(idofn) + rhs(idofn) /Nua(ipoin)
       end do
    end do
    call nsi_semi_implicit_dirichlet() 

    call nsi_semi_implicit_pressure_and_correction(1.0_rp,1_ip) 

    call local_basis_local_to_global(kfl_fixrs_nsi,uu)
    call local_basis_local_to_global(kfl_fixrs_nsi,uu0)


    !call memory_deallo(mem_modul(1:2,modul),'RHS'  ,'mod_nsi_semi_implicit',rhs)
    !call memory_deallo(mem_modul(1:2,modul),'BUN'  ,'mod_nsi_semi_implicit',bun)
    !call memory_deallo(mem_modul(1:2,modul),'BU1'  ,'mod_nsi_semi_implicit',bu1)
    !call memory_deallo(mem_modul(1:2,modul),'BU2'  ,'mod_nsi_semi_implicit',bu2)
    !call memory_deallo(mem_modul(1:2,modul),'MUN'  ,'mod_nsi_semi_implicit',Mun)
    !call memory_deallo(mem_modul(1:2,modul),'AUUU0','mod_nsi_semi_implicit',Auuu0)
    !call memory_deallo(mem_modul(1:2,modul),'AUUU1','mod_nsi_semi_implicit',Auuu1)
    !call memory_deallo(mem_modul(1:2,modul),'AUUU2','mod_nsi_semi_implicit',Auuu2)
    !call memory_deallo(mem_modul(1:2,modul),'AUUU3','mod_nsi_semi_implicit',Auuu3)
    !call memory_deallo(mem_modul(1:2,modul),'U1'   ,'mod_nsi_semi_implicit',u1)
    !call memory_deallo(mem_modul(1:2,modul),'U2'   ,'mod_nsi_semi_implicit',u2)
    !call memory_deallo(mem_modul(1:2,modul),'U3'   ,'mod_nsi_semi_implicit',u3)
    !call memory_deallo(mem_modul(1:2,modul),'CM'   ,'mod_nsi_semi_implicit',cm)

  end subroutine nsi_semi_implicit_viscous_matrix2

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Impose Dirichlet
  !> @details Impose Dirichlet conditions
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_semi_implicit_dirichlet(xarray)

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

  end subroutine nsi_semi_implicit_dirichlet

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Pressure and correction
  !> @details Solve pressure equation and correct velocity
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_semi_implicit_pressure_and_correction(fact_substep,eval_rc)

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

    call messages_live('S',ADVANCE='no')
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
             call nsi_semi_implicit_drho_dt(0_ip, fact_substep)
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
    call messages_live('C',ADVANCE='no')

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
       call nsi_semi_implicit_dirichlet(rt_tmp)    ! This obtains vafor - done here the master will not enter but it does not matter
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
    call nsi_semi_implicit_dirichlet()                            ! Prescribe bc 

    if( kfl_matdi_nsi == NSI_DIRICHLET_ALGORITHM ) then
       do idofn = 1,npoin*ndime
          rt_tmp(idofn) = rt(idofn)
       end do
#ifndef VAFOR_NO_PR
       call nsi_semi_implicit_dirichlet(rt_tmp)
#endif
       call norm2x(ndime,rt_tmp,resin_nsi(1))                        ! || rm ||
    else
       call norm2x(ndime,rt,resin_nsi(1))                            ! || rm ||       
    end if
    call cputim(time2)
    cputi_nsi(4) = cputi_nsi(4) + time2 - time1

  end subroutine nsi_semi_implicit_pressure_and_correction

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-03-20
  !> @brief   Density
  !> @details Time derivative of density
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_semi_implicit_drho_dt(i_ord, fact_substep)
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

  end subroutine nsi_semi_implicit_drho_dt

  subroutine nsi_semi_implicit_viscous_matrix5()

    integer(ip)       :: ipoin,idime,itotv,idofn,nn
    real(rp)          :: xfact
    real(rp)          :: b1_i,b2_i,b3_i,b4_i,b5_i,a21_i,a22_i,a32_i,a33_i,a43_i,a44_i
    real(rp)          :: b1_e,b2_e,b3_e,b4_e,a21_e,a32_e,a43_e
    real(rp)          :: c2,c3,c4

    b1_i = 0.0_rp
    b2_i = 4078465402807.0_rp/5118992086463.0_rp
    b3_i = -1068609889687.0_rp/1488061735778.0_rp
    b4_i = 94533336407.0_rp/126055292720.0_rp
    b5_i = 112416685574.0_rp/655665149019.0_rp
    a21_i = 0.0_rp
    a22_i = 9.0_rp/25.0_rp 
    a32_i = 379490756215.0_rp/588608184103.0_rp
    a33_i = 81921593785.0_rp/419520366036.0_rp
    a43_i = -39458195936308.0_rp/94684797045601.0_rp
    a44_i = 12.0_rp/25.0_rp

    b1_e = 67447694372.0_rp/739814670703.0_rp
    b2_e = 640712099409.0_rp/1099358471078.0_rp
    b3_e = -299809194319.0_rp/611386646053.0_rp
    b4_e = 878905218902.0_rp/1076559421011.0_rp
    a21_e = 9.0_rp/25.0_rp
    a32_e = 869434674241.0_rp/1161054947863.0_rp
    a43_e = 359201878931.0_rp/1930920984086.0_rp

    c2 = 9.0_rp/25.0_rp
    c3 = 21.0_rp/25.0_rp
    c4 = 43.0_rp/50.0_rp

    call solver_SpMV(solve_visc(1),cm,uu0,mun)

    !
    ! u1
    !
    call messages_live('M',ADVANCE='no')

    call nsi_solini(3_ip)                                     ! Initialize solver

    call local_basis_global_to_local(kfl_fixrs_nsi,uu)        ! Global to local
    call local_basis_global_to_local(kfl_fixrs_nsi,uu0)       ! Global to local

    solve(NSI_SOLVER_VISCOUS_TERM) % kfl_do_assembly = SOL_YES
    call nsi_matrix()
    
    if( ittim > 4 ) then 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bun(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = 0.5_rp*(3.0_rp*press(ipoin,3)-press(ipoin,4)) + 0.5_rp*c2*(press(ipoin,3)-press(ipoin,4))
       end do
    else 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bun(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = pp(ipoin)
       end do
    end if
    call rhsmod(ndime,bun)

    call solver_SpMV(solve_visc(1),Visc,uu0,Auuu0)

    call solver_copy_matrix(solve_visc(1),Visc,ViscM)
    xfact = (a22_i)* dtime
    call solver_add_matrix(solve_visc(1),cm,ViscM,MATRIX_FACTOR=xfact,ADD_FACTOR=1.0_rp)

    do itotv = 1,npoin*ndime
       rhs(itotv) = mun(itotv) - dtime*(a21_i)*Auuu0(itotv) + dtime*(a21_e)*bun(itotv) 
    end do
    call solver_solve(solve_visc,ViscM,rhs,uu,EXCHANGE_RHS=.false.,ONLY_SOLVE=.true.)

    !call nsi_semi_implicit_pressure_and_correction(c2,0_ip) 

    call nsi_aupvec(1_ip,Grad,aux,vv)  
    call rhsmod(ndime,vv)
    call solver_explicit(solve_mome,vv,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu(idofn) - c2*vv(idofn)*dtime /Nua(ipoin)
       end do
    end do
    call nsi_semi_implicit_dirichlet() 
    call local_basis_local_to_global(kfl_fixrs_nsi,uu)
    call local_basis_local_to_global(kfl_fixrs_nsi,uu0)

    call nsi_updunk(ITASK_INNITE)

    !
    ! u2
    !
    call messages_live('M',ADVANCE='no')

    call nsi_solini(3_ip)                                     ! Initialize solver

    call local_basis_global_to_local(kfl_fixrs_nsi,uu)        ! Global to local
    call local_basis_global_to_local(kfl_fixrs_nsi,uu0)       ! Global to local

    solve(NSI_SOLVER_VISCOUS_TERM) % kfl_do_assembly = SOL_NO
    call nsi_matrix()

    if( ittim > 4 ) then 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bu1(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = 0.5_rp*(3.0_rp*press(ipoin,3)-press(ipoin,4)) + 0.5_rp*c3*(press(ipoin,3)-press(ipoin,4)) 
       end do
    else 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bu1(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = pp(ipoin)
       end do
    end if
    call rhsmod(ndime,bu1)

    call solver_SpMV(solve_visc(1),Visc,uu,Auuu1)

    call solver_copy_matrix(solve_visc(1),Visc,ViscM)
    xfact = (a33_i)* dtime
    call solver_add_matrix(solve_visc(1),cm,ViscM,MATRIX_FACTOR=xfact,ADD_FACTOR=1.0_rp)

    do itotv = 1,npoin*ndime
       rhs(itotv) = mun(itotv) + (b1_e)*dtime * bun(itotv) + (a32_e)*dtime * bu1(itotv)& 
                  - dtime*(b1_i) * Auuu0(itotv) - dtime*(a32_i) * Auuu1(itotv) 
    end do
    call solver_solve(solve_visc,ViscM,rhs,uu,EXCHANGE_RHS=.false.,ONLY_SOLVE=.true.)

    !call nsi_semi_implicit_pressure_and_correction(c3,0_ip) 

    call nsi_aupvec(1_ip,Grad,aux,vv)  
    call rhsmod(ndime,vv)
    call solver_explicit(solve_mome,vv,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu(idofn) - c3*vv(idofn)*dtime /Nua(ipoin)
       end do
    end do
    call nsi_semi_implicit_dirichlet() 

    call local_basis_local_to_global(kfl_fixrs_nsi,uu)
    call local_basis_local_to_global(kfl_fixrs_nsi,uu0)

    call nsi_updunk(ITASK_INNITE)

    !
    ! u3
    !
    call messages_live('M',ADVANCE='no')
    call nsi_solini(3_ip)                                     ! Initialize solver

    call local_basis_global_to_local(kfl_fixrs_nsi,uu)        ! Global to local
    call local_basis_global_to_local(kfl_fixrs_nsi,uu0)       ! Global to local

    solve(NSI_SOLVER_VISCOUS_TERM) % kfl_do_assembly = SOL_NO
    call nsi_matrix()

    if( ittim > 4 ) then 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bu2(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = 0.5_rp*(3.0_rp*press(ipoin,3)-press(ipoin,4)) + 0.5_rp*c4*(press(ipoin,3)-press(ipoin,4)) 
       end do
    else 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bu2(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = pp(ipoin)
       end do
    end if
    call rhsmod(ndime,bu2)

    call solver_SpMV(solve_visc(1),Visc,uu,Auuu2)

    call solver_copy_matrix(solve_visc(1),Visc,ViscM)
    xfact = (a44_i)* dtime
    call solver_add_matrix(solve_visc(1),cm,ViscM,MATRIX_FACTOR=xfact,ADD_FACTOR=1.0_rp)

    do itotv = 1,npoin*ndime
       rhs(itotv) = mun(itotv) + (b1_e)*dtime * bun(itotv) + (b2_e)*dtime * bu1(itotv)+  (a43_e)*dtime * bu2(itotv)& 
                  - dtime*(b1_i) * Auuu0(itotv) - dtime*(b2_i) * Auuu1(itotv) - dtime*(a43_i) * Auuu2(itotv) 
    end do
    call solver_solve(solve_visc,ViscM,rhs,uu,EXCHANGE_RHS=.false.,ONLY_SOLVE=.true.)

    !call nsi_semi_implicit_pressure_and_correction(c4,0_ip) 

    call nsi_aupvec(1_ip,Grad,aux,vv)  
    call rhsmod(ndime,vv)
    call solver_explicit(solve_mome,vv,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu(idofn) - c4*vv(idofn)*dtime /Nua(ipoin)
       end do
    end do
    call nsi_semi_implicit_dirichlet() 

    call local_basis_local_to_global(kfl_fixrs_nsi,uu)
    call local_basis_local_to_global(kfl_fixrs_nsi,uu0)

    call nsi_updunk(ITASK_INNITE)

    !
    ! u4
    !
    call messages_live('M',ADVANCE='no')
    call nsi_solini(3_ip)                                     ! Initialize solver

    call local_basis_global_to_local(kfl_fixrs_nsi,uu)        ! Global to local
    call local_basis_global_to_local(kfl_fixrs_nsi,uu0)       ! Global to local

    solve(NSI_SOLVER_VISCOUS_TERM) % kfl_do_assembly = SOL_NO
    call nsi_matrix()

    if( ittim > 4 ) then 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bu3(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = 0.5_rp*(3.0_rp*press(ipoin,3)-press(ipoin,4)) + 0.5_rp*(press(ipoin,3)-press(ipoin,4))
       end do
    else 
       do ipoin = 1,npoin
          do idime = 1,ndime
             itotv = (ipoin-1)*ndime + idime
             bu3(itotv) = rhsid(itotv)
          end do
          aux(ipoin) = pp(ipoin)
       end do
    end if
    call rhsmod(ndime,bu3)

    call solver_SpMV(solve_visc(1),Visc,uu,Auuu3)

    call solver_copy_matrix(solve_visc(1),Visc,ViscM)
    xfact = (b5_i)* dtime
    call solver_add_matrix(solve_visc(1),cm,ViscM,MATRIX_FACTOR=xfact,ADD_FACTOR=1.0_rp)

    do itotv = 1,npoin*ndime
       rhs(itotv) = mun(itotv) + (b1_e)*dtime*bun(itotv) + (b2_e)*dtime*bu1(itotv)+  (b3_e)*dtime*bu2(itotv)+ (b4_e)*dtime*bu3(itotv) &
                  - dtime*(b1_i) * Auuu0(itotv) - dtime*(b2_i) * Auuu1(itotv) - dtime*(b3_i) * Auuu2(itotv)- dtime*(b4_i) * Auuu3(itotv)
    end do
    call solver_solve(solve_visc,ViscM,rhs,uu,EXCHANGE_RHS=.false.,ONLY_SOLVE=.true.)

    call nsi_aupvec(1_ip,Grad,aux,vv)  
    call rhsmod(ndime,vv)
    call solver_explicit(solve_mome,vv,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu(idofn) - vv(idofn)*dtime /Nua(ipoin)
       end do
    end do
    call nsi_semi_implicit_dirichlet() 

    call local_basis_local_to_global(kfl_fixrs_nsi,uu)
    call local_basis_local_to_global(kfl_fixrs_nsi,uu0)

    call nsi_updunk(ITASK_INNITE)
    !
    ! u+1
    !

    call messages_live('M',ADVANCE='no')
    call nsi_solini(3_ip)                                     ! Initialize solver

    call local_basis_global_to_local(kfl_fixrs_nsi,uu)        ! Global to local
    call local_basis_global_to_local(kfl_fixrs_nsi,uu0)       ! Global to local

    call solver_SpMV(solve_visc(1),Visc,uu,Auuu4)

    do itotv = 1,npoin*ndime
       rhs(itotv) = (b1_e)*dtime*bun(itotv) + (b2_e)*dtime*bu1(itotv)+  (b3_e)*dtime*bu2(itotv)+ (b4_e)*dtime*bu3(itotv) &
                  - dtime*(b1_i) * Auuu0(itotv) - dtime*(b2_i) * Auuu1(itotv) - dtime*(b3_i) * Auuu2(itotv) &
                  - dtime*(b4_i) * Auuu3(itotv) - dtime*(b5_i)*Auuu4(itotv)
    end do
    call solver_explicit(solve_mome,rhs,EXCHANGE=.false.,solve_consistent=solve_mass,x=uu)
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn     = (ipoin-1)*ndime + idime
          uu(idofn) = uu0(idofn) + rhs(idofn) /Nua(ipoin)
       end do
    end do
    call nsi_semi_implicit_dirichlet() 

    call nsi_semi_implicit_pressure_and_correction(1.0_rp,1_ip) 

    call local_basis_local_to_global(kfl_fixrs_nsi,uu)
    call local_basis_local_to_global(kfl_fixrs_nsi,uu0)

    !call nsi_updunk(ITASK_INNITE)

    !call memory_deallo(mem_modul(1:2,modul),'RHS'  ,'mod_nsi_semi_implicit',rhs)
    !call memory_deallo(mem_modul(1:2,modul),'BUN'  ,'mod_nsi_semi_implicit',bun)
    !call memory_deallo(mem_modul(1:2,modul),'BU1'  ,'mod_nsi_semi_implicit',bu1)
    !call memory_deallo(mem_modul(1:2,modul),'BU2'  ,'mod_nsi_semi_implicit',bu2)
    !call memory_deallo(mem_modul(1:2,modul),'BU3'  ,'mod_nsi_semi_implicit',bu3)
    !call memory_deallo(mem_modul(1:2,modul),'MUN'  ,'mod_nsi_semi_implicit',Mun)
    !call memory_deallo(mem_modul(1:2,modul),'AUUU0','mod_nsi_semi_implicit',Auuu0)
    !call memory_deallo(mem_modul(1:2,modul),'AUUU1','mod_nsi_semi_implicit',Auuu1)
    !call memory_deallo(mem_modul(1:2,modul),'AUUU2','mod_nsi_semi_implicit',Auuu2)
    !call memory_deallo(mem_modul(1:2,modul),'AUUU3','mod_nsi_semi_implicit',Auuu3)
    !call memory_deallo(mem_modul(1:2,modul),'AUUU4','mod_nsi_semi_implicit',Auuu4)
    !call memory_deallo(mem_modul(1:2,modul),'U1'   ,'mod_nsi_semi_implicit',u1)
    !call memory_deallo(mem_modul(1:2,modul),'U2'   ,'mod_nsi_semi_implicit',u2)
    !call memory_deallo(mem_modul(1:2,modul),'U3'   ,'mod_nsi_semi_implicit',u3)
    !call memory_deallo(mem_modul(1:2,modul),'U4'   ,'mod_nsi_semi_implicit',u4)
    !call memory_deallo(mem_modul(1:2,modul),'CM'   ,'mod_nsi_semi_implicit',cm)

  end subroutine nsi_semi_implicit_viscous_matrix5
end module mod_nsi_semi_implicit
!> @}

