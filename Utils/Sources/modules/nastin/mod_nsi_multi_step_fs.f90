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
  use mod_postpr
  use mod_memory,                      only : memory_alloca
  use mod_memory,                      only : memory_deallo
  use mod_nsi_schur_operations
  use mod_matrices,                    only : matrices_gradient_divergence
  use mod_matrices,                    only : matrices_laplacian
  use mod_nsi_dirichlet_global_system, only : nsi_dirichlet_laplacian
  use mod_solver,                      only : solver_periodicity
  use mod_solver,                      only : solver_periodicity_rhs
  use mod_solver,                      only : solver_lumped_mass_system
  use mod_messages,                    only : livinf

  implicit none
  private

  real(rp),            save :: gamma1

  real(rp),   pointer, save :: vv(:)       ! auxiliar vector
  real(rp),   pointer, save :: xx(:)       ! auxiliar vector
  real(rp),   pointer, save :: aux(:)      ! auxiliar vector
  real(rp),   pointer, save :: uu0(:)      ! auxiliar vector
  real(rp),   pointer, save :: rm(:,:)     ! residual of momentum at time 
  real(rp),   pointer, save :: rc(:)       ! residual of continuity 
  real(rp),   pointer, save :: divup(:)    ! div of up 
  real(rp),   pointer, save :: rt(:)       ! Total residual of momentum 
  real(rp),   pointer, save :: rt_tmp(:)   ! Total residual of momentum 
  real(rp),   pointer, save :: bu(:)         
  real(rp),   pointer, save :: bp(:)
  real(rp),   pointer, save :: uu(:)
  real(rp),   pointer, save :: pp(:)
  real(rp),   pointer, save :: Q(:)
  real(rp),   pointer, save :: Dp(:)
  real(rp),   pointer, save :: Aupp(:)
  real(rp),   pointer, save :: cmasm(:)    ! Consistent mass matrix
  real(rp),   pointer, save :: masbc(:)    ! RHS for consistent mass matrix bc

  real(rp),   pointer, save :: Tim(:)      ! Nodal projection of dt/rho
  real(rp),   pointer, save :: Nu(:)       ! Nodal projection of mu/rho
  real(rp),   pointer, save :: Mass(:)     ! Mass matrix
  real(rp),   pointer, save :: Grad(:,:)   ! Gradient matrix   G_ij=-\int (div Ni) Nj
  real(rp),   pointer, save :: Div(:,:)    ! Divergence matrix D_ij= \int Ni (div Nj) => D=-G^t
  real(rp),   pointer       :: Lapl(:)     ! Laplacian matrix

  real(rp),            save :: c3 
  real(rp),            save :: a(5)
  real(rp),            save :: b(5,5)
  real(rp),            save :: sigma(5)
  integer(ip), parameter    :: conservative = 0
  integer(ip), parameter    :: n_moin_steps = 4   ! Number of steps as Moin would do it. That is, without the enhanced extrapolation proposed by Capuano for the case kfl_fscon_nsi==0 

  public :: nsi_multi_step_fs_solution

contains 

  !-----------------------------------------------------------------------
  !
  !> @date    03/10/2016
  !> @author  Guillaume Houzeaux
  !> @brief   Allocate
  !> @details Allocate
  !
  !-----------------------------------------------------------------------

  subroutine nsi_multi_step_fs_allocate()

    integer(ip), save :: ipass = 0
    integer(ip)       :: idofn,ipoin,idime
    real(rp)          :: auxra, epsle = 0.05_rp


    if( ipass == 0 ) then

       ipass = 1
       nullify(vv)
       nullify(aux)
       nullify(uu0)
       nullify(rm)
       nullify(rt)
       nullify(rt_tmp)
       nullify(rc)
       nullify(divup)
       nullify(Dp)
       nullify(Grad)
       nullify(Div)
       nullify(cmasm)
       nullify(masbc)
       nullify(xx)

       bu     => rhsid
       bp     => rhsid(ndbgs_nsi+1:)
       uu     => unkno
       pp     => unkno(ndbgs_nsi+1:)   
       Tim    => dt_rho_nsi
       Nu     => mu_rho_nsi
       Mass   => vmass
       cmasm  => amatr
       gamma1 =  1.0_rp - gamma_nsi

       if( INOTMASTER ) then

          call memory_alloca(mem_modul(1:2,modul),'VV',   'mod_nsi_solsch',vv,    ndime*npoin)
          call memory_alloca(mem_modul(1:2,modul),'XX',   'mod_nsi_solsch',xx,    ndime*npoin)
          call memory_alloca(mem_modul(1:2,modul),'UU0',  'mod_nsi_solsch',uu0,   ndime*npoin)
          call memory_alloca(mem_modul(1:2,modul),'RM',   'mod_nsi_solsch',rm,    ndime*npoin,5_ip) 
          call memory_alloca(mem_modul(1:2,modul),'RT',   'mod_nsi_solsch',rt,    ndime*npoin) 
          call memory_alloca(mem_modul(1:2,modul),'RTTMP','mod_nsi_solsch',rt_tmp,ndime*npoin) 
          call memory_alloca(mem_modul(1:2,modul),'RC',   'mod_nsi_solsch',rc,    npoin)
          call memory_alloca(mem_modul(1:2,modul),'AUX',  'mod_nsi_solsch',aux,   npoin)
          call memory_alloca(mem_modul(1:2,modul),'DIVUP','mod_nsi_solsch',divup, npoin)
          call memory_alloca(mem_modul(1:2,modul),'DP',   'mod_nsi_solsch',Dp,    npoin)
          call memory_alloca(mem_modul(1:2,modul),'AUPP', 'mod_nsi_solsch',Aupp,  ndime*npoin)
          call memory_alloca(mem_modul(1:2,modul),'MASBC','mod_nsi_solsch',masbc, ndime*npoin) 

       else

          call memory_alloca(mem_modul(1:2,modul),'VV',   'mod_nsi_solsch',vv,    1_ip)
          call memory_alloca(mem_modul(1:2,modul),'XX',   'mod_nsi_solsch',xx,    1_ip)
          call memory_alloca(mem_modul(1:2,modul),'AUX',  'mod_nsi_solsch',aux,   1_ip)
          call memory_alloca(mem_modul(1:2,modul),'UUO',  'mod_nsi_solsch',uu0,   1_ip)
          call memory_alloca(mem_modul(1:2,modul),'RM',   'mod_nsi_solsch',rm,    1_ip,1_ip)
          call memory_alloca(mem_modul(1:2,modul),'RT',   'mod_nsi_solsch',rt,    1_ip)
          call memory_alloca(mem_modul(1:2,modul),'RTTMP','mod_nsi_solsch',rt_tmp,1_ip)
          call memory_alloca(mem_modul(1:2,modul),'RC',   'mod_nsi_solsch',rc,    1_ip)
          call memory_alloca(mem_modul(1:2,modul),'DIVUP','mod_nsi_solsch',divup, 1_ip)
          call memory_alloca(mem_modul(1:2,modul),'DP',   'mod_nsi_solsch',Dp,    1_ip)
          call memory_alloca(mem_modul(1:2,modul),'AUPP', 'mod_nsi_solsch',Aupp,  1_ip)
          call memory_alloca(mem_modul(1:2,modul),'MASBC','mod_nsi_solsch',masbc, 1_ip) 

       end if

       if( kfl_grad_div_nsi /= 0 )                  call matrices_gradient_divergence(Grad,Div)
       if( kfl_grad_div_nsi /= 0 )                  call matrices_laplacian(Lapl)
       if( kfl_grad_div_nsi /= 0 .and. INOTMASTER ) call nsi_dirichlet_laplacian(Lapl)  ! I added INOTMASTER as in nsi_matrix
       if( kfl_grad_div_nsi /= 0 ) then
          Q => Lapl
       else
          Q => lapla_nsi
       end if
       !
       ! Periodicity has already been imposed in Q in nsi_dirichlet_laplacian
       !
       if( kfl_grad_div_nsi /= 0 ) then
          solve_sol => solve(2:)
          !call solver_periodicity('MATRIX',solve_sol,1_ip, 1_ip,Lapl,message='KEEP DIAGONAL')
          call solver_periodicity('MATRIX',solve_sol,ndime,1_ip,Grad)
          call solver_periodicity('MATRIX',solve_sol,ndime,1_ip,Div)
       end if
       
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
    end if

  end subroutine nsi_multi_step_fs_allocate

  subroutine nsi_multi_step_fs_eval_c_mass()

    use mod_solver,                      only : solver_initialize_matrix_and_rhs
    use mod_solver,                      only : solver_periodicity
    use mod_matrix,                      only : matrix_copy_matrix_to_block
    use mod_nsi_dirichlet_global_system, only : nsi_dirichlet_momentum
    use mod_nsi_dirichlet_global_system, only : nsi_dirichlet_rotate
    
    implicit none

    solve_sol => solve(1:)
    
    if( INOTMASTER ) then
       !
       ! Copy weighted mass matrix
       !
       cmasm(1:solve(1) % nzmat) = 0.0_rp
       masbc = 0.0_rp
       call matrix_copy_matrix_to_block(npoin,ndime,r_dom,c_dom,cmass_weighted,cmasm)
       !
       ! Periodicity
       !
       call solver_periodicity('MATRIX',solve_sol,ndime,ndime,cmasm)!,masbc)
       !
       ! Rotate matrix if necessary
       !
       call nsi_dirichlet_rotate(cmasm,masbc)
       !
       ! Impose boundary conditions
       !
       call nsi_dirichlet_momentum(Auu=cmasm,bu=masbc)
       
    end if
    
  end subroutine nsi_multi_step_fs_eval_c_mass

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-17
  !> @brief   Compute RHS when using consistent matrix
  !> @details Compute RHS and impose boundary conditions,
  !>          already caclulated
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_multi_step_fs_rhs(uu01,my_sign)
    
    real(rp), intent(inout) :: uu01(*)
    real(rp), intent(in)    :: my_sign    
    integer(ip)             :: idime,ipoin,idofn
    !
    ! Rotate if necessary
    ! 
    call nsi_rotsch(NSI_GLOBAL_TO_LOCAL,rt_tmp) ! Local to global
    !
    ! r = r / dt + M u^{k-1} 
    !
    call nsi_auuvec(1_ip,cmasm,uu01,xx)
    if( my_sign > 0.0_rp ) then
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn         = (ipoin-1)*ndime + idime
             rt_tmp(idofn) = rt_tmp(idofn) / dtinv_nsi + xx(idofn)
          end do
       end do
    else
       do ipoin = 1,npoin
          do idime = 1,ndime
             idofn         = (ipoin-1)*ndime + idime
             rt_tmp(idofn) = - rt_tmp(idofn) / dtinv_nsi + xx(idofn)
          end do
       end do       
    end  if
    !
    ! Periodicity
    !
    !call solver_periodicity_rhs(ndime,rt_tmp)
    !
    ! Boundary and initial conditions
    !
    do ipoin = 1,npoin
       do idime = 1,ndime
          idofn = (ipoin-1)*ndime + idime
          uu(idofn) = uu01(idofn)
          if( kfl_fixno_nsi(idime,ipoin) > 0 ) then
             rt_tmp(idofn) = masbc(idofn)
          end if
       end do
    end do

  end subroutine nsi_multi_step_fs_rhs
  !
  ! Implements : Explicit Runge–Kutta schemes for incompressible flow with improved energy-conservation properties
  !              F.Capuano, G.Coppola, L.Rández, L.deLuca, Journal ofComputationalPhysics328(2017)86–94
  !
  ! Schemme 3p5q(4) from the paper 3o time, 5o energy and 4 steps
  !
  subroutine nsi_multi_step_fs_solution()
     integer(ip)       :: idofn,ipoin,idime

     !-----------------------------------------------------------------
     !
     ! Allocate memory if necessary
     !
     !-----------------------------------------------------------------

     call nsi_multi_step_fs_allocate()

     ! guardar la uu de n
     if( INOTMASTER ) then
        do ipoin = 1,npoin
           do idime = 1,ndime
              idofn = (ipoin-1)*ndime + idime
              uu0(idofn)   = uu(idofn)
           end do
        end do
     end if

     if(kfl_tiacc_nsi == 2) then
        call nsi_multi_step_fs_eval(3_ip)
        call nsi_multi_step_fs_eval(4_ip)
     else if(kfl_tiacc_nsi == 3) then  
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
     else ! order 4
        call nsi_multi_step_fs_eval(1_ip)
        call nsi_multi_step_fs_eval(2_ip)
        call nsi_multi_step_fs_eval(3_ip)
        call nsi_multi_step_fs_eval(4_ip)
     end if


  end subroutine nsi_multi_step_fs_solution

  subroutine nsi_multi_step_fs_eval(istep)

     integer(ip) , intent(in) :: istep
     integer(ip)              :: idofn,ipoin,idime
     real(rp)                 :: timea,timeb


     !-----------------------------------------------------------------
     !
     ! Assemble equations
     !
     !-----------------------------------------------------------------

     call nsi_solini(3_ip)                                        ! Initialize solver
     if( kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM ) then
        call nsi_rotsch(NSI_GLOBAL_TO_LOCAL,uu)                   ! Global to local
        call nsi_rotsch(NSI_GLOBAL_TO_LOCAL,uu0)                  ! Global to local
     end if

     call nsi_matrix()                                            ! Assemble equation

     !do ipoin = 1,npoin
     !   write(90,*) ipoin,': ',Q(r_dom(ipoin):r_dom(ipoin+1)-1)
     !end do
     !call flush(90)
     !stop

     if( kfl_massm_nsi == NSI_CONSISTENT_MASS ) then
        call nsi_multi_step_fs_eval_c_mass()
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
     ! Go back to global system
     !
     !-----------------------------------------------------------------

     call livinf(165_ip,'C',0_ip) 
     if( kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM ) then
        call nsi_rotsch(NSI_LOCAL_TO_GLOBAL,uu) ! Local to global
     end if
     do ipoin = 1,npoin
        do idime = 1,ndime
           idofn = (ipoin-1)*ndime + idime
           veloc(idime,ipoin,1) = uu(idofn)  
        end do
     end do
     if( kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM ) then
        call nsi_rotsch(NSI_LOCAL_TO_GLOBAL,uu0) ! Local to global
     end if

     call livinf(164_ip,' ',1_ip) 
     call cputim(timeb)
     cpu_ass_sol_nsi(2) = timeb - timea

  end subroutine nsi_multi_step_fs_eval

  subroutine nsi_multi_step_fs_solution_sj(istep)
     integer(ip) , intent(in) :: istep
     integer(ip)       :: idofn,ipoin,idime,jstep
     real(rp) :: aux2
     real(rp) :: time1,time2

     !-----------------------------------------------------------------
     !
     ! Solve Momentum
     !
     !-----------------------------------------------------------------

     call cputim(time1)
     if( INOTMASTER ) then
        !
        ! Momentum residual rm^k
        !
        if(kfl_fscon_nsi == 0) then 
           if(ittim>n_moin_steps) then 
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
                 aux(ipoin) = pp(ipoin)
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
        
        if( kfl_grad_div_nsi /= 0 )  then
           call nsi_aupvec(1_ip,Grad,aux,vv)                          
        else
           call nsi_aupvec(1_ip,amatr(poaup_nsi:),aux,vv)         
        end if
        
        if( kfl_massm_nsi == NSI_LUMPED_MASS  ) then
           call rhsmod(ndime,rm(:,istep))
           call rhsmod(ndime,vv)
        else
           call solver_periodicity_rhs(ndime,rm(:,istep))
           call solver_periodicity_rhs(ndime,vv)
        end if
             
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

     if( kfl_massm_nsi == NSI_CONSISTENT_MASS ) then
        !
        ! Consistent mass
        !
        solve_sol => solve(1:)
        call nsi_rotsch(NSI_GLOBAL_TO_LOCAL,uu0)                      ! Global to local
        call nsi_multi_step_fs_rhs(uu0,1.0_rp)                        ! Compute RHS
        call solver(rt_tmp,uu,cmasm,pmatr)                            ! Solve system
        call nsi_rotsch(NSI_LOCAL_TO_GLOBAL,uu)                       ! Local to global
        call nsi_rotsch(NSI_LOCAL_TO_GLOBAL,uu0)                      ! Local to global
         
     else
        !
        ! Lumped mass
        !
        call solver_lumped_mass_system(ndime,rt_tmp,EXCHANGE=.false.) ! Lumped system
        do ipoin = 1,npoin
           do idime = 1,ndime
              idofn     = (ipoin-1)*ndime + idime
              uu(idofn) = uu0(idofn) + rt_tmp(idofn)*Tim(ipoin)       ! u = u0 + r * dt/rho
           end do
        end do
        call nsi_multi_step_fs_dirichlet()                            ! Prescribe Drichlet bc on  u'^{k+1}
     end if

     call livinf(165_ip,'M',0_ip)
     call cputim(time2)
     cputi_nsi(1) = cputi_nsi(1) + time2 - time1

     if( kfl_fscon_nsi == 1 ) then 
        call nsi_multi_step_fs_pressure_and_correction(a(istep+1),0_ip)
     end if

  end subroutine nsi_multi_step_fs_solution_sj

  subroutine nsi_multi_step_fs_solution_final()
     integer(ip)       :: idofn,ipoin,idime
     real(rp)          :: aux
     real(rp) :: time1,time2

     !-----------------------------------------------------------------
     !
     ! Solve Momentum
     !
     !-----------------------------------------------------------------
     call cputim(time1)

     if( INOTMASTER ) then
        !
        ! Momentum residual rm^k
        !
        do idofn = 1,ndime*npoin                                        ! rm = bu 
           rm(idofn,4) = bu(idofn)  
        end do
        if( kfl_massm_nsi == NSI_LUMPED_MASS ) then
           call rhsmod(ndime,rm(:,4))
        else
           call solver_periodicity_rhs(ndime,rm(:,4))
        end if

        do ipoin = 1,npoin
           do idime = 1,ndime
              idofn     = (ipoin-1)*ndime + idime
              aux       = rm(idofn,1)*sigma(1) + rm(idofn,2)*sigma(2) + rm(idofn,3)*sigma(3) + rm(idofn,4)*sigma(4)
              rt(idofn) = aux
              !uu(idofn) = uu0(idofn) + (aux)*Tim(ipoin)/Mass(ipoin)
              !rt_tmp(idofn) = (aux)*Tim(ipoin)
              rt_tmp(idofn) = (aux)
           end do
        end do
     end if
          
     if( kfl_massm_nsi == NSI_CONSISTENT_MASS ) then
        !
        ! Consistent mass
        !
        solve_sol => solve(1:)
        call nsi_rotsch(NSI_GLOBAL_TO_LOCAL,uu0)                      ! Global to local
        call nsi_multi_step_fs_rhs(uu0,1.0_rp)
        call solver(rt_tmp,uu,cmasm,pmatr)                            ! Solve system
        call nsi_rotsch(NSI_LOCAL_TO_GLOBAL,uu)                       ! Local to global
        call nsi_rotsch(NSI_LOCAL_TO_GLOBAL,uu0)                      ! Local to global
     else
        !
        ! Lumped mass
        !
        call solver_lumped_mass_system(ndime,rt_tmp,EXCHANGE=.false.) ! CAMBIO ORIOL
        do ipoin = 1,npoin
           do idime = 1,ndime
              idofn     = (ipoin-1)*ndime + idime
              uu(idofn) = uu0(idofn) + rt_tmp(idofn)*Tim(ipoin)
           end do
        end do
        call nsi_multi_step_fs_dirichlet()                            ! Prescribe Dirichlet bc on u'^{k+1}
     end if
     
     call livinf(165_ip,'M',0_ip)
     call cputim(time2)
     cputi_nsi(1) = cputi_nsi(1) + time2 - time1

     call nsi_multi_step_fs_pressure_and_correction(1.0_rp,1_ip)

  end subroutine nsi_multi_step_fs_solution_final

  subroutine nsi_multi_step_fs_dirichlet(xarray)

     real(rp),   intent(out), optional :: xarray(:)
     integer(ip)                       :: idofn,ipoin,idime

     if( kfl_matdi_nsi == NSI_DIRICHLET_ALGORITHM ) then

        if( present(xarray) ) then
           call nsi_rotsch(NSI_GLOBAL_TO_LOCAL,xarray) ! Global to local
           idofn = 0
           do ipoin = 1,npoin
              do idime = 1,ndime
                 idofn = idofn + 1
                 if( kfl_fixno_nsi(idime,ipoin) > 0 .and. kfl_noslw_ker /= 0) &
                    !                     write(1000*idime+kfl_paral,*),xarray(idofn) 
                 vafor_nsi(idime,ipoin) = xarray(idofn) ! actualy in the multistep case we could think of averaging all the substeps
                 xarray(idofn) = 0.0_rp
              end do
           end do
           call nsi_rotsch(NSI_LOCAL_TO_GLOBAL,xarray) ! Local to global
        else
           call nsi_rotsch(NSI_GLOBAL_TO_LOCAL,uu) ! Global to local
           idofn = 0
           do ipoin = 1,npoin
              do idime = 1,ndime
                 idofn = idofn + 1
                 if( kfl_fixno_nsi(idime,ipoin) > 0 ) &
                    uu(idofn) = bvess_nsi(idime,ipoin,1)
              end do
           end do
           call nsi_rotsch(NSI_LOCAL_TO_GLOBAL,uu) ! Local to global
        end if

     end if

  end subroutine nsi_multi_step_fs_dirichlet

  subroutine nsi_multi_step_fs_pressure_and_correction(factor,eval_rc)
    
     real(rp),    intent(in) :: factor
     integer(ip), intent(in) :: eval_rc
     integer(ip)             :: idofn,ipoin,idime
     real(rp)                :: aux2
     real(rp)                :: time1,time2

     !-----------------------------------------------------------------
     !
     ! Solve pressure
     !
     !-----------------------------------------------------------------

     call cputim(time1)
     if( INOTMASTER ) then
        !
        ! Continuity residual rc'^{k+1}
        !
        if( kfl_grad_div_nsi /= 0 ) then
           call nsi_apuvec(1_ip,Div,uu,divup)
        else
           call nsi_apuvec(1_ip,amatr(poapu_nsi:),uu,divup)             ! rc' =  Apu u*
        end if
        do ipoin = 1,npoin                                               
           rc(ipoin) = (bp(ipoin) - divup(ipoin))/(factor*Tim(ipoin))
        end do
        call nsi_appvec(0_ip,Q,pp,rc,-1.0_rp)                           ! rc <= rc - Q p^k   
        Dp = 0.0_rp

     end if
     !
     ! Periodicity (only rc is modified, not the initial guess Dp - then the solver obiouslly behaves diferently when GRAD_DIV is added)
     ! I tried to add call solver_periodicity('SOLUTION'..)  but it introduced no advantage
     !
     !if( kfl_grad_div_nsi /= 0 .and. INOTMASTER ) call solver_periodicity('RHS',solve(2:),1_ip,1_ip,Q,rc,Dp)
     !
     ! Solve: Q z = rc'^k
     !
     call nsi_solini(2_ip)                                              ! Initialize solver

     call solver(rc,Dp,Q,pmatr)                                         ! Solve system Q Dp = rc^k

     if( INOTMASTER ) then
        !
        ! Dp      = Dp' + p^k
        ! p^{k+1} = Dp   p^k
        !
        do ipoin = 1,npoin
           Dp(ipoin) = Dp(ipoin) + pp(ipoin)  
           pp(ipoin) = Dp(ipoin) - fsrot_nsi*Nu(ipoin)*divup(ipoin)
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
        if( kfl_massm_nsi == NSI_LUMPED_MASS ) then
           call rhsmod(ndime,vv)
        else
           call solver_periodicity_rhs(ndime,vv)
        end if

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
              rt_tmp(idofn) = factor*vv(idofn)
              rt(idofn)     = rt(idofn) - vv(idofn) 
           end do
        end do
     end if
   
     if( kfl_massm_nsi == NSI_CONSISTENT_MASS ) then
        !
        ! Consisten mass
        !
        solve_sol => solve(1:)
        call nsi_rotsch(NSI_GLOBAL_TO_LOCAL,uu)                       ! Global to local
        call nsi_multi_step_fs_rhs(uu,-1.0_rp)                        ! Compute RHS
        call solver(rt_tmp,uu,cmasm,pmatr)                            ! Solve system
        call nsi_rotsch(NSI_LOCAL_TO_GLOBAL,uu)                       ! Local to global
     else
        call solver_lumped_mass_system(ndime,rt_tmp,EXCHANGE=.false.) ! CAMBIO ORIOL
        do ipoin = 1,npoin
           do idime = 1,ndime
              idofn     = (ipoin-1)*ndime + idime
              uu(idofn) = uu(idofn) - rt_tmp(idofn)*Tim(ipoin)
           end do
        end do
        call nsi_multi_step_fs_dirichlet()                            ! Prescribe bc 
     end if

     if( kfl_matdi_nsi == NSI_DIRICHLET_ALGORITHM ) then
        if( INOTMASTER ) rt_tmp = rt
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



end module mod_nsi_multi_step_fs
