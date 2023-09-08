!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_dirichlet_global_system
!> @author  Guillaume Houzeaux
!> @date    10/10/2016
!> @brief   Dirichlet conditions on Navier-Stokes
!> @details Impose Dirichlet conditions in global system
!>          For fractional, Auu and App are not assembled.
!>
!------------------------------------------------------------------------

module mod_nsi_dirichlet_global_system

  use def_kintyp
  use def_domain
  use def_master
  use def_nastin
  use mod_matrix, only : matrix_rotate_system

  implicit none
  
  private

  public :: nsi_dirichlet_global_system  
  public :: nsi_dirichlet_laplacian
  public :: nsi_dirichlet_momentum
  public :: nsi_dirichlet_rotate
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-17
  !> @brief   Impose Dirichlet on global system
  !> @details Impose Dirichlet boundary conditions on all the matrices
  !>          of the global system
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_dirichlet_global_system(Aup,Apu,Q,bu,bp,uu,pp,Auu,App)

    real(rp),    intent(out)           :: Aup(ndime,nzdom)
    real(rp),    intent(out)           :: Apu(ndime,nzdom)
    real(rp),    intent(out)           :: Q(nzdom)
    real(rp),    intent(out)           :: bu(ndime,npoin)
    real(rp),    intent(out)           :: bp(npoin)
    real(rp),    intent(out)           :: uu(ndime,npoin)
    real(rp),    intent(out)           :: pp(npoin)
    real(rp),    intent(out), optional :: Auu(ndime,ndime,nzdom)
    real(rp),    intent(out), optional :: App(nzdom)

    if( kfl_matdi_nsi == NSI_DIRICHLET_MATRIX .and. INOTMASTER ) then
       if( NSI_FRACTIONAL_STEP ) then
          call nsi_dirichlet_global_system_go(Aup,Apu,Q,bu,bp,uu,pp)
       else
          call nsi_dirichlet_global_system_go(Aup,Apu,Q,bu,bp,uu,pp,Auu,App)
       end if
    end if

  end subroutine nsi_dirichlet_global_system

  subroutine nsi_dirichlet_global_system_go(Aup,Apu,Q,bu,bp,uu,pp,Auu,App)

    real(rp),    intent(out)           :: Aup(ndime,nzdom)
    real(rp),    intent(out)           :: Apu(ndime,nzdom)
    real(rp),    intent(out)           :: Q(nzdom)
    real(rp),    intent(out)           :: bu(ndime,npoin)
    real(rp),    intent(out)           :: bp(npoin)
    real(rp),    intent(out)           :: uu(ndime,npoin)
    real(rp),    intent(out)           :: pp(npoin)
    real(rp),    intent(out), optional :: Auu(ndime,ndime,nzdom)
    real(rp),    intent(out), optional :: App(nzdom)
    real(rp)                           :: Auud(3),Qd,Appd
    integer(ip)                        :: ipoin,jzdom,idime,jdime,izdom,jpoin
    integer(ip)                        :: izdod,ibopo,jbopo,kpoin,iroty,kdime
    integer(ip)                        :: kfl_addma,idofn
    integer(4)                         :: istat
    real(rp)                           :: worma(ndime,ndime)
    real(rp),    pointer               :: rotma(:,:)


    !----------------------------------------------------------------------
    !
    ! Save matrix to compute internal force
    !
    !----------------------------------------------------------------------

    if( kfl_intfo_nsi >= 1 ) then

       if( kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM ) then    ! in the DIRICHLET_ALGORITHM case I am saving directly vafor_nsi but I am not sure what happens with skews

          if( NSI_FRACTIONAL_STEP ) call runend('PROGRAMALO VAGO')
          ! if( kfl_local_nsi == 1 ) then
          !    call runend('INTERNAL FORCE NOT CODED FOR LOCAL AXES')
          ! end if
          do ipoin = 1,npoin
             if( lpoty(ipoin) > 0 ) then
                if( kfl_intfo_nsi == 2 ) then
                   if( present(Auu) ) deallocate( intfo_nsi(ipoin) % Auu , stat = istat )
                   deallocate( intfo_nsi(ipoin) % Aup , stat = istat )
                   deallocate( intfo_nsi(ipoin) % bu  , stat = istat )
                end if
                jzdom = r_dom(ipoin+1)-r_dom(ipoin)
                if( present(Auu) ) allocate( intfo_nsi(ipoin) % Auu(ndime,ndime,jzdom) , stat = istat )
                allocate( intfo_nsi(ipoin) % Aup(ndime,jzdom)       , stat = istat )
                allocate( intfo_nsi(ipoin) % bu(ndime)              , stat = istat )
                jzdom = 0
                do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                   jzdom = jzdom + 1
                   if( present(Auu) ) then
                      do idime = 1,ndime
                         do jdime = 1,ndime
                            intfo_nsi(ipoin) % Auu(jdime,idime,jzdom) = Auu(jdime,idime,izdom)
                         end do
                      end do
                   end if
                   do idime = 1,ndime
                      intfo_nsi(ipoin) % Aup(idime,jzdom) = Aup(idime,izdom)
                      intfo_nsi(ipoin) % bu(idime)        = bu(idime,ipoin)
                   end do
                end do
             end if
          end do
       end if
    end if

    !----------------------------------------------------------------------
    !
    ! Rotate matrix
    !
    !----------------------------------------------------------------------
    
    if( kfl_local_nsi == 1 ) then

       do ipoin = 1,npoin
          ibopo = lpoty(ipoin)
          if( ibopo > 0 ) then
             iroty =  kfl_fixrs_nsi(ibopo)
             if( iroty /= 0 ) then
                if( iroty == -1 ) then                                    ! Tangent system
                   rotma => exnor(:,:,ibopo)
                else if( iroty >= 1 ) then                                ! Given system
                   rotma =>  skcos(:,:,iroty)
                else if( iroty == -2 ) then                               ! Given system
                   rotma => skcos_nsi(:,:,ibopo)
                else if( iroty == -3 ) then                               ! Geometrical normal
                   rotma => skcos(:,:,ibopo)
                end if
                !
                ! Modifies column number IPOIN of AMATR ( A_j,imodi <-- A_j,imodi R )
                !
                do jpoin = 1,npoin
                   do izdom = r_dom(jpoin),r_dom(jpoin+1)-1
                      kpoin = c_dom(izdom)
                      if( kpoin == ipoin ) then
                         if( present(Auu) ) then
                            do idime = 1,ndime
                               do jdime = 1,ndime
                                  worma(idime,jdime) = 0.0_rp
                                  do kdime = 1,ndime
                                     worma(idime,jdime) = worma(idime,jdime) &
                                          + Auu(kdime,idime,izdom) * rotma(kdime,jdime)
                                  end do
                               end do
                            end do
                            do idime = 1,ndime
                               do jdime = 1,ndime
                                  Auu(jdime,idime,izdom) = worma(idime,jdime)
                               end do
                            end do
                         end if

                         do jdime = 1,ndime
                            worma(1,jdime) = 0.0_rp
                            do kdime = 1,ndime
                               worma(1,jdime) = worma(1,jdime)&
                                    + Apu(kdime,izdom) * rotma(kdime,jdime)
                            end do
                         end do
                         do jdime = 1,ndime
                            Apu(jdime,izdom) = worma(1,jdime)
                         end do

                      end if
                   end do
                end do

                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   !
                   ! Modifies row number IPOIN of AMATR ( A_imodi,j <-- R^t A_imodi,j )
                   !
                   if( present(Auu) ) then
                      jpoin = c_dom(izdom)
                      do idime = 1,ndime
                         do jdime = 1,ndime
                            worma(idime,jdime) = 0.0_rp
                            do kdime = 1,ndime
                               worma(idime,jdime) = worma(idime,jdime) &
                                    + Auu(jdime,kdime,izdom) * rotma(kdime,idime)
                            end do
                         end do
                      end do
                      do idime = 1,ndime
                         do jdime = 1,ndime
                            Auu(jdime,idime,izdom) = worma(idime,jdime)
                         end do
                      end do
                   end if
                   !
                   ! Modify the part corresponding to a scalar unknown
                   !
                   do idime = 1,ndime
                      worma(idime,1) = 0.0_rp
                      do kdime = 1,ndime
                         worma(idime,1) = worma(idime,1) &
                              + rotma(kdime,idime) * Aup(kdime,izdom)
                      end do
                   end do
                   do idime = 1,ndime
                      Aup(idime,izdom) = worma(idime,1)
                   end do

                end do
                !
                ! Rotate RHS: bu
                !
                do idime = 1,ndime
                   worma(idime,1) = 0.0_rp
                   do kdime = 1,ndime
                      worma(idime,1) = worma(idime,1) &
                           + rotma(kdime,idime) * bu(kdime,ipoin)
                   end do
                end do
                do idime = 1,ndime
                   bu(idime,ipoin) = worma(idime,1)
                end do

             end if
          end if
       end do
    end if

    !----------------------------------------------------------------------
    !
    ! Impose velocity
    !
    !----------------------------------------------------------------------

    do ipoin = 1,npoin

       do idime = 1,ndime

          if( kfl_fixno_nsi(idime,ipoin) > 0 ) then
             !
             ! Eliminate dof of IPOIN from other equations (JPOIN)
             !
             do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                jpoin = c_dom(izdom)
                if( ipoin /= jpoin ) then

                   do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                      kpoin = c_dom(jzdom)
                      if( kpoin == ipoin ) then
                         if( present(Auu) ) then
                            do jdime = 1,ndime
                               bu(jdime,jpoin) = bu(jdime,jpoin) - Auu(idime,jdime,jzdom) * bvess_nsi(idime,ipoin,1)
                               Auu(idime,jdime,jzdom) = 0.0_rp
                            end do
                         end if
                         if( kfl_grad_div_nsi == 0 ) then
                            bp(jpoin) = bp(jpoin) - Apu(idime,jzdom) * bvess_nsi(idime,ipoin,1)
                            Apu(idime,jzdom) = 0.0_rp
                         end if
                      end if
                   end do

                end if
             end do
             !
             ! IZDOD: Diagonal
             !
             izdod = r_dom(ipoin) - 1
             jpoin = 0
             do while( jpoin /= ipoin )
                izdod = izdod + 1
                jpoin = c_dom(izdod)
             end do
             if( present(Auu) ) then
                Auud(idime) = Auu(idime,idime,izdod)
             else
                Auud(idime) = 1.0_rp
             end if
             if( abs(Auud(idime)) < zeror ) Auud(idime) = 1.0_rp
             !
             ! Set line to zero
             !
             do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                if( present(Auu) ) then
                   do jdime = 1,ndime
                      Auu(jdime,idime,izdom) = 0.0_rp
                   end do
                end if
                Aup(idime,izdom) = 0.0_rp
             end do
             !
             ! Prescribe value
             !
             idofn                  = (ipoin-1)*ndime + idime
             if( present(Auu) ) then
                Auu(idime,idime,izdod) = Auud(idime)
                bu(idime,ipoin)        = bvess_nsi(idime,ipoin,1) * Auud(idime)
                uu(idime,ipoin)        = bvess_nsi(idime,ipoin,1)
             else
                !
                ! In the case when Auu is eliminated (bu <= bu - Auu*u), bu
                ! should be zero on boundary condition
                !
                bu(idime,ipoin)        = 0.0_rp
                uu(idime,ipoin)        = bvess_nsi(idime,ipoin,1)
             end if

          end if

       end do

    end do

    !----------------------------------------------------------------------
    !
    ! Impose pressure Schur complement preconditioner
    !
    !----------------------------------------------------------------------

    kfl_addma = 0

    if( solve(2) % kfl_symme == 1 ) then

       do ipoin = 1,npoin
          do izdom = r_sym(ipoin),r_sym(ipoin+1) - 2
             jpoin = c_sym(izdom)
             jbopo = lpoty(jpoin)
             if( jbopo /= 0 ) then
                if( kfl_fixpr_nsi(1,jpoin) > 0 ) then
                   Q(izdom) = 0.0_rp
                end if
             end if
          end do
       end do

       do ipoin = 1,npoin

          ibopo = lpoty(ipoin)
          if( ibopo > 0 ) then

             if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                !
                ! IZDOD: Diagonal
                !
                izdod = r_sym(ipoin+1) - 1
                Qd = Q(izdod)
                if( abs(Qd) < zeror ) Qd = 1.0_rp
                !
                ! Set line to zero
                !
                do izdom = r_sym(ipoin),r_sym(ipoin+1) - 1
                   Q(izdom) = 0.0_rp
                end do
                !
                ! Presrcibe value
                !
                Q(izdod) = Qd

             end if

          end if

       end do

    else

       if( kfl_addma == 0 ) then
          do ipoin = 1,npoin
             do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                jpoin = c_dom(izdom)
                if( ipoin /= jpoin ) then
                   jbopo = lpoty(jpoin)
                   if( jbopo /= 0 ) then
                      if( kfl_fixpr_nsi(1,jpoin) > 0 ) then
                         Q(izdom) = 0.0_rp
                      end if
                   end if
                end if
             end do
          end do
       end if

       do ipoin = 1,npoin

          ibopo = lpoty(ipoin)
          if( ibopo > 0 ) then

             if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                !
                ! IZDOD: Diagonal
                !
                izdod = r_dom(ipoin) - 1
                jpoin = 0
                do while( jpoin /= ipoin )
                   izdod = izdod + 1
                   jpoin = c_dom(izdod)
                end do
                Qd = Q(izdod)
                if( abs(Qd) < zeror ) Qd = 1.0_rp

                if( kfl_addma == 0 ) then
                   !
                   ! Set line to zero
                   !
                   do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                      Q(izdom) = 0.0_rp
                   end do
                   !
                   ! Presrcibe value
                   !
                   Q(izdod) = Qd
                else
                   !
                   ! Put mass
                   !
                   do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                      !Q(izdom) = 0.0_rp
                   end do
                   call runend('NSI_MATDIR: NOT CODED')
!!!!Q(izdod) =  vmass(ipoin) / ( visco_nsi(1,1) + turmu(ipoin) )
                   if( ISLAVE ) call runend('MUST USE NOT EXCHANGED MASS')
                end if

             end if

          end if
       end do

    end if

    !----------------------------------------------------------------------
    !
    ! Impose pressure
    !
    !----------------------------------------------------------------------

    !if ( kfl_confi_nsi == 1 .and. nodpr_nsi > 0 ) then

    do ipoin = 1,npoin

       if( solve(1) % block_array(2) % kfl_fixno(1,ipoin) > 0 ) then

          if( solve(2) % kfl_symme == 1 ) then
             call runend('NOT CODED: CHECK IT')
             !
             ! Eliminate pressure at IPOIN on all lines
             !
             do jpoin = 1,npoin
                if( ipoin /= jpoin ) then
                   izdom = r_sym(jpoin)
                   do while( izdom < r_sym(jpoin+1) )
                      kpoin = c_sym(izdom)
                      if( kpoin == ipoin ) then
                         if( present(App) ) then
                            bp(jpoin)  = bp(jpoin) - App(izdom) * solve(1) % block_array(2) % bvess(1,ipoin)
                            App(izdom) = 0.0_rp
                         end if
                         Q(izdom)   = 0.0_rp
                         izdom      = r_sym(jpoin+1)
                      end if
                      izdom = izdom + 1
                   end do
                end if
             end do
             !
             ! Diagonal
             !
             izdod = r_sym(ipoin+1) - 1
             Qd    = Q(izdod)
             if( present(App) ) then
                Appd = App(izdod)
             else
                Appd = 0.0_rp
             end if
             if( abs(Qd)   < zeror )   Qd = 1.0_rp
             if( abs(Appd) < zeror ) Appd = 1.0_rp
             !
             ! Set line to zero
             !
             do izdom = r_sym(ipoin),r_sym(ipoin+1) - 1
                Q(izdom)   = 0.0_rp
                if( present(App) ) App(izdom) = 0.0_rp
             end do
             do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                do idime = 1,ndime
                   Apu(idime,izdom) = 0.0_rp
                end do
             end do
             !
             ! Prescribe value
             !
             Q(izdod)   = Qd
             if( present(App) ) App(izdod) = Appd
             bp(ipoin)  = Appd * solve(1) % block_array(2) % bvess(1,ipoin)
             pp(ipoin)  = solve(1) % block_array(2) % bvess(1,ipoin)

          else
             !
             ! Eliminate pressure at IPOIN on all lines
             !
             do jpoin = 1,npoin
                if( ipoin /= jpoin ) then
                   izdom = r_dom(jpoin)
                   do while( izdom < r_dom(jpoin+1) )
                      kpoin = c_dom(izdom)
                      if( kpoin == ipoin ) then
                         if( present(App) ) then
                            bp(jpoin)  = bp(jpoin) - App(izdom) * solve(1) % block_array(2) % bvess(1,ipoin)
                            App(izdom) = 0.0_rp
                         end if
                         Q(izdom)   = 0.0_rp
                         izdom      = r_dom(jpoin+1)
                      end if
                      izdom = izdom + 1
                   end do
                end if
             end do
             !
             ! Diagonal
             !
             izdod = r_dom(ipoin) - 1
             jpoin = 0
             do while( jpoin /= ipoin )
                izdod = izdod + 1
                jpoin = c_dom(izdod)
             end do
             if( present(App) ) then
                Appd = App(izdod)
             else
                Appd = 0.0_rp
             end if
             Qd   = Q(izdod)
             if( abs(Appd) < zeror ) Appd = 1.0_rp
             !
             ! Set line to zero
             !
             do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                if( present(App) ) App(izdom) = 0.0_rp
                Q(izdom)   = 0.0_rp
                do idime = 1,ndime
                   Apu(idime,izdom) = 0.0_rp
                end do
             end do
             !
             ! Presrcibe value
             !
             Q(izdod)   = Qd
             if( present(App) ) App(izdod) = Appd
             bp(ipoin)  = Appd * solve(1) % block_array(2) % bvess(1,ipoin)
             pp(ipoin)  = solve(1) % block_array(2) % bvess(1,ipoin)

          end if
          !
          ! Eliminate prescribed pressure from IPOIN momentum equation
          !
          do jpoin = 1,npoin
             izdom = r_dom(jpoin)
             do while( izdom < r_dom(jpoin+1) )
                kpoin = c_dom(izdom)
                if( kpoin == ipoin ) then
                   do idime = 1,ndime
                      bu(idime,jpoin)  = bu(idime,jpoin) - Aup(idime,izdom) * solve(1) % block_array(2) % bvess(1,ipoin)
                      Aup(idime,izdom) = 0.0_rp
                   end do
                   izdom = r_dom(jpoin+1)
                end if
                izdom = izdom + 1
             end do
          end do

       end if

    end do

    !end if


  end subroutine nsi_dirichlet_global_system_go

  subroutine nsi_dirichlet_laplacian(Q,bu)

    use mod_solver, only : solver_periodicity
    
    real(rp),    intent(out)           :: Q(nzdom)
    real(rp),    intent(in), optional  :: bu(ndime,*)
    real(rp)                           :: Auud(3),Qd,Appd
    integer(ip)                        :: ipoin,jzdom,izdom,jpoin,idime
    integer(ip)                        :: izdod,ibopo,jbopo,kpoin
    real(rp)                           :: dummu(2),dummb(2)
    
    !----------------------------------------------------------------------
    !
    ! Save matrix to compute internal force
    !
    !----------------------------------------------------------------------

    if( kfl_intfo_nsi >= 1 .and. present(bu) ) then

       if( kfl_matdi_nsi /= NSI_DIRICHLET_ALGORITHM ) then ! in the DIRICHLET_ALGORITHM case I am saving directly vafor_nsi
          if( NSI_FRACTIONAL_STEP ) call runend('PROGRAMALO VAGO')
          do ipoin = 1,npoin
             if( lpoty(ipoin) > 0 ) then
                if( kfl_intfo_nsi == 2 ) then
                   deallocate( intfo_nsi(ipoin) % bu  )
                end if
                allocate( intfo_nsi(ipoin) % bu(ndime) )
                do idime = 1,ndime
                   intfo_nsi(ipoin) % bu(idime) = bu(idime,ipoin)
                end do
             end if
          end do
       end if
    end if

    !----------------------------------------------------------------------
    !
    ! Impose pressure Schur complement preconditioner
    !
    !----------------------------------------------------------------------

    if( kfl_grad_div_nsi == 1 ) then

       !call solver_preprocess(solve,Q,b,u)
       call solver_periodicity('MATRIX',solve(2:),1_ip,1_ip,Q,dummb,dummu,'KEEP DIAGONAL')
       
       if( solve(2) % kfl_symme == 1 ) then

          do ipoin = 1,npoin
             do izdom = r_sym(ipoin),r_sym(ipoin+1) - 2
                jpoin = c_sym(izdom)
                jbopo = lpoty(jpoin)
                if( jbopo /= 0 ) then
                   if( kfl_fixpr_nsi(1,jpoin) > 0 ) then
                      Q(izdom) = 0.0_rp
                   end if
                end if
             end do
          end do

          do ipoin = 1,npoin

             ibopo = lpoty(ipoin)
             if( ibopo > 0 ) then

                if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                   !
                   ! IZDOD: Diagonal
                   !
                   izdod = r_sym(ipoin+1) - 1
                   Qd = Q(izdod)
                   if( abs(Qd) < zeror ) Qd = 1.0_rp
                   !
                   ! Set line to zero
                   !
                   do izdom = r_sym(ipoin),r_sym(ipoin+1) - 1
                      Q(izdom) = 0.0_rp
                   end do
                   !
                   ! Presrcibe value
                   !
                   Q(izdod) = Qd

                end if

             end if

          end do

       else

          do ipoin = 1,npoin
             do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                jpoin = c_dom(izdom)
                if( ipoin /= jpoin ) then
                   jbopo = lpoty(jpoin)
                   if( jbopo /= 0 ) then
                      if( kfl_fixpr_nsi(1,jpoin) > 0 ) then
                         Q(izdom) = 0.0_rp
                      end if
                   end if
                end if
             end do
          end do

          do ipoin = 1,npoin

             ibopo = lpoty(ipoin)
             if( ibopo > 0 ) then

                if( kfl_fixpr_nsi(1,ipoin) > 0 ) then
                   !
                   ! IZDOD: Diagonal
                   !
                   izdod = r_dom(ipoin) - 1
                   jpoin = 0
                   do while( jpoin /= ipoin )
                      izdod = izdod + 1
                      jpoin = c_dom(izdod)
                   end do
                   Qd = Q(izdod)
                   if( abs(Qd) < zeror ) Qd = 1.0_rp
                   !
                   ! Set line to zero
                   !
                   do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                      Q(izdom) = 0.0_rp
                   end do
                   !
                   ! Presrcibe value
                   !
                   Q(izdod) = Qd

                end if

             end if
          end do

       end if

    end if

  end subroutine nsi_dirichlet_laplacian

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-17
  !> @brief   Impose bc in momentum equations
  !> @details Impose bc in momentum equations, usefule for FS technique
  !> 
  !-----------------------------------------------------------------------

  subroutine nsi_dirichlet_momentum(Auu,bu,uu,UNIT_DIAGONAL)

    real(rp),    optional, intent(inout) :: Auu(ndime,ndime,*)
    real(rp),    optional, intent(inout) :: bu(ndime,*)
    real(rp),    optional, intent(inout) :: uu(ndime,*)
    logical(lg), optional, intent(in)    :: UNIT_DIAGONAL
    integer(ip)                          :: ipoin,idime,izdom,jpoin
    integer(ip)                          :: jzdom,kpoin,izdod,jdime
    integer(ip)                          :: idofn
    real(rp)                             :: Auud(ndime)
    logical(lg)                          :: if_unit_diagonal

    if_unit_diagonal = .false.
    if( present(UNIT_DIAGONAL ) ) if_unit_diagonal = UNIT_DIAGONAL

    do ipoin = 1,npoin

       do idime = 1,ndime

          if( kfl_fixno_nsi(idime,ipoin) > 0 ) then
             !
             ! Eliminate dof of IPOIN from other equations (JPOIN)
             !
             do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                jpoin = c_dom(izdom)
                if( ipoin /= jpoin ) then
                   do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                      kpoin = c_dom(jzdom)
                      if( kpoin == ipoin ) then
                         if( present(Auu) ) then
                            do jdime = 1,ndime
                               if( present(bu) ) bu(jdime,jpoin) = bu(jdime,jpoin) - Auu(idime,jdime,jzdom) * bvess_nsi(idime,ipoin,1)
                               Auu(idime,jdime,jzdom) = 0.0_rp
                            end do
                         end if
                      end if
                   end do

                end if
             end do
             !
             ! IZDOD: Diagonal
             !
             izdod = r_dom(ipoin) - 1
             jpoin = 0
             do while( jpoin /= ipoin )
                izdod = izdod + 1
                jpoin = c_dom(izdod)
             end do
             if( present(Auu) .and. (.not. if_unit_diagonal) ) then                
                Auud(idime) = Auu(idime,idime,izdod)
             else
                Auud(idime) = 1.0_rp
             end if
             if( abs(Auud(idime)) < zeror ) Auud(idime) = 1.0_rp
             !
             ! Set line to zero
             !
             if( present(Auu) ) then
                do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                   do jdime = 1,ndime
                      Auu(jdime,idime,izdom) = 0.0_rp
                   end do
                end do
             end if
             !
             ! Prescribe value
             !
             idofn = (ipoin-1)*ndime + idime
             if( present(Auu) ) then
                Auu(idime,idime,izdod) = Auud(idime)
                if( present(bu) ) bu(idime,ipoin) = bvess_nsi(idime,ipoin,1) * Auud(idime)
                if( present(uu) ) uu(idime,ipoin) = bvess_nsi(idime,ipoin,1)
             else
                !
                ! In the case when Auu is eliminated (bu <= bu - Auu*u), bu
                ! should be zero on boundary condition
                !
                if( present(bu) ) bu(idime,ipoin) = 0.0_rp
                if( present(uu) ) uu(idime,ipoin) = bvess_nsi(idime,ipoin,1)
             end if

          end if

       end do

    end do

  end subroutine nsi_dirichlet_momentum

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-04-17
  !> @brief   Rotate momentum equations
  !> @details Rotate momentum equations
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_dirichlet_rotate(Auu,bu)

    real(rp),    optional, intent(inout) :: Auu(ndime,ndime,*)
    real(rp),    optional, intent(inout) :: bu(ndime,*)

    integer(ip)                          :: ipoin,ibopo,iroty
    real(rp),    pointer                 :: rotma(:,:)
    
    do ipoin = 1,npoin
       ibopo = lpoty(ipoin)
       if( ibopo > 0 ) then
          iroty =  kfl_fixrs_nsi(ibopo)
          if( iroty /= 0 ) then
             if( iroty == -1 ) then                                    ! Tangent system
                rotma => exnor(:,:,ibopo)
             else if( iroty >= 1 ) then                                ! Given system
                rotma =>  skcos(:,:,iroty)
             else if( iroty == -2 ) then                               ! Given system
                rotma => skcos_nsi(:,:,ibopo)
             else if( iroty == -3 ) then                               ! Geometrical normal
                rotma => skcos(:,:,ibopo)
             end if
             call matrix_rotate_system(npoin,ndime,ipoin,r_dom,c_dom,rotma,bu,Auu)
          end if
       end if
    end do

  end subroutine nsi_dirichlet_rotate
  
end module mod_nsi_dirichlet_global_system
!> @}
