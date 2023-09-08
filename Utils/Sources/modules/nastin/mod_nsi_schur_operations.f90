module mod_nsi_schur_operations

  implicit none
  private

  public :: nsi_auuvec
  public :: nsi_aupvec
  public :: nsi_appvec
  public :: nsi_apuvec
  public :: nsi_solini
  public :: nsi_inivec
  public :: nsi_rotsch
  public :: nsi_allvec

contains

  subroutine nsi_aupvec(itask,Aup,pp,uu,xfact_opt,message)
    !
    ! uu = Aup pp
    !
    use def_kintyp
    use def_domain
    use def_master
    use mod_solver, only : solver_smvp
    implicit none
    integer(ip),  intent(in)           :: itask
    real(rp),     intent(in)           :: Aup(ndime,nzdom)
    real(rp),     intent(in)           :: pp(npoin)
    real(rp),     intent(inout)        :: uu(ndime,npoin)
    real(rp),     intent(in), optional :: xfact_opt
    character(*), intent(in), optional :: message
    integer(ip)                        :: izdom,ipoin,jpoin
    real(rp)                           :: xfact
    logical(lg)                        :: use_transpose_Apu

    if( present(xfact_opt) ) then
       xfact = xfact_opt
    else
       xfact = 1.0_rp
    end if
    use_transpose_Apu = .false.
    if( present(message) ) then
       if( trim(message) == 'TRANSPOSE' ) then
          use_transpose_Apu = .true.
       end if
    end if

    if ( INOTMASTER ) then

       if( .not. use_transpose_Apu ) then

          if( itask == 1 ) then
             do ipoin = 1,npoin
                uu(1:ndime,ipoin) = 0.0_rp
             end do
          end if

          if( itask == 2 ) then

             !call solver_matrix_vector_product(npoin,ndime,1_ip,Aup,c_dom,r_dom,pp,uu)

             !$OMP PARALLEL DO SCHEDULE (STATIC)                &
             !$OMP DEFAULT  ( NONE )                            &
             !$OMP PRIVATE  ( ipoin, izdom, jpoin )             &
             !$OMP SHARED   ( r_dom, c_dom, uu, Aup, pp,        &
#ifndef NDIMEPAR
             !$OMP            ndime,                            &
#endif
             !$OMP            npoin, xfact )
             do ipoin = 1,npoin
                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   jpoin = c_dom(izdom)
                   uu(1:ndime,ipoin) = uu(1:ndime,ipoin) - xfact * Aup(1:ndime,izdom) * pp(jpoin)
                end do
             end do
             !$OMP END PARALLEL DO
          else
             !$OMP PARALLEL DO SCHEDULE (STATIC)                &
             !$OMP DEFAULT  ( NONE )                            &
             !$OMP PRIVATE  ( ipoin, izdom, jpoin )             &
             !$OMP SHARED   ( r_dom, c_dom, uu, Aup, pp,        &
#ifndef NDIMEPAR
             !$OMP            ndime,                            &
#endif
             !$OMP            npoin, xfact )
             do ipoin = 1,npoin
                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   jpoin = c_dom(izdom)
                   uu(1:ndime,ipoin) = uu(1:ndime,ipoin) + xfact * Aup(1:ndime,izdom) * pp(jpoin)
                end do
             end do
             !$OMP END PARALLEL DO
          end if

       else

          do ipoin = 1,npoin
             do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                jpoin = c_dom(izdom)
                uu(1:ndime,jpoin) = uu(1:ndime,jpoin) + xfact * Aup(1:ndime,izdom) * pp(ipoin)
             end do
          end do

       end if
    end if

  end subroutine nsi_aupvec

  subroutine nsi_apuvec(itask,Apu,uu,pp,xfact_opt)
    !
    ! ITASK = 0 ... pp = pp + Apu uu
    !         1 ... pp =      Apu uu
    !         2 ... pp = pp - Apu uu
    !
    !
    use def_kintyp
    use def_domain
    use def_master
    implicit none
    integer(ip), intent(in)           :: itask
    real(rp),    intent(in)           :: Apu(ndime,nzdom)
    real(rp),    intent(in)           :: uu(ndime,npoin)
    real(rp),    intent(inout)        :: pp(npoin)
    real(rp),    intent(in), optional :: xfact_opt
    integer(ip)                       :: idime,izdom,ipoin,jpoin
    real(rp)                          :: xfact

    if( present(xfact_opt) ) then
       xfact = xfact_opt
    else
       xfact = 1.0_rp
    end if

    if ( INOTMASTER ) then

       if( itask == 1 ) then
          do ipoin = 1,npoin
             pp(ipoin) = 0.0_rp
          end do
       end if

       if( itask == 2 ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC)               &
          !$OMP DEFAULT  ( NONE )                            &
          !$OMP PRIVATE  ( ipoin, izdom, jpoin, idime )      &
          !$OMP SHARED   ( r_dom, c_dom, uu, Apu, pp,        &
#ifndef NDIMEPAR
          !$OMP            ndime,                            &
#endif
          !$OMP            npoin, xfact )
          do ipoin = 1,npoin
             do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                jpoin     = c_dom(izdom)
                do idime = 1,ndime
                   pp(ipoin) = pp(ipoin) - xfact * Apu(idime,izdom) * uu(idime,jpoin)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO SCHEDULE (STATIC)               &
          !$OMP DEFAULT  ( NONE )                            &
          !$OMP PRIVATE  ( ipoin, izdom, jpoin, idime )      &
          !$OMP SHARED   ( r_dom, c_dom, uu, Apu, pp,        &
#ifndef NDIMEPAR
          !$OMP            ndime,                            &
#endif
          !$OMP            npoin, xfact )
          do ipoin = 1,npoin
             do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                jpoin = c_dom(izdom)
                do idime = 1,ndime
                   pp(ipoin) = pp(ipoin) + xfact * Apu(idime,izdom) * uu(idime,jpoin)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end if

    end if

  end subroutine nsi_apuvec

  subroutine nsi_aputvec(itask,Apu,pp,uu)
    !
    ! ITASK = 0 ... pp = pp + Apu uu
    !         1 ... pp =      Apu uu
    !         2 ... pp = pp - Apu uu
    !
    !
    use def_kintyp
    use def_domain
    use def_master
    implicit none
    integer(ip), intent(in)    :: itask
    real(rp),    intent(in)    :: Apu(ndime,nzdom),pp(npoin)
    real(rp),    intent(inout) :: uu(ndime,npoin)
    integer(ip)                :: idime,izdom,ipoin,jpoin

    if ( INOTMASTER ) then

       if( itask == 1 ) then
          uu(1:ndime,1:npoin) = 0.0_rp
       end if

       if( itask == 2 ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC)               &
          !$OMP DEFAULT  ( NONE )                            &
          !$OMP PRIVATE  ( ipoin, izdom, jpoin, idime )      &
          !$OMP SHARED   ( r_dom, c_dom, uu, Apu, pp,        &
#ifndef NDIMEPAR
          !$OMP            ndime,                            &
#endif
          !$OMP            npoin )
          do ipoin = 1,npoin
             do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                jpoin = c_dom(izdom)
                do idime = 1,ndime
                   uu(idime,jpoin) = uu(idime,jpoin) - Apu(idime,izdom) * pp(ipoin)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO SCHEDULE (STATIC)               &
          !$OMP DEFAULT  ( NONE )                            &
          !$OMP PRIVATE  ( ipoin, izdom, jpoin, idime )      &
          !$OMP SHARED   ( r_dom, c_dom, uu, Apu, pp,        &
#ifndef NDIMEPAR
          !$OMP            ndime,                            &
#endif
          !$OMP            npoin )
          do ipoin = 1,npoin
             do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                jpoin = c_dom(izdom)
                do idime = 1,ndime
                   uu(idime,jpoin) = uu(idime,jpoin) + Apu(idime,izdom) * pp(ipoin)
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end if

    end if

  end subroutine nsi_aputvec

  subroutine nsi_appvec(itask,App,pp,qq,xfact_opt)
    !
    ! qq = App pp
    !
    use def_kintyp
    use def_domain
    use def_nastin
    use def_master
    implicit none
    integer(ip), intent(in)              :: itask
    real(rp),    intent(in)              :: App(*)
    real(rp),    intent(in)              :: pp(npoin)
    real(rp),    intent(inout), target   :: qq(npoin)
    real(rp),    intent(in),    optional :: xfact_opt
    integer(ip)                          :: kk,ipoin,jpoin,izsym,izdom
    real(rp)                             :: raux,raux2,xfact


    if( present(xfact_opt) ) then
       xfact = xfact_opt
    else
       xfact = 1.0_rp
    end if


    if ( INOTMASTER ) then

       if( itask == 1 ) then
          do ipoin = 1,npoin
             qq(ipoin) = 0.0_rp
          end do
       end if

       if( itask == 2 ) then
          if( solve(2) % kfl_symme == 1 ) then
             do ipoin = 1, npoin
                raux  = 0.0_rp
                raux2 = pp(ipoin)
                do izsym= r_sym(ipoin), r_sym(ipoin+1)-2
                   jpoin     = c_sym(izsym)
                   raux      = raux + App(izsym) * pp(jpoin)
                   qq(jpoin) = qq(jpoin) - App(izsym) * raux2
                end do
                kk        = r_sym(ipoin+1)-1
                qq(ipoin) = qq(ipoin) - raux - App(kk) * pp(c_sym(kk))
             end do
          else
             !$OMP PARALLEL DO SCHEDULE (STATIC)              &
             !$OMP DEFAULT  ( NONE )                           &
             !$OMP PRIVATE  ( ipoin, izdom, jpoin )            &
             !$OMP SHARED   ( r_dom, c_dom, App, pp, qq,       &
             !$OMP            npoin, xfact  )
             do ipoin = 1,npoin
                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   qq(ipoin) = qq(ipoin) - xfact * App(izdom) * pp(c_dom(izdom))
                end do
             end do
             !$OMP END PARALLEL DO
          end if
       else
          if( solve(2) % kfl_symme == 1 ) then
             do ipoin = 1,npoin
                raux  = 0.0_rp
                raux2 = pp(ipoin)
                do izsym= r_sym(ipoin), r_sym(ipoin+1)-2
                   jpoin     = c_sym(izsym)
                   raux      = raux + App(izsym) * pp(jpoin)
                   qq(jpoin) = qq(jpoin) + App(izsym) * raux2
                end do
                kk        = r_sym(ipoin+1)-1
                qq(ipoin) = qq(ipoin) + raux + App(kk) * pp(c_sym(kk))
             end do
          else
             !$OMP PARALLEL DO SCHEDULE (STATIC)              &
             !$OMP DEFAULT  ( NONE )                           &
             !$OMP PRIVATE  ( ipoin, izdom, jpoin )            &
             !$OMP SHARED   ( r_dom, c_dom, App, pp, qq ,      &
             !$OMP            npoin, xfact )
             do ipoin = 1,npoin
                do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                   qq(ipoin) = qq(ipoin) + xfact * App(izdom) * pp(c_dom(izdom))
                end do
             end do
             !$OMP END PARALLEL DO
          end if
       end if

    end if

  end subroutine nsi_appvec

  subroutine nsi_auuvec(itask,Auu,vv,uu)
    !
    ! uu = Auu vv
    !
    use def_kintyp
    use def_domain
    use def_master
    implicit none
    integer(ip), intent(in)    :: itask
    real(rp),    intent(in)    :: Auu(ndime,ndime,nzdom)
    real(rp),    intent(in)    :: vv(ndime,npoin)
    real(rp),    intent(inout) :: uu(ndime,npoin)
    integer(ip)                :: idime,izdom,ipoin,jpoin,jdime

    if ( INOTMASTER ) then

       if( itask == 1 ) then
          do ipoin = 1,npoin
             uu(1:ndime,ipoin) = 0.0_rp
          end do
       end if

       if( itask == 2 ) then
          !$OMP PARALLEL DO SCHEDULE (STATIC)                 &
          !$OMP DEFAULT  ( NONE )                              &
          !$OMP PRIVATE  ( ipoin, izdom, jpoin, idime, jdime ) &
          !$OMP SHARED   ( r_dom, c_dom, Auu, uu, vv,          &
#ifndef NDIMEPAR
          !$OMP            ndime,                              &
#endif
          !$OMP            npoin )
          do ipoin = 1,npoin
             do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                jpoin = c_dom(izdom)
                do idime = 1,ndime
                   do jdime = 1,ndime
                      uu(idime,ipoin) = uu(idime,ipoin) - Auu(jdime,idime,izdom) * vv(jdime,jpoin)
                   end do
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       else
          !$OMP PARALLEL DO SCHEDULE (STATIC)                 &
          !$OMP DEFAULT  ( NONE )                              &
          !$OMP PRIVATE  ( ipoin, izdom, jpoin, idime, jdime ) &
          !$OMP SHARED   ( r_dom, c_dom, Auu, uu, vv,          &
#ifndef NDIMEPAR
          !$OMP            ndime,                              &
#endif
          !$OMP            npoin )
          do ipoin = 1,npoin
             do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                jpoin = c_dom(izdom)
                do idime = 1,ndime
                   do jdime = 1,ndime
                      uu(idime,ipoin) = uu(idime,ipoin) + Auu(jdime,idime,izdom) * vv(jdime,jpoin)
                   end do
                end do
             end do
          end do
          !$OMP END PARALLEL DO
       end if

    end if

  end subroutine nsi_auuvec

  subroutine nsi_solini(itask)
    !-----------------------------------------------------------------------
    !****f* Nastin/nsi_solini
    ! NAME
    !    nsi_solini
    ! DESCRIPTION
    !    This routine loads the solver data for the incomp. NS equations.
    !    In general, it may change from time step to time step or even
    !    from iteration to iteration.
    ! USED BY
    !    nsi_begite
    !***
    !-----------------------------------------------------------------------
    use def_parame
    use def_master
    use def_domain
    use def_nastin
    use def_solver
    implicit none
    integer(ip), intent(in) :: itask

    if( itask == 1 ) then
       !
       ! Momentum
       !
       ivari_nsi            =  ivari_nsi_mom
       !solve(1) % ndofn     =  ndime
       !solve(1) % ndof2     =  solve(1) % ndofn**2
       !solve(1) % nzmat     =  solve(1) % ndof2 * nzdom
       !solve(1) % nunkn     =  solve(1) % ndofn * npoin
       !solve(1) % nzrhs     =  solve(1) % ndofn * npoin
       solve(1) % bvess     => bvess_nsi(:,:,1)
       solve_sol            => solve(1:)

       if( solve(1) % kfl_algso == SOL_SOLVER_RICHARDSON .or. solve(1) % kfl_algso == SOL_SOLVER_MATRIX_RICHARDSON ) then
          solve(1) % xdiag = 1.0_rp / dtinv_nsi
       end if

    else if( itask == 2 ) then
       !
       ! Continuity
       !
       ivari_nsi          =  ivari_nsi_cont
       solve_sol          => solve(2:)

    else if( itask == 3 ) then
       !
       ! Momentum + continuity (used for example to initialize matrix)
       !
       return
       
       ivari_nsi        =  ivari_nsi_mom
       solve(1) % ndofn =  ndime + 1
       solve(1) % ndof2 =  solve(1) % ndofn**2
       solve(1) % nzmat =  nmauu_nsi + nmaup_nsi + nmapu_nsi + nmapp_nsi
       solve(1) % nunkn =  solve(1) % ndofn * npoin
       solve(1) % nzrhs =  solve(1) % ndofn * npoin
       solve_sol        => solve(1:)

       if( IMASTER ) then
          Auu_nsi => nul1r
          Aup_nsi => nul1r
          Apu_nsi => nul1r
          App_nsi => nul1r
       else
          Auu_nsi => amatr(poauu_nsi:poauu_nsi+nmauu_nsi-1)
          Aup_nsi => amatr(poaup_nsi:poaup_nsi+nmaup_nsi-1)
          Apu_nsi => amatr(poapu_nsi:poapu_nsi+nmapu_nsi-1)
          App_nsi => amatr(poapp_nsi:poapp_nsi+nmapp_nsi-1)
       end if

    else if( itask == 4 ) then
       !
       ! Momentum variation (second solve in Du)
       !
       ivari_nsi            =  ivari_nsi_mom
       !solve(1) % ndofn     =  ndime
       !solve(1) % ndof2     =  solve(1) % ndofn**2
       !solve(1) % nzmat     =  solve(1) % ndof2 * nzdom
       !solve(1) % nunkn     =  solve(1) % ndofn * npoin
       !solve(1) % nzrhs     =  solve(1) % ndofn * npoin
       solve(1) % bvess     => null()
       solve_sol            => solve(1:)

    end if

  end subroutine nsi_solini

  subroutine nsi_inivec(ncomp,xx)
    !
    ! uu = Aup pp
    !
    use def_kintyp
    use def_domain
    use def_master
    implicit none
    integer(ip), intent(in)  :: ncomp
    real(rp),    intent(out) :: xx(ncomp)
    integer(ip)              :: icomp

    if( INOTMASTER ) then
       do icomp = 1,ncomp
          xx(icomp) = 0.0_rp
       end do
    end if

  end subroutine nsi_inivec

  subroutine nsi_rotsch(itask,uu,given_exnor)
    !-----------------------------------------------------------------------
    !****f* Nastin/nsi_rotunk
    ! NAME
    !    nsi_rotunk
    ! DESCRIPTION
    !    This routine rotates the nodal velocities using the appropiate
    !    rotation matrix:
    !    ITASK=1 ... From global to local
    !    ITASK=2 ... From local to global
    !    Modifications need to be done only if there exist image nodes or
    !    boundary conditions in skew systems.
    ! USES
    !    mbvab0
    ! USED BY
    !    nsi_solite
    !***
    !-----------------------------------------------------------------------
    use def_domain
    use def_master
    use def_nastin
    use mod_maths,  only     :  maths_local_orthonormal_basis
    implicit none
    integer(ip), intent(in)              :: itask
    real(rp)   , intent(in), optional    :: given_exnor(3)
    real(rp),    intent(inout) :: uu(*)
    integer(ip)                :: ipoin,ibopo,iroty,itotv,jerro
    real(rp)                   :: worma(3),worve(3),given_rotma(3,3)

    if( kfl_local_nsi /= 0 .and. INOTMASTER ) then

       if( itask == 1 ) then
          !
          ! Global to local
          !
          do ipoin = 1,npoin
             ibopo = lpoty(ipoin)
             if( ibopo > 0 ) then
                iroty = kfl_fixrs_nsi(ibopo)

                if( iroty == -1 ) then
                   !
                   ! Boundary conditions in the tangent skew system
                   !
                   itotv=(ipoin-1)*ndime
                   worve(1:ndime)=uu(itotv+1:itotv+ndime)
                   if (present(given_exnor)) then
                      jerro = 0
                      given_rotma= 0.0_rp
                      given_rotma(1:3,1)=given_exnor(1:3)
                      call maths_local_orthonormal_basis(ndime,given_rotma,jerro)
                      if( jerro /= 0 ) then
                         write(*,*) 'MOD_NSI_SCHUR_OPERATIONS: TANGENT COULD NOT BE COMPUTED FOR NODE ',lninv_loc(ipoin)
                      end if                      
                      call mbvatb(worma,given_rotma(1,1),worve,ndime,ndime)
                   else
                      call mbvatb(worma,exnor(1,1,ibopo),worve,ndime,ndime)
                   end if

                   uu(itotv+1:itotv+ndime)=worma(1:ndime)

                else if( iroty == -2 ) then
                   !
                   ! Boundary conditions in the NSI tangent skew system
                   !
                   itotv=(ipoin-1)*ndime
                   worve(1:ndime)=uu(itotv+1:itotv+ndime)
                   call mbvatb(worma,skcos_nsi(1,1,ibopo),worve,ndime,ndime)
                   uu(itotv+1:itotv+ndime)=worma(1:ndime)

                else if( iroty == -3 ) then
                   !
                   ! Boundary conditions in geometrical system
                   !
                   itotv=(ipoin-1)*ndime
                   worve(1:ndime)=uu(itotv+1:itotv+ndime)
                   call mbvatb(worma,skcos(1,1,ibopo),worve,ndime,ndime)
                   uu(itotv+1:itotv+ndime)=worma(1:ndime)

                else if( iroty >= 1 ) then
                   !
                   ! Boundary conditions in a given skew system
                   !
                   itotv=(ipoin-1)*ndime
                   worve(1:ndime)=uu(itotv+1:itotv+ndime)
                   call mbvatb(worma,skcos(1,1,iroty),worve,ndime,ndime)
                   uu(itotv+1:itotv+ndime)=worma(1:ndime)

                end if
             end if
          end do

       else if( itask == 2 ) then
          !
          ! Local to global
          !
          do ipoin = 1,npoin
             ibopo = lpoty(ipoin)
             if( ibopo > 0 ) then
                iroty = kfl_fixrs_nsi(ibopo)
                if( iroty == -1 ) then
                   !
                   ! Boundary conditions in the tangent skew system
                   !
                   itotv=(ipoin-1)*ndime
                   worve(1:ndime)=uu(itotv+1:itotv+ndime)

                   if (present(given_exnor)) then
                      jerro = 0
                      given_rotma= 0.0_rp
                      given_rotma(1:3,1)=given_exnor(1:3)
                      call maths_local_orthonormal_basis(ndime,given_rotma,jerro)
                      if( jerro /= 0 ) then
                         write(*,*) 'MOD_NSI_SCHUR_OPERATIONS: TANGENT COULD NOT BE COMPUTED FOR NODE ',lninv_loc(ipoin)
                      end if                      
                      call mbvab0(worma,given_rotma(1,1),worve,ndime,ndime)
                   else
                      call mbvab0(worma,exnor(1,1,ibopo),worve,ndime,ndime)
                   end if

                   uu(itotv+1:itotv+ndime)=worma(1:ndime)

                else if( iroty == -2 ) then
                   !
                   ! Boundary conditions in the NSI skew system
                   !
                   itotv=(ipoin-1)*ndime
                   worve(1:ndime)=uu(itotv+1:itotv+ndime)
                   call mbvab0(worma,skcos_nsi(1,1,ibopo),worve,ndime,ndime)
                   uu(itotv+1:itotv+ndime)=worma(1:ndime)

                else if( iroty == -3 ) then
                   !
                   ! Boundary conditions in geometrical system
                   !
                   itotv=(ipoin-1)*ndime
                   worve(1:ndime)=uu(itotv+1:itotv+ndime)
                   call mbvab0(worma,skcos(1,1,ibopo),worve,ndime,ndime)
                   uu(itotv+1:itotv+ndime)=worma(1:ndime)

                else if( iroty >= 1 ) then
                   !
                   ! Boundary conditions in a given skew system
                   !
                   itotv=(ipoin-1)*ndime
                   worve(1:ndime)=uu(itotv+1:itotv+ndime)
                   call mbvab0(worma,skcos(1,1,iroty),worve,ndime,ndime)
                   uu(itotv+1:itotv+ndime)=worma(1:ndime)
                end if
             end if
          end do

       end if

    end if

  end subroutine nsi_rotsch

  subroutine nsi_allvec(itask,Auu,Aup,App,Apu,vv,uu)
    !
    ! uu = Auu vv
    !
    use def_kintyp
    use def_domain
    use def_master
    implicit none
    integer(ip), intent(in)  :: itask
    real(rp),    intent(in)  :: Auu(ndime,ndime,nzdom)
    real(rp),    intent(in)  :: Aup(ndime,nzdom)
    real(rp),    intent(in)  :: App(nzdom)
    real(rp),    intent(in)  :: Apu(ndime,nzdom)
    real(rp),    intent(in)  :: vv(ndime+1,npoin)
    real(rp),    intent(out) :: uu(ndime+1,npoin)
    integer(ip)              :: idime,izdom,ipoin,jpoin,jdime,ndofn

    if ( INOTMASTER ) then

       ndofn = ndime + 1

       do ipoin = 1,npoin
          uu(1:ndime,ipoin) = 0.0_rp
       end do

       !$OMP PARALLEL DO SCHEDULE (STATIC)                 &
       !$OMP DEFAULT  ( NONE )                              &
       !$OMP PRIVATE  ( ipoin, izdom, jpoin, idime, jdime ) &
       !$OMP SHARED   ( r_dom, c_dom, Auu, Aup, Apu, App,   &
#ifndef NDIMEPAR
       !$OMP         ndime,                              &
#endif
       !$OMP            uu, vv, ndofn, npoin )
       do ipoin = 1,npoin
          do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
             jpoin = c_dom(izdom)
             do idime = 1,ndime
                do jdime = 1,ndime
                   uu(idime,ipoin) = uu(idime,ipoin) + Auu(jdime,idime,izdom) * vv(jdime,jpoin)
                end do
                uu(idime,ipoin) = uu(idime,ipoin) + Aup(idime,izdom) * vv(ndofn,jpoin)
                uu(ndofn,ipoin) = uu(ndofn,ipoin) + Apu(idime,izdom) * vv(idime,jpoin)
             end do
             uu(ndofn,ipoin) = uu(ndofn,ipoin) + App(izdom) * vv(ndofn,jpoin)
          end do
       end do
       !$OMP END PARALLEL DO

       if( itask == 1 ) call pararr('SLX',NPOIN_TYPE,npoin*ndofn,uu)

    end if

  end subroutine nsi_allvec

end module mod_nsi_schur_operations



