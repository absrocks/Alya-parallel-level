!-----------------------------------------------------------------------
!> @addtogroup Algebraic_Solver
!> @{
!> @file    presol.f90
!> @date    22/10/2013
!> @author  Guillaume Houzeaux
!> @brief   Operations on matrix and solvers
!> @details Perform some modifications on the matrix and
!>          or RHS to account for:
!>          \verbatim
!>
!>          ---------------------------------------------------
!>          ITASK     Periodicity             Dirichlet
!>          ---------------------------------------------------
!>                    amatr                   amatr
!>                      +        rhsid          +        rhsid
!>                    rhsid                   rhsid
!>          ---------------------------------------------------
!>            1        X           X            X          X
!>            2        -           X            -          -
!>            3        X           X            -          -
!>            4        -           X            -          X
!>          ---------------------------------------------------
!>
!>         \endverbatim
!> @}
!-----------------------------------------------------------------------
subroutine presol(itask,ndof1,ndof2,rhsid,amatr,unkno)
  use def_kintyp,         only : ip,rp,lg
  use def_master,         only : kfl_paral,IMASTER,gisca,lninv_loc
  use def_master,         only : zeror,current_zone
  use def_solver,         only : solve_sol,memit
  use def_domain,         only : r_dom,c_dom,nperi,lperi,npoin,&
       &                         nzone
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_memory,         only : memory_alloca,memory_deallo

  use def_domain, only : kfl_codno
  use def_master, only : ID_TEMPER,modul,current_code
  implicit none
  integer(ip), intent(in)    :: itask
  integer(ip), intent(in)    :: ndof1
  integer(ip), intent(in)    :: ndof2
  real(rp),    intent(inout) :: rhsid(ndof1,*)
  real(rp),    intent(inout) :: amatr(ndof2,ndof1,*)
  real(rp),    intent(inout) :: unkno(ndof1,*)
  integer(ip)                :: izdom,jzdom,ifoun
  integer(ip)                :: idofn,ipoin,jpoin,kpoin,jj,ii
  integer(ip)                :: iperi,kperi,kzdom,jdofn,izdod,ndofn,jtotn
  integer(ip)                :: lpoin,ndofn_i,ndofn_j,iblok,jblok,num_blocks
  real(rp)                   :: amatrd(ndof1)
  integer(ip), pointer       :: kfl_fixno(:,:)
  real(rp),    pointer       :: bvess(:,:)
  logical(lg), pointer       :: in_my_zone(:)
  real(rp),    pointer       :: reaction(:,:)

  nullify( kfl_fixno )
  nullify( bvess )
  nullify( in_my_zone )
  nullify( reaction )

  if( IMASTER ) return

  !----------------------------------------------------------------------
  !
  ! Periodicity
  !
  !----------------------------------------------------------------------

  if( nperi /= 0 ) then

     if( itask == 0 .or. itask == 4 ) then
        !
        ! Only RHS
        !
        if( nzone > 1 ) then
           call memory_alloca(memit,'IN_MY_ZONE','presol',in_my_zone,npoin)
           do ipoin = 1,npoin
              in_my_zone(ipoin) = .true.
           end do
           do iperi = 1,nperi
              jpoin = lperi(2,iperi)
              ipoin = lperi(1,iperi)
              if ( ipoin > 0 .and. jpoin > 0 ) then
                 if( in_my_zone(ipoin) .and. in_my_zone(jpoin) ) then
                    rhsid(1:ndof1,ipoin) = rhsid(1:ndof1,ipoin) + rhsid(1:ndof1,jpoin)
                    rhsid(1:ndof1,jpoin) = rhsid(1:ndof1,ipoin)
                 end if
              end if
           end do
           call memory_deallo(memit,'IN_MY_ZONE','presol',in_my_zone)
        else
           do iperi = 1,nperi
              jpoin = lperi(2,iperi)
              ipoin = lperi(1,iperi)
              if( ipoin > 0 .and. jpoin > 0 ) then
                 rhsid(1:ndof1,ipoin) = rhsid(1:ndof1,ipoin) + rhsid(1:ndof1,jpoin)
                 rhsid(1:ndof1,jpoin) = rhsid(1:ndof1,ipoin)
              end if
           end do
        end if

     else if( itask == 1 .or. itask == 3 ) then
        !
        ! Check errors
        !
        if( solve_sol(1) % kfl_algso == 0 ) then
           call runend('PERIODICITY ONLY POSSIBLE WITH SPARSE SOLVERS')
        else if( solve_sol(1) % kfl_symme == 1 ) then
           call runend('MATRIX PERIODICITY ONLY POSSIBLE WITH ASYMMETRIC SOLVERS')
        end if
        !
        ! Put all slaves' columns into master's column
        ! KPOIN is neither master nor slave
        ! JPOIN is slave
        ! IPOIN is JPOIN's master
        !
        call memory_alloca(memit,'IN_MY_ZONE','presol',in_my_zone,npoin)
        do ipoin = 1,npoin
           in_my_zone(ipoin) = .true.
        end do

        do iperi = 1,nperi
           jpoin = lperi(2,iperi)
           ipoin = lperi(1,iperi)
           if ( ipoin <= 0 .or. jpoin <= 0 ) then
              in_my_zone(abs(ipoin)) = .false.
           end if
        end do

        call memgen(1_ip,npoin,0_ip)
        do iperi = 1,nperi
           kperi = size(lperi,1,KIND=ip)
           ipoin = lperi(1,iperi)
           do ii = 2,kperi
              jpoin = lperi(ii,iperi)
              if ( ipoin>0 .and. jpoin>0 ) then
                 if( in_my_zone(ipoin) .and. in_my_zone(jpoin) ) gisca(jpoin) = ipoin
                 gisca(jpoin) = ipoin
              end if
           end do
        end do

        do kpoin = 1,npoin

           if( gisca(kpoin) == 0 ) then
              do kzdom = r_dom(kpoin),r_dom(kpoin+1)-1
                 jpoin = c_dom(kzdom)

                 if( gisca(jpoin) > 0 .and. in_my_zone(jpoin) ) then
                    ipoin = gisca(jpoin)

                    ifoun = 0
                    izdom2: do izdom = r_dom(kpoin),r_dom(kpoin+1)-1
                       if( c_dom(izdom) == ipoin ) then
                          ifoun = 1

                          do ii = 1,ndof2
                             do jj = 1,ndof1
                                amatr(ii,jj,izdom) = amatr(ii,jj,izdom) + amatr(ii,jj,kzdom)
                                amatr(ii,jj,kzdom) = 0.0_rp
                             end do
                          end do
                          exit izdom2
                       end if
                    end do izdom2
                    if( ifoun == 0 ) then
                       print*,'NOT FOUND NODE, SLAVE, MASTER A=',kpoin,jpoin,ipoin
                       call runend('PRESOL 2: NODE NOT FOUND')
                    end if
                 end if
              end do
           end if
        end do
        !
        ! Put all slaves' rows into master's row
        ! Put slave row to zero
        ! JPOIN:  slave
        ! IPOIN:  master
        ! Slave:  coef. JZDOM for JPOIN-KPOIN
        ! Master: coef. IZOMD for IPOIN-KPOIN
        !
        do jpoin = 1,npoin
           if( gisca(jpoin) > 0 ) then
              ipoin = gisca(jpoin)
              if( itask == 1 ) then
                 rhsid(1:ndof1,ipoin) = rhsid(1:ndof1,ipoin) + rhsid(1:ndof1,jpoin)
                 rhsid(1:ndof1,jpoin) = 0.0_rp
              end if
              do jzdom = r_dom(jpoin),r_dom(jpoin+1)-1
                 kpoin = c_dom(jzdom)
                 if( in_my_zone(kpoin) ) then
                    if( gisca(kpoin) > 0 ) kpoin = gisca(kpoin)
                    ifoun = 0
                    izdom1: do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                       if( c_dom(izdom) == kpoin ) then
                          ifoun = 1
                          do ii = 1,ndof2
                             do jj = 1,ndof1
                                amatr(ii,jj,izdom) = amatr(ii,jj,izdom) + amatr(ii,jj,jzdom)
                                amatr(ii,jj,jzdom) = 0.0_rp
                             end do
                          end do
                          exit izdom1
                       end if
                    end do izdom1
                    if( ifoun == 0 .and. ipoin /= kpoin ) then
                       print*,'NOT FOUND NODE, SLAVE, MASTER B=',kfl_paral,lninv_loc(kpoin),lninv_loc(jpoin),lninv_loc(ipoin)
                       call runend('PRESOL 1: NODE NOT FOUND')
                    end if
                 end if
              end do

           end if
        end do

        do jpoin = 1,npoin
           if( gisca(jpoin) > 0 ) then
              izdom3: do jzdom = r_dom(jpoin),r_dom(jpoin+1)-1
                 if( c_dom(jzdom) == jpoin ) then
                    do ii = 1,ndof1
                       ! amatr(ii,ii,jzdom) = 1.0_rp
                    end do
                    exit izdom3
                 end if
              end do izdom3
           end if
        end do

        call memgen(3_ip,npoin,0_ip)
        call memory_deallo(memit,'IN_MY_ZONE','presol',in_my_zone)

     else if( itask == 2 ) then
        !
        ! Solution
        !
        if( nzone > 1 ) then
           call memory_alloca(memit,'IN_MY_ZONE','presol',in_my_zone,npoin)
           do ipoin = 1,npoin
              in_my_zone(ipoin) = .false.
           end do
           do ipoin = 1,npoin
              in_my_zone(ipoin) = .true.
           end do
           do iperi = 1,nperi
              jpoin = lperi(2,iperi)
              ipoin = lperi(1,iperi)
              if( ipoin > 0 .and. jpoin > 0 ) then
                 if( in_my_zone(ipoin) .and. in_my_zone(jpoin) ) then
                    unkno(1:ndof1,jpoin) = unkno(1:ndof1,ipoin)
                 end if
              end if
           end do
           call memory_deallo(memit,'IN_MY_ZONE','presol',in_my_zone)
        else
           do iperi = 1,nperi
              jpoin = lperi(2,iperi)
              ipoin = lperi(1,iperi)
              if( ipoin > 0 .and. jpoin > 0 ) then
                 unkno(1:ndof1,jpoin) = unkno(1:ndof1,ipoin)
              end if
           end do
        end if

     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Reaction on Dirichlet nodes: save matrix and RHS before imposing Dirichlet b.c.s
  !
  !----------------------------------------------------------------------

  if( itask == 1 .or. itask == 4 ) then

     if( associated(solve_sol(1) % lpoin_reaction) .and. solve_sol(1) % block_num == 1 ) then
        !
        ! Allocate
        !
        num_blocks = solve_sol(1) % num_blocks

        if( num_blocks == 1 ) then
           !
           ! Monolithic system
           !
           ndofn = solve_sol(1) % ndofn
           do ipoin = 1,npoin
              if( solve_sol(1) % lpoin_reaction(ipoin) ) then
                 solve_sol(1) % lpoin_block(ipoin) % block1_num(1) % rhs(1:ndofn) = rhsid(1:ndofn,ipoin)
                 jzdom = 0
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jzdom = jzdom + 1
                    solve_sol(1) % lpoin_block(ipoin) % block2_num(1,1) % matrix(1:ndofn,1:ndofn,jzdom) = amatr(1:ndofn,1:ndofn,izdom)
                 end do
              end if
           end do

        else

           call runend('PRESOL: NOT CODED FOR THIS NUMBER OF BLOCKS')

        end if
     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Impose Neumann b.c.
  !
  !----------------------------------------------------------------------

  if( itask == 1 ) then
     if( solve_sol(1) % kfl_bvnat == 1 ) then
        if( solve_sol(1) % kfl_iffix == 1 ) then
           kfl_fixno => solve_sol(1) % kfl_fixno
           do ipoin = 1,npoin
              do idofn = 1,ndof1
                 if( solve_sol(1) % kfl_fixno(idofn,ipoin) <= 0 ) then
                    rhsid(idofn,ipoin) = rhsid(idofn,ipoin) + solve_sol(1) % bvnat(idofn,ipoin)
                 end if
              end do
           end do
        else
           do ipoin = 1,npoin
              rhsid(1:ndof1,ipoin) = rhsid(1:ndof1,ipoin) + solve_sol(1) % bvnat(1:ndof1,ipoin)
           end do
        end if
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! Impose Dirichlet b.c.
  !
  !----------------------------------------------------------------------

  if( itask == 4 ) then
     !
     ! Richardson residual-based solvers: put RHS=0 on Dirichlet nodes
     !
     kfl_fixno => solve_sol(1) % kfl_fixno

     if( solve_sol(1) % kfl_iffix == 1 ) then
        do ipoin = 1,npoin
           do idofn = 1,ndof1
              if( kfl_fixno(idofn,ipoin) > 0 ) then
                 rhsid(idofn,ipoin) = 0.0_rp
              end if
           end do
        end do
     end if

  else if( itask == 1 ) then

     kfl_fixno => solve_sol(1) % kfl_fixno
     bvess     => solve_sol(1) % bvess

     if( solve_sol(1) % kfl_iffix == 1 .and. associated(solve_sol(1) % bvess) ) then
        !
        ! Dirichlet value is not given in BVESS
        !
        do ipoin = 1,npoin

           do idofn = 1,ndof1

              if( kfl_fixno(idofn,ipoin) > 0 ) then
                 !
                 ! Eliminate dof of IPOIN from other equations (JPOIN)
                 ! Keep rows unchanged in order to compute the reaction force
                 !
                 unkno(idofn,ipoin) = bvess(idofn,ipoin)
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    jpoin = c_dom(izdom)
                    if( ipoin /= jpoin ) then
                       do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                          kpoin = c_dom(jzdom)
                          if( kpoin == ipoin ) then
                             do jdofn = 1,ndof1
                                rhsid(jdofn,jpoin)       = rhsid(jdofn,jpoin) - amatr(idofn,jdofn,jzdom) * bvess(idofn,ipoin)
                                amatr(idofn,jdofn,jzdom) = 0.0_rp
                             end do
                          end if
                       end do
                    end if
                 end do
                 !
                 ! => Uncomment to cancel out rows with Dirichlet b.c.
                 !
                 !
                 ! IZDOD: Diagonal
                 !
                 izdod = r_dom(ipoin) - 1
                 jpoin = 0
                 do while( jpoin /= ipoin )
                    izdod = izdod + 1
                    jpoin = c_dom(izdod)
                 end do
                 amatrd(idofn) = amatr(idofn,idofn,izdod)
                 if( abs(amatrd(idofn)) < zeror ) amatrd(idofn) = 1.0_rp
                 !
                 ! Set line to zero
                 !
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    do jdofn = 1,ndof1
                       amatr(jdofn,idofn,izdom) = 0.0_rp
                    end do
                 end do
                 !
                 ! Presrcibe value
                 !
                 amatr(idofn,idofn,izdod) = amatrd(idofn)
                 rhsid(idofn,ipoin)       = bvess(idofn,ipoin) * amatrd(idofn)
                 !
                 ! <= Uncomment
                 !
              end if

           end do

        end do

     else if( solve_sol(1) % kfl_iffix == 2 .or. ( solve_sol(1) % kfl_iffix == 1 .and. .not. associated(solve_sol(1) % bvess) ) ) then
        !
        ! Dirichlet value is not given: Impose zero
        !
        do ipoin = 1,npoin

           do idofn = 1,ndof1

              if( kfl_fixno(idofn,ipoin) > 0 ) then
                 !
                 ! Eliminate dof of IPOIN from other equations (JPOIN)
                 !
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    jpoin = c_dom(izdom)
                    if( ipoin /= jpoin ) then
                       do jzdom = r_dom(jpoin),r_dom(jpoin+1) - 1
                          kpoin = c_dom(jzdom)
                          if( kpoin == ipoin ) then
                             do jdofn = 1,ndof1
                                amatr(idofn,jdofn,jzdom) = 0.0_rp
                             end do
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
                 amatrd(idofn) = amatr(idofn,idofn,izdod)
                 if( abs(amatrd(idofn)) < zeror ) amatrd(idofn) = 1.0_rp
                 !
                 ! Set line to zero
                 !
                 do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                    do jdofn = 1,ndof1
                       amatr(jdofn,idofn,izdom) = 0.0_rp
                    end do
                 end do
                 !
                 ! Prescribe value
                 !
                 amatr(idofn,idofn,izdod) = amatrd(idofn)
                 if( solve_sol(1) % kfl_iffix /= 2 ) &
                      rhsid(idofn,ipoin)       = 0.0_rp

              end if

           end do

        end do

     end if

  end if

  !----------------------------------------------------------------------
  !
  ! Recover Dirichlet condition
  !
  !----------------------------------------------------------------------

!!$  if( itask == 5 ) then
!!$
!!$     if( solve_sol(1) % kfl_iffix == 1 ) then
!!$        if( associated(solve_sol(1) % bvess) ) then
!!$           do lpoin = 1,npoiz(current_zone)
!!$              ipoin = lpoiz(current_zone) % l(lpoin)
!!$              do idofn = 1,ndof1
!!$                 if( solve_sol(1) % kfl_fixno(idofn,ipoin) > 0 ) unkno(idofn,ipoin) = solve_sol(1) % bvess(idofn,ipoin)
!!$              end do
!!$           end do
!!$        else
!!$           do lpoin = 1,npoiz(current_zone)
!!$              ipoin = lpoiz(current_zone) % l(lpoin)
!!$              do idofn = 1,ndof1
!!$                 if( solve_sol(1) % kfl_fixno(idofn,ipoin) > 0 ) unkno(idofn,ipoin) = 0.0_rp
!!$              end do
!!$           end do
!!$        end if
!!$     end if
!!$
!!$  end if

  !----------------------------------------------------------------------
  !
  ! Reaction on Dirichlet nodes: save matrix and RHS before imposing Dirichlet b.c.s
  !
  !----------------------------------------------------------------------

  if( itask == 5 ) then

     if( solve_sol(1) % kfl_react == 1 .and. associated(solve_sol(1) % lpoin_reaction) ) then

        num_blocks = solve_sol(1) % num_blocks

        if( num_blocks == 1 ) then
           !
           ! Monolithic system
           !
           ndofn = solve_sol(1) % ndofn
           do ipoin = 1,npoin
              solve_sol(1) % reaction(1:ndofn,ipoin) = 0.0_rp
              if( solve_sol(1) % lpoin_reaction(ipoin) ) then
                 jzdom = 0
                 solve_sol(1) % reaction(1:ndofn,ipoin) = solve_sol(1) % lpoin_block(ipoin) % block1_num(1) % rhs(1:ndofn)
                 do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                    jpoin = c_dom(izdom)
                    jzdom = jzdom + 1
                    do idofn = 1,ndofn
                       do jdofn = 1,ndofn
                          solve_sol(1) % reaction(idofn,ipoin) = solve_sol(1) % reaction(idofn,ipoin) &
                               & - solve_sol(1) % lpoin_block(ipoin) % block2_num(1,1) % matrix(jdofn,idofn,jzdom) &
                               & * unkno(jdofn,jpoin)
                       end do
                    end do
                 end do
              end if
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(solve_sol(1) % reaction,'SUM','IN MY ZONE')
        else
           call runend('PRESOL: NOT CODED')
        end if

     end if

  end if

end subroutine presol
