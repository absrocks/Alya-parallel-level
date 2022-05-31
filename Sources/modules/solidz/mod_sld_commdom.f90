!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    mod_sld_commdom.f90
!> @author  Miguel Zavala and Matias Rivero
!> @date
!> @brief
!> @details
!
!> < 2015May08 -> add 'aux_sld'
!> < 2015May11 -> unstick
!> < 2015Jul03
!> < 2015JUL17 -> commdom_sld_nominal_stress_dot_n, commdom_sld_trace_sigma, commdom_sld_n_sigma_n
!>                n_relax -> 0.125 'on'
!> < 2015Jul30 -> commdom_sld_space_funcion
!> < 2015Ago30  > ifixbo gcc compilation error
!> < 2016MAR29 ->
!> < 2019SEP04 -> COMMDOM==2
!>
!>   --------------------------------------------------------| ITERATIONS |---!
!>   +
!>   |_Alya
!>     |_call Turnon()
!>     |_call Iniunk()
!>     |_time: do while
!>       |_call Timste()
!>       |_reset: do                            [i]             [j]
!>         |_call Begste()
!>           |_block: do while
!>             |_coupling: do while
!>               |_call Begzon()         (7)   [ +r-]          [ +r-]   (3)
!>               |_modules: do while
!>                 |_call Doiter()       (8)   [+dU-]          [ +U-]   (4)
!>                 |_call Concou()       (1)   [ +U-]          [+dU-]   (5)
!>               |_call Endzon()         (2)   [ -s+]          [ -s+]   (6)
!>             |_call Conblk()
!>       |_call Endste()
!>   |_call Turnof()
!>
!>
!> @}
!-----------------------------------------------------------------------

module mod_sld_commdom

  use def_kintyp,          only : ip, rp, lg
  use def_kintyp,          only : soltyp
  use def_master,          only : IMASTER, INOTMASTER, ISEQUEN, npoi3
  use def_master,          only : momod, modul
  use def_domain,          only : npoin, nboun, ndime, coord
#ifdef COMMDOM
  use mod_commdom_dynamic, only : commdom_dynamic_check_fixno
  use mod_commdom_dynamic, only : commdom_dynamic_set_vals
  use mod_commdom_dynamic, only : commdom_dynamic_reduce_sum
#endif
  use def_solidz,          only : frxid_sld, bvess_sld, bvnat_sld
  use def_solidz,          only : kfl_fixbo_sld, kfl_fixno_sld, ncomp_sld, aux_sld
  use mod_sld_pdn_contact, only : commdom_sld_plugin_pdn, sld_rbo_contact_plugin_rbo

  !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  implicit none
  !
  type(soltyp), pointer :: solve(:)
  !
  public  :: commdom_sld_plugin
  public  :: commdom_sld_sigma_dot_n          ! Postprocess
  public  :: commdom_sld_n_sigma_n            ! Postprocess
  public  :: commdom_sld_nominal_stress_dot_n ! Postprocess
#ifdef COMMDOM
  public  :: commdom_sld_interior_list        ! Postprocess
#endif
  !
#ifdef COMMDOM
  private :: nodes2bvnat                      !
  private :: commdom_sld_residual             !
#if COMMDOM==4
  private :: commdom_sld_set_unkno            ! plugin4
  private :: commdom_sld_set_bvnat            ! plugin4
#endif
#endif
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
contains

  !-----------------------------------------------------------------------
  !> @author
  !> @date
  !> @brief
  !> @details Subroutine plugin manager
  !-----------------------------------------------------------------------

  subroutine commdom_sld_plugin()

#ifdef COMMDOM
    use def_solidz, only : kfl_rigid_sld
    use def_solidz, only : kfl_conta_sld
    use def_solidz, only : SLD_PDN_RBO_DEFORMABLE, SLD_PDN_UNILATERAL, SLD_PDN_BILATERAL

#if  COMMDOM==2
    !
    ! PDN contact
    !
    if ( (kfl_conta_sld == SLD_PDN_UNILATERAL .or. kfl_conta_sld == SLD_PDN_BILATERAL) .and. kfl_rigid_sld == 0_ip ) then
       !
       ! Unilateral and Bilateral (Rivero Phd 2018)
       !
       call commdom_sld_plugin_pdn()

    else if ( kfl_conta_sld == SLD_PDN_RBO_DEFORMABLE ) then
       !
       ! Rigid body - Deformable
       !
       call sld_rbo_contact_plugin_rbo()

    end if
#elif  COMMDOM==4
    call commdom_sld_plugin4()
#endif
#endif

  end subroutine commdom_sld_plugin

  !-----------------------------------------------------------------------
  !> @author
  !> @date
  !> @brief
  !> @details Subroutine plugin 4
  !-----------------------------------------------------------------------

#ifdef COMMDOM
#if   COMMDOM==4
  subroutine commdom_sld_plugin4()                                               ! < 2016MAR29

    use mod_commdom_alya,    only: COMMDOM_COUPLING
    use mod_commdom_driver,  only: CNT_SMS, CNT_CPLNG, CNT_SENDRECV, commdom_driver_exchange02
    use def_master,          only: displ, unkno, title

    implicit none

    integer(ip) :: idime, ipoin, idof, icomp
    real(rp) :: d_relax = 1.0_rp
    real(rp) :: n_relax = 1.0_rp
    !
    real(rp) :: residual2(3,2) = 0.0_rp
    !
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    if( any(CNT_SENDRECV) ) then
       !-----------------------------------------------------------| code_i==1 |---!
       code_i: if(CNT_CPLNG%current_code==CNT_CPLNG%code_i) then
          if( CNT_SENDRECV(7) ) then
             !-----------------------------------------------------------| U--> |---! ! (2) SOLID --> UNKNO -->  ALEFOR
             call commdom_sld_set_unkno( CNT_CPLNG%var_ij(1:ndime,1:npoin) )
             !-----------------------------------------------------------------||---!
             !
             print *, "["//trim(title)//"] ", "SOLID <==> UNKNO/REACT", CNT_SMS
             call commdom_driver_exchange02( CNT_CPLNG, debug=.true.)                !<
             !
             !--------------------------------------------------------| dUdn<-- |---! ! (2) SOLID <-- REACT <-- NASTIN
             call commdom_sld_set_bvnat( CNT_CPLNG%var_ji(1:ndime,1:npoin), n_relax=0.1_rp, algebraic=.true. )
             !-----------------------------------------------------------------||---!
          endif
       endif code_i
       !---------------------------------------------------------------------||---!
       !----------------------------------------------------------| code_j==2 |---!
       code_j: if(CNT_CPLNG%current_code==CNT_CPLNG%code_j) then
          if( CNT_SENDRECV(7) ) then
             print *, "ERROR: ", "["//trim(title)//"]", CNT_SMS, "!!"
             stop
             !--------------------------------------------------------| dUdn--> |---!
             to_send02: if(inotmaster) then
                CNT_CPLNG%var_ij(1:ndime,1:npoin) = 0.0_rp
             end if to_send02
          endif
          !-----------------------------------------------------------------||---!
          !
          call commdom_driver_exchange02( CNT_CPLNG , debug=.true.)
          !
          !-----------------------------------------------------------| U<-- |---!
          to_recv02: if(inotmaster) then
             !call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(idime,1:npoin),     veloc(idime,1:npoin,1), relax_op=n_relax, res2=residual2(3,1:2), debug=.false. )
             CNT_CPLNG%var_ji(1:ndime,1:npoin) = 0.0_rp
          end if to_recv02
       endif code_j
       !-----------------------------------------------------------------||---!
    endif
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
  end subroutine commdom_sld_plugin4
#endif
#endif

  !----------------------------------------------------------------------------
  !> @author
  !> @date
  !> @brief   Interior list of nodes
  !> @details
  !----------------------------------------------------------------------------

#ifdef COMMDOM
  subroutine commdom_sld_interior_list( prop )

    use mod_commdom_plepp,   only: PLEPP_CPLNG
    use def_master,          only: ittim

    implicit none

    real(rp), intent(inout) :: prop(npoin)
    integer(ip)             :: ipoin

    if ( INOTMASTER .and. ittim > 0_ip ) then
       prop(:) = 0_ip
       do ipoin = 1,size(PLEPP_CPLNG%interior_list_j)
          prop(PLEPP_CPLNG%interior_list_j(ipoin)) = 1_ip
       end do
    end if

  end subroutine commdom_sld_interior_list
#endif

  !----------------------------------------------------------------------------
  !> @author
  !> @date
  !> @brief   Postprocess contact pressure
  !> @details
  !----------------------------------------------------------------------------

  subroutine commdom_sld_sigma_dot_n( prop )
    use def_domain,          only: lpoty, exnor
    use def_solidz,          only: caust_sld, nvgij_inv_sld
    use def_solidz,          only: kfl_gdepo
    use def_master,          only: gdepo
    implicit none
    real(rp), intent(inout) :: prop(ndime,npoin)
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    integer(ip) :: ipoin, ibopo, i, j, ivoig
    real(rp)    :: G_ij(ndime, ndime), sigma_dot_n(ndime), nor(ndime)
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    if( kfl_gdepo==0) then
       print *, "[commdom_sld_sigma_dot_n] ", "add NUMERICAL_TREATMENT -> ROTAT:  ON"
       call runend("EXIT!!")
    endif
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    if(inotmaster) then
       call sld_fotens()
       prop(1:ndime,1:npoin) =  0.0_rp
       !
       do ipoin = 1,npoin
          sigma_dot_n(1:ndime) = 0.0_rp
          !
          ibopo = lpoty(ipoin)
          if(ibopo > 0) then
             !
             !     | G_11 G_12 G_13 |
             ! G = | G_21 G_22 G_23 |
             !     | G_31 G_23 G_33 |
             !
             !     |  g_1  g_6  g_5 |              [ g_1  g_2  g_3  g_4  g_5  g_6]
             !   = | G_21  g_2  g_4 |  <---> G^T = [G_11 G_22 G_33 G_23 G_13 G_12]
             !     | G_31 G_23  g_3 |
             !
             G_ij = 0.0_rp
             !
             ! Stresses
             !
             do i = 1,ndime
                do j = 1,ndime
                   ivoig = nvgij_inv_sld(i, j)
                   G_ij(i,j) = caust_sld(ivoig,ipoin)
                end do
             end do
             !
             nor(1:ndime)         = matmul( gdepo(1:ndime,1:ndime,ipoin), exnor(1:ndime,1,ibopo) )
             sigma_dot_n(1:ndime) = matmul( G_ij, nor )
          endif
          !
          prop(1:ndime,ipoin) = sigma_dot_n(1:ndime)
          !
       end do
       !
    endif

  end subroutine commdom_sld_sigma_dot_n

  !----------------------------------------------------------------------------
  !> @author
  !> @date
  !> @brief   Postprocess
  !> @details
  !----------------------------------------------------------------------------

  subroutine commdom_sld_n_sigma_n(prop)

    use def_domain,          only: lpoty, exnor
    use def_solidz,          only: caust_sld, nvgij_inv_sld
    use def_master,          only: gdepo
    use def_solidz,          only: kfl_gdepo

    implicit none

    real(rp), intent(inout) :: prop(npoin)
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    integer(ip) :: ipoin, ibopo, i, j, ivoig
    real(rp)    :: G_ij(ndime,ndime), sigma_dot_n(ndime), nor(ndime)
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    if( kfl_gdepo==0) then
       print *, "[commdom_sld_n_sigma_n] ", "add NUMERICAL_TREATMENT -> ROTAT:  ON"
       call runend("EXIT!!")
    endif
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    if(inotmaster) then
       call sld_fotens()
       prop(1:npoin) =  0.0_rp
       !
       do ipoin = 1,npoin
          sigma_dot_n(1:ndime) = 0.0_rp
          !
          ibopo = lpoty(ipoin)
          if(ibopo > 0) then
             !
             !     | G_11 G_12 G_13 |
             ! G = | G_21 G_22 G_23 |
             !     | G_31 G_23 G_33 |
             !
             !     |  g_1  g_6  g_5 |              [ g_1  g_2  g_3  g_4  g_5  g_6]
             !   = | G_21  g_2  g_4 |  <---> G^T = [G_11 G_22 G_33 G_23 G_13 G_12]
             !     | G_31 G_23  g_3 |
             !
             G_ij = 0.0_rp
             !
             ! Stresses
             !
             do i = 1,ndime
                do j = 1,ndime
                   ivoig = nvgij_inv_sld(i, j)
                   G_ij(i,j) = caust_sld(ivoig,ipoin)
                end do
             end do
             !
             nor(1:ndime)         = matmul( gdepo(1:ndime,1:ndime,ipoin), exnor(1:ndime,1,ibopo) )
             sigma_dot_n(1:ndime) = matmul( G_ij, nor )
          endif
          !
          prop(ipoin) = dot_product( nor(1:ndime), sigma_dot_n(1:ndime) )
          !
       end do
       !
    endif

  end subroutine commdom_sld_n_sigma_n

  !----------------------------------------------------------------------------
  !> @author
  !> @date
  !> @brief   Postprocess nomial stress
  !> @details
  !----------------------------------------------------------------------------

  subroutine commdom_sld_nominal_stress_dot_n( vprop , sprop)

    use def_domain,          only: lpoty, exnor
    use def_solidz,          only: caust_sld, nvgij_inv_sld
    use def_master,          only: gdepo
    use def_solidz,          only: kfl_gdepo

    implicit none

    real(rp), optional, intent(inout) :: vprop(ndime,npoin)
    real(rp), optional, intent(inout) :: sprop(      npoin)

    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    integer(ip) :: ipoin, ibopo, i, j, ivoig
    real(rp)    :: G_ij(ndime,ndime)
    !
    real(rp)  :: M(ndime,ndime), Minv(ndime,ndime)
    real(rp)  :: detM
    real(rp)  :: P_ij(ndime,ndime), P_dot_n(ndime), n_P_n
    !
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    if( kfl_gdepo==0) then
       print *, "[commdom_sld_nominal_stress_dot_n] ", "add NUMERICAL_TREATMENT -> ROTAT:  ON"
       call runend("EXIT!!")
    endif
    !-----------------------------------------------------------------------||---!
    !
    !       P = J F^{-1} \cdot \sigma         Nominal stress tensor
    !       F = dx/dX                         Deformation gradient
    !       J = det(F)                        Jacobian determinant
    !  \sigma                                 Cauchy stress
    !
    !-----------------------------------------------------------------------||---!
    if(inotmaster) then
       call sld_fotens()
       !
       do ipoin = 1,npoin
          P_dot_n(1:ndime) = 0.0_rp
          n_P_n            = 0.0_rp
          !
          ibopo = lpoty(ipoin)
          if(ibopo > 0) then
             !
             G_ij = 0.0_rp
             do i = 1,ndime
                do j = 1,ndime
                   ivoig = nvgij_inv_sld(i, j)
                   G_ij(i,j) = caust_sld(ivoig,ipoin)
                end do
             end do
             !
             detM = 0.0_rp                              !! det(F)
             Minv(1:ndime,1:ndime) = 0.0_rp                              !! F^{-1}
             M(1:ndime,1:ndime) = gdepo(1:ndime,1:ndime,ipoin) !! F = dx/dX
             !
             call invmtx( M(1:ndime,1:ndime), Minv(1:ndime,1:ndime), detM, ndime)
             P_ij(1:ndime,1:ndime) = detM * matmul( Minv(1:ndime,1:ndime),  G_ij(1:ndime,1:ndime) )   !!  P  = J F^{-1} \sigma
             P_dot_n(1:ndime)      =        matmul( P_ij(1:ndime,1:ndime), exnor(1:ndime,1,ibopo) )   !!  Pn = P \cdot n
             n_P_n                 =        dot_product( exnor(1:ndime,1,ibopo), P_dot_n(1:ndime) )   !! nPn = n \cdot P  \cdot n
             !
          endif
          !
          if( present(vprop) ) vprop(1:ndime,ipoin) = P_dot_n(1:ndime)
          if( present(sprop) ) sprop(        ipoin) = n_P_n
          !
       end do
       !
    endif

  end subroutine commdom_sld_nominal_stress_dot_n

  !----------------------------------------------------------------------------
  !> @author  Matias Rivero
  !> @date
  !> @brief
  !> @details
  !----------------------------------------------------------------------------
#ifdef COMMDOM
  subroutine nodes2bvnat(prop, h_flux) !< ok
    use def_parame
    use def_master
    use def_domain
    use def_elmtyp, only: TRI03, TRI06, TET04, HEX08
    implicit none
    real(rp),   intent( in   ) :: prop(npoin)
    real(rp),   intent( inout) :: h_flux(nboun)
    !
    real(rp)    :: bprop(mnodb)
    real(rp)    :: xbprop(mgaus)
    real(rp)    :: bocod(ndime,mnodb)
    real(rp)    :: xbocod(ndime,mgaus)
    integer(ip) :: boidx(mnodb)
    integer(ip) :: igaub, iboun, pgaub, pnodb, pblty
    !
    !
    if(INOTMASTER) then
       !
       boundaries: do iboun = 1,nboun
          !
          h_flux(iboun) = 0.0_rp
          !
          if(kfl_fixbo_sld(iboun) == 6) then
             !
             pblty = ltypb(iboun)
             !
             !      tria03: &
             !      if(pblty == TRI03) then
             pnodb = nnode(pblty)
             pgaub = ngaus(pblty) !< pgaub == 1
             !
             boidx(1:pnodb)         = lnodb(1:pnodb,iboun)
             !
             bocod(1:ndime,1:pnodb) = coord(1:ndime, boidx(1:pnodb) )
             do igaub = 1,pgaub
                xbocod(1:ndime,igaub) = matmul( bocod(1:ndime,1:pnodb), elmar(pblty) % shape(1:pnodb,igaub) )
             enddo
             !
             bprop(1:pnodb)  = prop( boidx(1:pnodb) )
             do igaub = 1,pgaub
                xbprop(igaub) = dot_product( bprop(1:pnodb), elmar(pblty) % shape(1:pnodb,igaub) )
             enddo
             !
             h_flux(iboun) = sum( xbprop(1:pgaub) )/real(pgaub,rp)
             !
          endif
          !
       end do boundaries

    endif

  end subroutine nodes2bvnat
#endif

  !-----------------------------------------------------------------------
  !>
  !> @author
  !> @date
  !> @brief   plugin4
  !> @details plugin4
  !>
  !-----------------------------------------------------------------------
#ifdef COMMDOM
#if   COMMDOM==4
  subroutine commdom_sld_set_unkno( var_ij )                                   ! < 2016MAR29
    implicit none
    real(rp),               intent(inout) :: var_ij(ndime,npoin)
    !--------------------------------------------------------------| dUdn--> |---!
    if(inotmaster) then
       !
       var_ij(1:ndime,1:npoin) = 2.0_rp
       !var_ij(1:ndime,1:npoin) = displ(1:ndime,1:npoin,1) - displ(1:ndime,1:npoin,3)
       !var_ij(1:ndime,1:npoin) = displ(1:ndime,1:npoin,1) + coord(1:ndime,1:npoin)
       !
    endif
    !-----------------------------------------------------------------------||---!
  end subroutine commdom_sld_set_unkno
#endif
#endif

  !-----------------------------------------------------------------------
  !>
  !> @author
  !> @date
  !> @brief   plugin4
  !> @details plugin4
  !>
  !-----------------------------------------------------------------------
#ifdef COMMDOM
#if   COMMDOM==4
  subroutine commdom_sld_set_bvnat( var_ji, n_relax, algebraic )               ! < 2016MAR29
    use def_master,        only : gdepo
    use def_solidz,        only :  kfl_fixno_sld
    implicit none
    real(rp),               intent(inout) :: var_ji(ndime,npoin)
    real(rp),               intent(in   ) :: n_relax
    logical(lg), optional,  intent(in   ) :: algebraic
    !
    real(rp)    :: residual2(3,2) = 0.0_rp
    integer(ip) :: idime, ipoin
    !
    !--------------------------------------------------------------| dUdn<-- |---!
    if(inotmaster) then
       !
       ! U(x) = dx/dX.U(X)
       !
       do ipoin = 1,npoin
          where( kfl_fixno_sld(1:ndime,ipoin) /= 1 ) &
               var_ji(1:ndime,ipoin) = matmul( gdepo(1:ndime,1:ndime,ipoin), var_ji(1:ndime,ipoin) )
       enddo
       !
       if( present(algebraic).and.algebraic ) then
          do idime = 1,ndime
             call commdom_dynamic_set_vals( var_ji(idime,1:npoin), solve(1)%bvnat(idime,1:npoin), relax_op=n_relax, res2=residual2(idime,1:2) )
          enddo
       else
          !aux_sld(1:ndime,1:npoin,1) = aux_sld(1:ndime,1:npoin,2)
          !
          do idime = 1,ndime
             call commdom_dynamic_set_vals( var_ji(idime,1:npoin), aux_sld(idime,1:npoin,1), relax_op=n_relax, res2=residual2(idime,1:2), debug=.false. )
             call nodes2bvnat( aux_sld(idime,1:npoin,1), bvnat_sld(idime,1:nboun,1) )
          enddo
          !
          !aux_sld(1:ndime,1:npoin,2) = aux_sld(1:ndime,1:npoin,1)
          !aux_sld(1:ndime,1:npoin,1) = 0.0_rp
       endif
       !
       var_ji(1:ndime,1:npoin) = 0.0_rp
       !
    endif
    !-----------------------------------------------------------------------||---!
    call commdom_sld_residual( residual2(1:3,1:2), debug=.true.)
    !-----------------------------------------------------------------------||---!
  end subroutine commdom_sld_set_bvnat
#endif
#endif

  !----------------------------------------------------------------------------
  !> @author  Matias Rivero
  !> @date
  !> @brief   Residual
  !> @details
  !----------------------------------------------------------------------------

  subroutine commdom_sld_residual( res2, debug )

    use mod_commdom_driver, only: CNT_CPLNG
    use mod_communications, only: PAR_SUM
    use def_master,         only: title, ittim, iblok, inotslave
    use def_coupli,         only: coupling_driver_iteration, coupling_driver_max_iteration
    !use def_master,  only: kfl_gocou, title, INOTSLAVE
    use def_coupli,  only: kfl_gozon

    implicit none

    real(rp),               intent(inout) :: res2(3,2)
    logical(lg), optional,  intent(in   ) :: debug
    !
    real(rp)       :: daux(5)
    integer(ip)    :: iaux(5)
    character(128) :: saux(5)
    !
    real(rp)    ::  prop_in = 0.0_rp
    real(rp)    :: prop_out = 0.0_rp
    integer(ip) ::      dof
    !
    dof = CNT_CPLNG%n_dof
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    if( coupling_driver_max_iteration(iblok) /= 1) then
       !
       ! 'local' reduction
       call PAR_SUM(2_ip, res2(1,1:2), 'IN MY CODE')
       call PAR_SUM(2_ip, res2(2,1:2), 'IN MY CODE')
       call PAR_SUM(2_ip, res2(3,1:2), 'IN MY CODE')
       !
       res2(1:dof,2) = res2(1:dof,2)/res2(1:dof,1)
       !
       ! 'coupling' reduction  !< 2015Jul03
       prop_in  = 0.0_rp
       prop_out = 0.0_rp
       if( any(res2(1:dof,2) <= CNT_CPLNG%tolerance(1:dof)) ) prop_in = 1
#ifdef COMMDOM
       call commdom_dynamic_reduce_sum(prop_in, prop_out)
#endif
       if(prop_out > 0) then
          kfl_gozon = 0
          !!kfl_gocou = 0
          coupling_driver_iteration(iblok) = coupling_driver_max_iteration(iblok) - 1 !< !!???
       endif
       !
       if(INOTSLAVE.and.(present(debug).and.debug)) then
          iaux(1) = (ittim-1) * coupling_driver_max_iteration(iblok) + coupling_driver_iteration(iblok)
          daux(1) = res2(1,2)
          daux(2) = res2(2,2)
          daux(3) = res2(3,2)
          daux(4) = sum( res2(1:3,2) )
          write(saux(1), '(I0.11)') iaux(1)
          write(saux(2), '(E13.6)') daux(1)
          write(saux(3), '(E13.6)') daux(2)
          write(saux(4), '(E13.6)') daux(3)
          write(saux(5), '(E13.6)') daux(4)
          !
          print *, "["//trim(title)//".IT] "//trim(saux(1))//""//trim(saux(2))//""//trim(saux(3))//""//trim(saux(4))//""//trim(saux(5))
          !
       endif
       !
    endif

  end subroutine commdom_sld_residual

end module mod_sld_commdom
