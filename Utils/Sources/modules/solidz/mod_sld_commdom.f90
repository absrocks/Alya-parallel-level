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
!> < 2015JUL31 -> ommdom_sld_rotate_prop
!> < 2015Jul31 -> fixed testsuite 'TensileBar'
!>                commdom_sld_rotate_prop
!> < 2015Ago30  > ifixbo gcc compilation error
!> < 2016MAR29 ->
!> < 2016JUN02 -> (COMMDOM==2)||(COMMDOM==-2)
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
  public  :: commdom_sld_trace_sigma          ! Postprocess
  public  :: commdom_sld_nominal_stress_dot_n ! Postprocess
  public  :: commdom_sld_get_reaction         !
#ifdef COMMDOM
  public  :: commdom_sld_interior_list        ! Postprocess
#endif
  !
  private :: nodes2bvnat                      ! plugin2
  private :: commdom_sld_residual             ! plugin2
  !
  private :: commdom_sld_set_unkno            ! plugin4 (not used)
  private :: commdom_sld_set_bvnat            ! plugin4 (not used)
  !
  private :: commdom_sld_space_funcion        ! Not used

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

    !#if   (COMMDOM==2)||(COMMDOM==-2)
    !  implicit none
    !  integer(ip) :: pdn_flag
    !  pdn_flag = 0_ip !temporary fix
    !
    !  if (pdn_flag == 0_ip) then
    !    call commdom_sld_plugin2()
    !  else if (pdn_flag == 1_ip) then
    !    call commdom_sld_plugin_pdn()
    !  end if
#if    COMMDOM==-2
    call commdom_sld_plugin2()
#elif  COMMDOM==2
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

  !-----------------------------------------------------------------------
  !> @author
  !> @date
  !> @brief
  !> @details Subroutine plugin -2
  !-----------------------------------------------------------------------

#ifdef COMMDOM
#if   COMMDOM==-2
  subroutine commdom_sld_plugin2()
    use mod_commdom_alya,    only: COMMDOM_COUPLING
    use mod_commdom_driver,  only: CNT_SENDRECV, CNT_SMS
    use mod_commdom_driver,  only: CNT_CPLNG, commdom_driver_exchange02
    use def_master,          only: displ, unkno
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
       !
       if (CNT_CPLNG%current_code==CNT_CPLNG%code_i) then
          if ( CNT_SENDRECV(7) ) then
             !-----------------------------------------------------------| U--> |---!
             if (inotmaster) then
                CNT_CPLNG%var_ij(1:ndime,1:npoin) = 0.0_rp
                CNT_CPLNG%var_ij(1:ndime,1:npoin) = displ(1:ndime,1:npoin,1) + coord(1:ndime,1:npoin) !< 2015May08
             endif
             !-----------------------------------------------------------------||---!
             !
             call commdom_driver_exchange02( CNT_CPLNG, debug=.false.)               !< 2015Jul03 + debug
             !
             !--------------------------------------------------------| dUdn<-- |---!
             if(inotmaster) then
                !
                if( CNT_CPLNG%current_what(1_ip) ) then
                   !n_relax = 0.125_rp                                                    !< 2015JUL17
                   call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(1,1:npoin), solve(1)%bvnat(1,1:npoin), relax_op=n_relax, res2=residual2(1,1:2) )
                   call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(2,1:npoin), solve(1)%bvnat(2,1:npoin), relax_op=n_relax, res2=residual2(2,1:2) )
                   call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(3,1:npoin), solve(1)%bvnat(3,1:npoin), relax_op=n_relax, res2=residual2(3,1:2) )
                else                                                                 !< 2015JUL17
                   n_relax = 0.125_rp
                   !
                   aux_sld(1:ndime,1:npoin,1) = aux_sld(1:ndime,1:npoin,2)
                   !
                   call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(1,1:npoin), aux_sld(1,1:npoin,1), relax_op=n_relax, res2=residual2(1,1:2), debug=.false. )
                   call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(2,1:npoin), aux_sld(2,1:npoin,1), relax_op=n_relax, res2=residual2(2,1:2) )
                   call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(3,1:npoin), aux_sld(3,1:npoin,1), relax_op=n_relax, res2=residual2(3,1:2) )
                   !
                   call nodes2bvnat( aux_sld(1,1:npoin,1), bvnat_sld(1,1:nboun,1) )
                   call nodes2bvnat( aux_sld(2,1:npoin,1), bvnat_sld(2,1:nboun,1) )
                   call nodes2bvnat( aux_sld(3,1:npoin,1), bvnat_sld(3,1:nboun,1) )
                   !
                   aux_sld(1:ndime,1:npoin,2) = aux_sld(1:ndime,1:npoin,1)
                   aux_sld(1:ndime,1:npoin,1) = 0.0_rp
                endif
                CNT_CPLNG%var_ji(1:ndime,1:npoin) = 0.0_rp
                !
             endif
             !-----------------------------------------------------------------||---!
             call commdom_sld_residual( residual2(1:3,1:2), debug=.true.)           !< 2015Jul03
             !-----------------------------------------------------------------||---!
          endif
       endif
       if(CNT_CPLNG%current_code==CNT_CPLNG%code_j) then
          if( CNT_SENDRECV(7) ) then
             !--------------------------------------------------------| dUdn--> |---!
             if(inotmaster) then
                !
                CNT_CPLNG%var_ij(1:ndime,1:npoin) = 0.0_rp
                if( CNT_CPLNG%current_what(1_ip) ) then
                   CNT_CPLNG%var_ij(1:ndime,1:npoin) = solve(1)%reaction(1:ndime,1:npoin)
                else
                   !call commdom_sld_sigma_dot_n( CNT_CPLNG%var_ij(1:ndime,1:npoin) )                 !< 2015JUL17
                   call commdom_sld_nominal_stress_dot_n( vprop=CNT_CPLNG%var_ij(1:ndime,1:npoin) )  !< 2015JUL17
                   !
                   !call commdom_sld_rotate_prop( CNT_CPLNG%var_ij(1:ndime,1:npoin), n=.true., t=.false. ) !< 2015JUL31
                   !
                   CNT_CPLNG%var_ij(1:ndime,1:npoin) = -CNT_CPLNG%var_ij(1:ndime,1:npoin)            !< 2015JUL17
                endif
                !
             endif
             !-----------------------------------------------------------------||---!
             !
             call commdom_driver_exchange02( CNT_CPLNG , debug=.false.)
             !
             !-----------------------------------------------------------| U<-- |---!
             if(inotmaster) then
                !
                !                     Dirichlet nodes : fixval==0 ---V
                !
                ! solve(1) % kfl_fixno => kfl_fixno_sld
                !
                call commdom_dynamic_check_fixno(kfl_fixno_sld, 1_ip, 1_ip, .True.) ! fixno, idofn, fixval, ToDo
                call commdom_dynamic_check_fixno(kfl_fixno_sld, 2_ip, 1_ip, .True.) ! fixno, idofn, fixval, ToDo
                call commdom_dynamic_check_fixno(kfl_fixno_sld, 3_ip, 1_ip, .True.) ! fixno, idofn, fixval, ToDo
                !---------------------------------------------------------| U<-- |---!
                !
                CNT_CPLNG%var_ji(1:ndime,1:npoin) = CNT_CPLNG%var_ji(1:ndime,1:npoin) - coord(1:ndime,1:npoin)   !< 2015May08
                !  where( -0.9*CNT_CPLNG%var_ji(1:ndime,1:npoin) < bvess_sld(1:ndime,1:npoin,1) ) CNT_CPLNG%var_ji(1:ndime,1:npoin) = 0.0_rp ! < 2015May11 -> unstick
                !
                !call commdom_sld_rotate_prop( CNT_CPLNG%var_ji(1:ndime,1:npoin),  n=.true., t=.false. ) !< 2015JUL31
                !
                call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(1,1:npoin), bvess_sld(1,1:npoin,1), relax_op=d_relax )
                call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(2,1:npoin), bvess_sld(2,1:npoin,1), relax_op=d_relax )
                call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(3,1:npoin), bvess_sld(3,1:npoin,1), relax_op=d_relax, debug=.false.)
                !
                where(aux_sld(1:ndime,1:npoin,1) < bvess_sld(1:ndime,1:npoin,1)) &
                     aux_sld(1:ndime,1:npoin,1) = bvess_sld(1:ndime,1:npoin,1)
                !
             endif
             !-----------------------------------------------------------------||---!
             call commdom_sld_residual( residual2(1:3,1:2), debug=.false.)           !< 2015Jul03
             !-----------------------------------------------------------------||---!
          endif
       endif
       !
    endif
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
  end subroutine commdom_sld_plugin2
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
  !> @author  Matias Rivero
  !> @date
  !> @brief   Postprocess trace sigma
  !> @details
  !----------------------------------------------------------------------------

  subroutine commdom_sld_trace_sigma( prop )

    use def_domain,          only: lpoty
    use def_solidz,          only: caust_sld, nvgij_inv_sld

    implicit none

    real(rp), intent(inout) :: prop(npoin)

    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    integer(ip) :: ipoin, ibopo, i, j, ivoig
    real(rp)    :: G_ij(ndime,ndime)
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    if(inotmaster) then
       call sld_fotens()
       prop(1:npoin) =  0.0_rp
       !
       do ipoin = 1,npoin
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
          endif
          !
          ! -p = 1/3 sigma_ii = 1/3 trace( sigma )
          prop(ipoin) = -1.0_rp/3.0_rp * sum( (/ (G_ij(i,i), i=1,ndime) /) )
          !
       end do
       !
    endif

  end subroutine commdom_sld_trace_sigma

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

  subroutine commdom_sld_space_funcion( bou_prop, fixbo, idime )

    use def_parame
    use def_master
    use def_domain
    use def_kermod !, only : number_space_time_function, space_time_function
    use def_solidz
    use mod_ker_space_time_function

    implicit none

    real(rp),    intent(inout) ::  bou_prop( nboun )
    integer(ip), intent(in   ) ::  fixbo
    integer(ip), intent(in   ) ::  idime
    !
    integer(ip) :: igaub, iboun, pgaub, pnodb, pblty
    integer(ip) :: boidx(mnodb)
    real(rp)    :: bocod(ndime,mnodb)
    real(rp)    :: xbocod(ndime,mgaus)
    !
    integer(ip) :: ifunc, n_funbo
    real(rp)    :: y(mgaus), y_average

    integer(ip) :: idata
    real(rp)    :: t1, p1, t2, p2, b, m

    real(rp)    :: x

    logical(lg) :: ifixbo, codes_on_boundaries = .false.

    codes_on_boundaries = kfl_icodb > 0  !< 2015Jul31

    IF01: if(INOTMASTER.and.codes_on_boundaries) then
       !
       !< *.sld.dat
       !
       ! BOUNDARY_CONDITIONS, TRANSIENT
       !   ...                 |
       !   kfl_conbc_sld = 0 _ /
       !   ...
       !
       IF02: if( kfl_conbc_sld == 0 ) then

          !------------------------------------------------------------------------!
          n_funbo = sum( kfl_funbo(1:nboun) )
          !------------------------------------------------------------------------!
          !
          !<  *.ker.dat
          !
          ! SPACE_&_TIME_FUNCTION
          !   FUNCTION = XXX
          !    1 + x*y
          !   END_FUNCTION
          ! END_SPACE_&_TIME_FUNCTION
          !
          !< *.sld.dat
          ! CODES, BOUNDARIES
          !
          !   1   3  0.0  0.0  0.0  1, SPACE_TIME_FUNCTION=XXX
          !
          if( number_space_time_function > 0 ) then
             !----------------------------------------------------------------------!
             do iboun = 1,nboun
                ifunc     = kfl_funbo(iboun)
                y_average = 0.0_rp
                !
                if( ifunc < 0 ) then                                                !< ifunc < 0 FUNTION FROM ker ??
                   !
                   pblty = ltypb(iboun)
                   pnodb = nnode(pblty)
                   pgaub = ngaus(pblty)
                   !
                   boidx(1:pnodb)         = lnodb(1:pnodb,iboun)
                   bocod(1:ndime,1:pnodb) = coord(1:ndime, boidx(1:pnodb) )
                   !
                   y(1:mgaus) = 0.0_rp
                   do igaub = 1,pgaub
                      xbocod(1:ndime,igaub) = matmul( bocod(1:ndime,1:pnodb), elmar(pblty) % shape(1:pnodb,igaub) )
                      !--------------------------------------------------------------!
                      call ker_space_time_function(ifunc=-ifunc, x=xbocod(1:ndime,igaub), t=cutim, xvalu=y(igaub) )
                      !--------------------------------------------------------------!
                   enddo
                   !
                   y_average = sum( y(1:pgaub) )/real(pgaub,rp)
                   bou_prop(iboun) = y_average
                   !
                endif
             enddo
             !----------------------------------------------------------------------!
          endif
          !
          !< *.sld.dat
          ! BOUNDARY_CONDITIONS, TRANSIENT
          !  ...
          ! CODES, BOUNDARIES
          !  1   3  0.0  0.0  0.0  1
          !  ...
          !
          ! FUNCTIONS
          !   TOTAL_NUMBER = 1
          !    CONDITION
          !      FUNCTION_NUMBER: 1
          !      TIMESHAPE: DISCRETE
          !        SHAPEDEFINITION
          !   4
          !   0.0  0.0
          !   0.090909090909089 -90.756302521
          !   3.636363636363635 -88.4359562007
          !   7.09090909090909  -86.1176470588
          !        END_SHAPEDEFINITION
          !    END_CONDITION
          ! END_FUNCTIONS
          !----------------------------------------------------------------------!
          do iboun = 1,nboun
             call runend('QUE ES ESO?') ! he comentado la linea siguiente
             ifixbo = kfl_fixbo(iboun) == fixbo
             ifunc  = kfl_funbo(iboun)
             !
             call runend('QUE ES ESO?') ! he comentado la linea siguiente
             if( (ifunc>0).and.ifixbo ) then                                                  !< ifunc > 0 FUNTION FROM sld ??
                pblty = ltypb(iboun)
                pnodb = nnode(pblty)
                pgaub = ngaus(pblty)
                !
                boidx(1:pnodb)         = lnodb(1:pnodb,iboun)
                bocod(1:ndime,1:pnodb) = coord(1:ndime, boidx(1:pnodb) )
                !
                y(1:mgaus) = 0.0_rp
                do igaub = 1,pgaub
                   xbocod(1:ndime,igaub) = matmul( bocod(1:ndime,1:pnodb), elmar(pblty) % shape(1:pnodb,igaub) )
                   !--------------------------------------------------------------!
                   x = xbocod(idime,igaub)
                   do idata = 1,mtloa_sld(ifunc)-1
                      t1 = tload_sld(ifunc)%a(ndime+1,idata  )
                      t2 = tload_sld(ifunc)%a(ndime+1,idata+1)
                      if(x>= t1) then
                         if(x < t2) then
                            p1 = tload_sld(ifunc)%a(1,idata  )
                            p2 = tload_sld(ifunc)%a(1,idata+1)
                            m  = (p2-p1)/(t2-t1)
                            b  = p1 - m * t1
                            y(igaub) = m * x + b
                         end if
                      end if
                   end do
                   !--------------------------------------------------------------!
                enddo
                !
                y_average = sum( y(1:pgaub) )/real(pgaub,rp)
                bou_prop(iboun) = y_average
                !
                !call runend('QUE ES ESO?') ! he comentado la linea siguiente
             endif
          enddo
          !----------------------------------------------------------------------!
          !
       endif IF02
       !
    endif IF01

  end subroutine commdom_sld_space_funcion

  !----------------------------------------------------------------------------
  !> @author  Matias Rivero
  !> @date
  !> @brief
  !> @details
  !----------------------------------------------------------------------------

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

  subroutine commdom_sld_set_bvnat( var_ji, n_relax, algebraic )               ! < 2016MAR29
    use def_master,        only : gdepo
    use def_solidz,        only :  kfl_fixno_sld
    implicit none
    real(rp),               intent(inout) :: var_ji(ndime,npoin)
    real(rp),               intent(in   ) :: n_relax
    logical(lg), optional,  intent(in   ) :: algebraic
    !
    real(rp) :: residual2(3,2) = 0.0_rp
    integer(ip) :: idime, ipoin
    !
#ifdef COMMDOM
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
#endif
    call commdom_sld_residual( residual2(1:3,1:2), debug=.true.)
    !-----------------------------------------------------------------------||---!
  end subroutine commdom_sld_set_bvnat

  subroutine commdom_sld_get_reaction( var_ij, algebraic )                     ! < 2016MAR29
    implicit none
    real(rp),               intent(inout) :: var_ij(ndime,npoin)
    logical(lg), optional,  intent(in   ) :: algebraic
    !--------------------------------------------------------------| dUdn--> |---!
    if(inotmaster) then
       !
       var_ij(1:ndime,1:npoin) = 0.0_rp
       if( present(algebraic).and.algebraic ) then
          var_ij(1:ndime,1:npoin) = solve(1)%reaction(1:ndime,1:npoin)
       else
          call commdom_sld_nominal_stress_dot_n( vprop=var_ij(1:ndime,1:npoin) )
          var_ij(1:ndime,1:npoin) = -var_ij(1:ndime,1:npoin)
       endif
       !
    endif
    !-----------------------------------------------------------------------||---!
  end subroutine commdom_sld_get_reaction

  subroutine commdom_sld_rotate_prop( prop, n, t)
    use def_domain,          only: lpoty, exnor
    use def_master,          only: gdepo
    use def_solidz,          only: kfl_gdepo
    implicit none
    real(rp),              intent(inout) :: prop(ndime,npoin)
    logical(lg), optional, intent(in   ) :: n
    logical(lg), optional, intent(in   ) :: t
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    integer(ip) :: ipoin, ibopo
    real(rp)    :: v_old(ndime), v_new(ndime)
    real(rp)    :: Tref(ndime,ndime), Tdef(ndime,ndime)
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
       !
       do ipoin = 1,npoin
          !
          ibopo = lpoty(ipoin)
          if(ibopo > 0) then
             !
             !      REF     ->      DEF
             !  Tref==exnor       Tdef = Tref.F
             !
             Tref(1:ndime,1:ndime) = exnor(1:ndime,1:ndime,ibopo)
             Tdef(1:ndime,1:ndime) = matmul( gdepo(1:ndime,1:ndime,ipoin), Tref(1:ndime,1:ndime) )
             !
             !   DEF                        CURV
             !  Gn = <Gn_x,Gn_y,Gn_z>     Gn = <Gn_n,Gn_t1,Gn_t2>
             !
             if(present(n).or.present(t)) then
                !
                ! Cartesian   -> Curvilinear
                v_new(1:ndime) = 0.0_rp
                v_old(1:ndime) = prop(1:ndime,ipoin)
                call sld_rotvec(-1_ip, v_old, Tref(1:ndime,1:ndime), v_new, ndime)
                !
                ! Curvilinear -> Cartesian
                v_old(1:ndime) = v_new(1:ndime)
                v_new(1:ndime) = 0.0_rp
                if(present(t).and.t) v_old(1)       = 0.0_rp                     !<- Gn = <   0,Gn_t1,Gn_t2>
                if(present(n).and.n) v_old(2:ndime) = 0.0_rp                     !<- Gn = <Gn_n,    0,    0>
                call sld_rotvec( 1_ip, v_old, Tref(1:ndime,1:ndime), v_new, ndime)
                !
                prop(1:ndime,ipoin) = v_new(1:ndime)
                !
             endif
             !
          endif
          !
       end do
       !
    end if

  end subroutine commdom_sld_rotate_prop

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
