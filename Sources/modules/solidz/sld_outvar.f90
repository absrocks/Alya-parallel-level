!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_outvar.f90
!> @author  Mariano Vazquez
!> @date    August, 2006
!>          - Subroutine written
!> @brief   Field Output variables for post-process of results
!>
!> @details
!>
!>          \verbatim
!>          Post-processing vectors (gevec):
!>          --------------------------------
!>          Vectors shall be stored in gevec with ndime components
!>
!>          Post-processing scalars at element/point level (gesca):
!>          -------------------------------------------------------
!>          Scalars shall be stored in gesca
!>
!>          Post-processing symmetric tensors (gevec):
!>          ------------------------------------------
!>          Tensors are stored in gevec with a voigt notation. They are
!>          defined as "TENSO" in sld_inivar. Finally, this subroutine
!>          spread the tensor components in scalars, dumped in:
!>          "...XX", "...YY", "...XY", ...
!>
!>          Post-processing scalars at gauss point level (ger3p):
!>          -----------------------------------------------------
!>          Scalars shall be stored in ger3p
!>
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_outvar(ivari,imesh)

  use def_kintyp,          only: ip,rp
  use def_parame,          only: zero
  use def_master,          only: coupling
  use def_master,          only: gesca,gevec,ger3p
  use def_master,          only: INOTMASTER, INOTEMPTY, kfl_paral
  use def_master,          only: ITER_K, ITER_K_STATE
  use def_master,          only: nfacg,npoin_type,rhsid,lfacg,displ
  use def_master,          only: postp,solve,cutim,ittim,forcf,solve_sol,lninv_loc
  use def_domain,          only: ndime,nelem,npoin,ngaus,nboun,neset,nbset
  use def_master,          only: gdepo
  use def_domain,          only: nnodf,nnode
  use def_domain,          only: lmate,lnods,ltype,lface,lnoch,lnodb,ltypb,lpoty
  use def_domain,          only: leset,lbset
  use def_domain,          only: vmass,coord
  use def_domain,          only: kfl_icodb, kfl_codbo, kfl_icodn, kfl_codno
  use def_elmtyp,          only: NODE_CONTACT_SOLID, HEX08, QUA04
  use def_solidz,          only: bvess_sld
  use def_solidz,          only: veloc_sld, accel_sld
  use def_solidz,          only: frxid_sld, fexte_sld
  use def_solidz,          only: lmate_sld,lawst_sld,lawch_sld
  use def_solidz,          only: svegm_sld,svaux_sld
  use def_solidz,          only: kfote_sld
  use def_solidz,          only: caust_sld, cause_sld, caunp_sld, caulo_sld, caunn_sld
  use def_solidz,          only: green_sld, lepsi_sld, eprin_sld
  use def_solidz,          only: epsee_sld, celen_sld
  use def_solidz,          only: dxfem_sld,vxfem_sld,axfem_sld,lnenr_sld
  use def_solidz,          only: kfl_fixno_sld, kfl_fixbo_sld, kfl_fixrs_sld
  use def_solidz,          only: jacrot_du_dq_sld
  use def_solidz,          only: lcrkf_sld,isoch_sld
  use def_solidz,          only: nvoig_sld
  use def_solidz,          only: sigei_sld,rorig_sld
  use def_solidz,          only: kfl_fiber_sld,fibeg_sld,fibde_sld,seqvm_sld
  use def_solidz,          only: kfl_cohes_sld, kfl_damag_sld
  use def_solidz,          only: axis1_sld, axis2_sld, axis3_sld, orien_sld
  use mod_sld_csys,        only: sld_csys_midsurface_element_normal, sld_csys_rotuni
  use def_solidz,          only: kfl_conta_sld, fcont_sld
  use mod_communications,  only: PAR_INTERFACE_NODE_EXCHANGE
  use mod_sld_commdom,     only: commdom_sld_sigma_dot_n
  use mod_sld_commdom,     only: commdom_sld_nominal_stress_dot_n, commdom_sld_n_sigma_n !< 2015JUL17
  use mod_commdom_driver,  only: commdom_driver_n_fixno
  use mod_commdom_driver,  only: commdom_driver_get_total_flux
  use mod_exm_sld_eccoupling, only: calcium_ecc, kfl_exmsld_3Dcou_ecc
#ifdef COMMDOM
  use mod_commdom_dynamic, only: commdom_dynamic_outvar
  use mod_commdom_dynamic, only: commdom_dynamic_kdtree
  use mod_sld_commdom,     only: commdom_sld_interior_list
#endif
  use mod_sld_fe2
  use mod_exm_sld_eccoupling, only : has_exmsld_coupling
  use mod_outvar,             only : outvar
  !
  implicit none

  integer(ip), intent(in) :: ivari      !< checking variable ivari
  integer(ip), intent(in) :: imesh   !< Mesh to postprocess on
  integer(ip)             :: pnode,pgaus,pelty,pblty,pnodb
  integer(ip)             :: ivoig,icont,ipoin,jdime
  integer(ip)             :: jpoin,idime,iline,inodf,iface,ibopo
  integer(ip)             :: ielem,ifacg,inode,ielty,idofn,igaus,isdva
  integer(ip)             :: iboun,inodb
  real(rp)                :: xauxi(2),anpoi,andis,angul,dummr,foref
  real(rp)                :: u(3),dudx(9),d2udx2(27)
  real(rp)                :: norvec(ndime)
  real(rp)                :: svegm_aux
  real(rp)                :: bvess_aux(ndime,npoin), fexte_aux(ndime,npoin)
  real(rp), pointer       :: aux(:,:) => null()
  integer(ip)             :: e
  integer(ip)             :: cost
  logical                 :: non_linear, converged

  !
  ! Define postprocess variable
  !
  select case (ivari)

  case(0_ip)
     !
     ! Do not postprocess anything
     !
     return

  case(1_ip)
     !
     ! DISPL: Displacement
     !
     if ( INOTMASTER ) gevec => displ(:,:,ITER_K)

  case(2_ip)
     !
     ! VELOC: Velocity
     !
     if ( INOTMASTER ) gevec => veloc_sld(:,:,ITER_K)

  case(3_ip)
     !
     ! ACCEL: Acceleration
     !
     if ( INOTMASTER ) gevec => accel_sld(:,:,ITER_K)

  case(4_ip)
     !
     ! SIGMA: Cauchy stress tensor
     !
     if ( INOTMASTER ) then
        kfote_sld = kfote_sld + 1_ip
        ! fotens computes caust, green and lepsi, so do it only the first kfote
        if (kfote_sld == 1_ip) call sld_fotens
        call memgen(0_ip, nvoig_sld, npoin)
        gevec(1:nvoig_sld,1:npoin) = caust_sld(1:nvoig_sld,1:npoin)
     end if

  case(5_ip)
     !
     ! EPSIL: Strain tensor (Green)
     !
     if ( INOTMASTER ) then
        kfote_sld = kfote_sld + 1_ip
        ! fotens computes caust, green and lepsi, so do it only the first kfote
        if (kfote_sld == 1_ip) call sld_fotens
        call memgen(0_ip, nvoig_sld, npoin)
        gevec(1:nvoig_sld,1:npoin) = green_sld(1:nvoig_sld,1:npoin)
     end if

  case(6_ip)
     !
     ! LNEPS: Logarithmic strain tensor
     !
     if ( INOTMASTER ) then
        kfote_sld = kfote_sld + 1_ip
        ! fotens computes caust, green and lepsi, so do it only the first kfote
        if (kfote_sld == 1_ip) call sld_fotens
        call memgen(0_ip, nvoig_sld, npoin)
        gevec(1:nvoig_sld,1:npoin) = lepsi_sld(1:nvoig_sld,1:npoin)

     end if

  case(7_ip)
     !
     ! SEQVM: Von Mises Stress
     !
     if ( INOTMASTER ) then
        kfote_sld = kfote_sld + 1_ip
        ! fotens computes caust, green and lepsi, so do it only the first kfote
        if (kfote_sld == 1_ip) call sld_fotens
        gesca => seqvm_sld
     end if

  case(8_ip)
     !
     ! Fibres
     !
     if (kfl_fiber_sld == 0) call runend('SLD_OUTVAR: CANNOT POSTPROCESS A FIBER FIELD WHEN UNDEFINED')
     call sld_fidbes()
     gevec => fibde_sld

  case(9_ip)
     !
     ! BVESS: Essential boundary conditions
     !
     if ( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        bvess_aux = 0.0_rp
        do ipoin = 1, npoin
           ibopo = lpoty(ipoin)
           bvess_aux(1:ndime,ipoin) = bvess_sld(1:ndime,ipoin,ITER_K)
           if ( ibopo > 0 ) then
              if ( kfl_fixrs_sld(ipoin) /= 0_ip .and. kfl_fixno_sld(1,ipoin) /=3_ip) then
                 ! Local --> Global
                 call sld_csys_rotuni(2_ip,ndime,ipoin,bvess_aux(1:ndime,ipoin))
             ! else if (kfl_fixno_sld(1,ipoin) == 3_ip) then
             !    ! Me interesa ver en el contato lo que aplico realmente
             !    call sld_csys_rotuni(1_ip,ndime,ipoin,bvess_aux(1:ndime,ipoin))
              end if
           end if
        end do
        gevec(1:ndime,1:npoin) = bvess_aux(1:ndime,1:npoin)
     end if

  case(10_ip)
     !
     ! Fibres on elements
     !
     call sld_elmope(4_ip)
     ger3p => fibeg_sld

  case(11_ip)
     !
     ! XFEM Displacement
     !
     gevec => dxfem_sld(:,:,1)

  case(12_ip)
     !
     ! XFEM Velocity
     !
     gevec => vxfem_sld(:,:,1)

  case(13_ip)
     !
     ! XFEM Acceleration
     !
     gevec => axfem_sld(:,:,1)

  case(14_ip)
     !
     ! LINEL: Linelets of preconditioner CG
     !
     if( INOTMASTER ) then
        icont=0
        do ipoin = 1,npoin
           rhsid(ipoin)=0.0_rp
        end do
        do iline=1,solve(1)%nline
           icont=icont+1
           do ipoin=solve(1)%lline(iline),solve(1)%lline(iline+1)-1
              jpoin=solve(1)%lrenup(ipoin)
              rhsid(jpoin)=real(icont,rp)
           end do
        end do
        gesca => rhsid
     end if

  case(15_ip)
     !
     ! Stress tensor
     !
     if (INOTMASTER) then
        kfote_sld= kfote_sld + 1
        if (kfote_sld == 1) then   ! fotens computes caust, green and lepsi, so do it only the first kfote
           call sld_fotens
        end if
     end if

     gesca => caunn_sld

  case(16_ip)
     !
     ! Rotation on plane XY
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin

           anpoi= 0.0_rp  ! coordinates angle in xy
           if (abs(coord(1,ipoin)) > 0) anpoi= atan(coord(2,ipoin)/coord(1,ipoin))
           andis= 0.0_rp  ! displacement angle in xy
           if (abs(displ(1,ipoin,1)) > 0) anpoi= atan(displ(2,ipoin,1)/displ(1,ipoin,1))
           !!acaaaaaaaaa
           angul= andis - anpoi     ! displacement angle with respect to coord

           gesca(ipoin)= angul

        end do
     end if

  case(17_ip)
     !
     ! Rotation on plane XZ
     !
     if( INOTMASTER .and. ndime==3) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           xauxi(1)= displ(3,ipoin,1) - rorig_sld(3)
           xauxi(2)= displ(2,ipoin,1) - rorig_sld(2)
           gesca(ipoin)= 0.0_rp
           if (xauxi(1) > 0.0_rp) gesca(ipoin)= atan(xauxi(2)/xauxi(1))
        end do
     end if

  case(18_ip)
     !
     ! DCOHE: Damage variable for cohesive elements
     !
     if( kfl_cohes_sld == 0_ip ) return
     if( INOTMASTER ) then
        !
        ! Averaged value from GP (Scalar per element)
        !
        call memgen(0_ip,nelem,0_ip)
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           pnode = nnode(ltype(ielem))
           if ( lawch_sld(lmate(ielem)) == 904_ip .or. &
                lawch_sld(lmate(ielem)) == 905_ip  ) then
              pgaus = pnode/2
              gesca(ielem) = sum(svegm_sld(ielem)%a(1,1:pgaus,ITER_K_STATE))/real(pgaus,rp)
           end if
        end do
     end if

  case(19_ip)
     !
     ! LNENR_sld
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(lnenr_sld(ipoin),rp)
        end do
     end if

  case(20_ip)
     !
     ! Error
     !
     if( INOTMASTER) then
        call memgen(zero,ndime,npoin)
        do ipoin = 1,npoin
           call sld_exacso(1_ip,coord(1,ipoin),u,dudx,d2udx2,dummr,dummr)
           do idime = 1,ndime
              gevec(idime,ipoin) = displ(idime,ipoin,1)-u(idime)
           end do
        end do
     end if

  case(21_ip)
     !
     ! Cracked nodes
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ifacg = 1,nfacg
           if( lcrkf_sld(ifacg) /= 0 ) then
              iface = lfacg(3,ifacg)
              ielem = lfacg(1,ifacg)
              ielty = ltype(ielem)
              do inodf = 1,nnodf(ielty) % l(iface)
                 inode = lface(ielty) % l(inodf,iface)
                 ipoin = lnods(inode,ielem)
                 gesca(ipoin) = 1.0_rp
              end do
           end if
        end do
        call pararr('SLX',NPOIN_TYPE,npoin,gesca)
        do ipoin = 1,npoin
           gesca(ipoin) = min(1.0_rp,gesca(ipoin))
        end do
     end if

  case(22_ip)
     !
     ! Groups
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real( solve_sol(1) % lgrou(ipoin) ,rp)
        end do
     end if

  case(23_ip)
     !
     ! Sigma largest eigenvalue
     !
     if (INOTMASTER) then
        kfote_sld= kfote_sld + 1
        if (kfote_sld == 1) then   ! fotens computes caust, green and lepsi, so do it only the first kfote
           call sld_fotens
        end if
     end if

     gesca => sigei_sld

  case(24_ip)
     !
     ! SDV1E: d1 (sm152)
     !
     if( INOTMASTER ) then
        !
        ! Averaged value from GP (Scalar per element)
        !
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           svegm_aux = 0.0_rp
           if ( lawst_sld(lmate(ielem)) == 152_ip ) then
              do igaus = 1,pgaus
                 svegm_aux = svegm_aux + svegm_sld(ielem)%a(4,igaus,2)
              end do
           end if
           gesca(ielem) = svegm_aux/real(pgaus,rp)
        end do
     end if

  case(25_ip)
     !
     ! SDV2E: dG (sm152)
     !
     if( INOTMASTER ) then
        !
        ! Averaged value from GP (Scalar per element)
        !
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           svegm_aux = 0.0_rp
           if ( lawst_sld(lmate(ielem)) == 152_ip ) then
              do igaus = 1,pgaus
                 svegm_aux = svegm_aux + svegm_sld(ielem)%a(5,igaus,2)
              end do
           end if
           gesca(ielem) = svegm_aux/real(pgaus,rp)
        end do
     end if

  case(26_ip)
     !
     ! SDV3E: dK (sm152)
     !
     if( INOTMASTER ) then
        !
        ! Averaged value from GP (Scalar per element)
        !
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           svegm_aux = 0.0_rp
           if ( lawst_sld(lmate(ielem)) == 152_ip ) then
              do igaus = 1,pgaus
                 svegm_aux = svegm_aux + svegm_sld(ielem)%a(6,igaus,2)
              end do
           end if
           gesca(ielem) = svegm_aux/real(pgaus,rp)
        end do
     end if

  case(27_ip)
     !
     ! SDV4E: d6 (sm152)
     !
     if( INOTMASTER ) then
        !
        ! Averaged value from GP (Scalar per element)
        !
        call memgen(zero,nelem,zero)
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           svegm_aux = 0.0_rp
           if ( lawst_sld(lmate(ielem)) == 152_ip ) then
              do igaus = 1,pgaus
                 svegm_aux = svegm_aux + svegm_sld(ielem)%a(7,igaus,2)
              end do
           end if
           gesca(ielem) = svegm_aux/real(pgaus,rp)
        end do
     end if

  case(28_ip)
     !
     ! FIXRS: Fixity code for local axes
     !
     if ( INOTMASTER ) then
        call memgen(0_ip,npoin,0_ip)
        do ipoin = 1,npoin
           if ( lpoty(ipoin) > 0 ) then
              gesca(ipoin) = real(kfl_fixrs_sld(ipoin),rp)
           else
              gesca(ipoin) = 0.0_rp
           end if
        end do
     end if

  case(29_ip)
     !
     ! SDV1G: d1 (sm152)
     !
     if( INOTMASTER ) then
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           do igaus = 1,pgaus
              if( lawst_sld(lmate_sld(ielem)) == 152_ip ) then
                 svaux_sld(ielem)%a(1,igaus,1) = svegm_sld(ielem)%a(4,igaus,2)
              else
                 svaux_sld(ielem)%a(1,igaus,1) = 0.0_rp
              end if
           end do
        end do
        ger3p => svaux_sld
     end if

  case(30_ip)
     !
     ! SDV2G: dG (sm152)
     !
     if( INOTMASTER ) then
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           do igaus = 1,pgaus
              if( lawst_sld(lmate_sld(ielem)) == 152_ip ) then
                 svaux_sld(ielem)%a(1,igaus,1) = svegm_sld(ielem)%a(5,igaus,2)
              else
                 svaux_sld(ielem)%a(1,igaus,1) = 0.0_rp
              end if
           end do
        end do
        ger3p => svaux_sld
     end if

  case(31_ip)
     !
     ! SDV3G: dG (sm152)
     !
     if( INOTMASTER ) then
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           do igaus = 1,pgaus
              if( lawst_sld(lmate_sld(ielem)) == 152_ip ) then
                 svaux_sld(ielem)%a(1,igaus,1) = svegm_sld(ielem)%a(6,igaus,2)
              else
                 svaux_sld(ielem)%a(1,igaus,1) = 0.0_rp
              end if
           end do
        end do
        ger3p => svaux_sld
     end if

  case(32_ip)
     !
     ! SDV4G: d6 (sm152)
     !
     if( INOTMASTER ) then
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           do igaus = 1,pgaus
              if( lawst_sld(lmate_sld(ielem)) == 152_ip ) then
                 svaux_sld(ielem)%a(1,igaus,1) = svegm_sld(ielem)%a(7,igaus,2)
              else
                 svaux_sld(ielem)%a(1,igaus,1) = 0.0_rp
              end if
           end do
        end do
        ger3p => svaux_sld
     end if

  case(33_ip)
     !
     ! FIXNO: Fixity codes for each DoF
     !
     if ( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        gevec(1:ndime,1:npoin) = real(kfl_fixno_sld(1:ndime,1:npoin),rp)
     end if

  case(34_ip)
     !
     ! Exact force
     !
     if( INOTMASTER) then
        call memgen(zero,ndime,npoin)
        do idofn = 1,npoin*ndime
           rhsid(idofn) = 0.0_rp
        end do
        call sld_elmope(9_ip)
        call PAR_INTERFACE_NODE_EXCHANGE(ndime,rhsid,'SUM')
        do ipoin = 1,npoin
           idofn = (ipoin-1)*ndime
           do idime = 1,ndime
              idofn = idofn + 1
              gevec(idime,ipoin) = rhsid(idofn) / vmass(ipoin)
           end do
        end do
     end if

  case(35_ip)
     !
     ! Fluid force
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        do ipoin = 1,npoin
           if( lnoch(ipoin) ==  NODE_CONTACT_SOLID ) then
              do idime = 1,ndime
                 if( kfl_fixno_sld(idime,ipoin) /= 1 ) then
                    foref = 0.0_rp
                    do jdime = 1,ndime
                       foref = foref + gdepo(idime,jdime,ipoin) * forcf(jdime,ipoin)
                    end do
                    gevec(idime,ipoin) = foref
                 end if
              end do
           end if
        end do
     end if

  case(36_ip)
     !
     ! Element-wise Cauchy stress
     !
     if( INOTMASTER ) then
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           cause_sld(ielem) % a(1:ndime,1:pgaus,1) = 0.0_rp
        end do
        call sld_elmope(10_ip)
        ger3p => cause_sld
     end if

  case(37_ip)
     !
     ! FRXID: Nodal Reactions Forces
     !
     if ( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        idofn = 0_ip
        do ipoin = 1,npoin
           do idime = 1,ndime
              idofn = idofn + 1
              gevec(idime,ipoin) = -frxid_sld(idofn) ! opposite sign
           end do
        end do
     end if

  case(38_ip)
     !
     ! Recovered stresses
     !
     if( INOTMASTER ) then
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           caunp_sld(ielem) % a(1:ndime,1:pgaus,1) = 0.0_rp
        end do
        call sld_elmope(11_ip)
        ger3p => caunp_sld
     end if

  case(39_ip)
     !
     ! Element-wise Cauchy stress according to the local coordinate system        ( *AQU* )
     !
     if( INOTMASTER ) then
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           caulo_sld(ielem) % a(:,:,1) = 0.0_rp
        end do
        call sld_elmope(12_ip)
        ger3p => caulo_sld
     end if

  case(40_ip)
     !
     ! AXIS1: (Material CSYS)
     !
     if ( INOTMASTER ) then
        call memgen(0_ip, ndime, nelem)
        if ( kfl_fiber_sld > 3 ) then
           gevec(1:ndime,1:nelem) = axis1_sld(1:ndime,1:nelem)
        else
           gevec(1:ndime,1:nelem) = 0.0_rp
        end if
     end if

  case(41_ip)
     !
     ! AXIS2: (Material CSYS)
     !
     if ( INOTMASTER ) then
        call memgen(0_ip, ndime, nelem)
        if ( kfl_fiber_sld > 3 ) then
           gevec(1:ndime,1:nelem) = axis2_sld(1:ndime,1:nelem)
        else
           gevec(1:ndime,1:nelem) = 0.0_rp
        end if
     end if

  case(42_ip)
     !
     ! AXIS3: (Material CSYS)
     !
     if ( INOTMASTER ) then
        call memgen(0_ip, ndime, nelem)
        if ( kfl_fiber_sld > 3 .and. ndime == 3_ip ) then
           gevec(1:ndime,1:nelem) = axis3_sld(1:ndime,1:nelem)
        else
           gevec(1:ndime,1:nelem) = 0.0_rp
        end if
     end if

  case(43_ip)
     !
     ! ORIEN: Orientation angle
     !
     if ( INOTMASTER ) then
        call memgen(0_ip,nelem,0_ip)
        if ( kfl_fiber_sld > 3 ) then
           gesca(1:nelem) = orien_sld(1:nelem)
        else
           gesca(1:nelem) = 0.0_rp
        end if
     end if

  case(44_ip)
     !
     ! SREAC: Solver reaction
     !
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        if ( solve_sol(1)%kfl_react == 1_ip ) then
           gevec(1:ndime,1:npoin) = solve_sol(1)%reaction(1:ndime,1:npoin)
        end if
     end if

  case ( 45_ip )
     ! 'NSIGM'
     if( INOTMASTER ) then
        call memgen(0_ip, ndime, npoin)
        call commdom_sld_sigma_dot_n( gevec(1:ndime,1:npoin) )
     endif

  case ( 46_ip )
     !
     ! PMATE: Material Numbering
     !
     if ( INOTMASTER ) then
        call memgen(0_ip,nelem,0_ip)
        gesca(1:nelem) = real(lmate(1:nelem),rp)
     end if

  case ( 47_ip )
#ifdef COMMDOM
     call commdom_dynamic_outvar( ) ! TOUCH
#else
     call memgen(0_ip,npoin,0_ip)
     gesca(1:npoin) = -69.69_rp
#endif

  case ( 48_ip )
     call memgen(0_ip,npoin,0_ip)
     gesca(1:npoin) = 0.0_rp
     call commdom_driver_n_fixno( gesca ) ! NFIXN

  case( 49_ip )
     ! 'TSIGN'
     if( INOTMASTER ) then
        allocate( aux(ndime,npoin) )
        call commdom_sld_sigma_dot_n( aux(1:ndime,1:npoin) )
        !
        call memgen(0_ip, ndime, npoin)
        do idime = 1,ndime
           call commdom_driver_get_total_flux( aux(idime,1:npoin), gevec(idime,1:npoin) )
        enddo
        !
        deallocate( aux )
     endif

  case ( 50_ip )
     ! 'DIISO'
     if( INOTMASTER ) then
        gesca => isoch_sld(:,1)
     endif

  case ( 51_ip )
     ! 'VMISO'
     if( INOTMASTER ) then
        gesca => isoch_sld(:,2)
     endif

  case ( 52_ip )
     !
     ! BOSET: Boundary set
     !
     if( nbset < 1 ) return
     if( INOTMASTER ) then
        call memgen(0_ip,npoin,0_ip)
        gesca(1:npoin) = 0.0_rp
        do iboun = 1,nboun
           do inodb = 1,nnode(abs(ltypb(iboun)))
              ipoin = lnodb(inodb,iboun)
              gesca(ipoin) = max(gesca(ipoin),real(lbset(iboun),rp))
           end do
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
     end if

  case ( 53_ip )
     !
     ! 'NOMIN'
     !
     if( INOTMASTER ) then
        call memgen(0_ip, ndime, npoin)
        call commdom_sld_nominal_stress_dot_n( vprop=gevec(1:ndime,1:npoin) )
     end if

  case ( 54_ip )
     ! 'NNSIG'
     if( INOTMASTER ) then
        call memgen(0_ip, npoin, 0_ip)
        call commdom_sld_n_sigma_n( gesca(1:npoin) )
     endif

  case ( 55_ip )
     !
     ! ELSET: Element set
     !
     if( neset < 1 ) return
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,0_ip)
        gesca(1:nelem) = real(leset(1:nelem),rp)
     end if

  case ( 56_ip )
     ! 'KDTRE'
     if( INOTMASTER ) then
        call memgen(zero, ndime, npoin)
        gevec(1:ndime,1:npoin) = 0.0_rp
#ifdef COMMDOM
        call commdom_dynamic_kdtree( gevec(1:ndime,1:npoin) )
#endif
     endif

  case( 57_ip )
     !
     ! Infinitessimal strain tensor
     !
     if (INOTMASTER) then
        !green_sld = 0.0_rp
        call sld_elmope(13_ip)
        call memgen(0.0_rp, 6_ip, npoin)
        do ipoin = 1,npoin
           do ivoig = 1,nvoig_sld
              gevec(ivoig,ipoin) = green_sld(ivoig,ipoin)
           end do
        end do
     end if

  case( 58_ip )
     !
     ! Infinitesimal stress tensor
     !
     if (INOTMASTER) then
        !caust_sld = 0.0_rp
        call sld_elmope(14_ip)
        call memgen(0.0_rp, 6_ip, npoin)
        do ipoin = 1,npoin
           do ivoig = 1,nvoig_sld
              gevec(ivoig,ipoin) = caust_sld(ivoig,ipoin)
           end do
        end do
     end if

  case( 59_ip )
     ! 'INTLI'
     if( INOTMASTER ) then
        call memgen(0_ip, npoin, 0_ip)
#ifdef COMMDOM
        call commdom_sld_interior_list( gesca(1:npoin) )
#endif
     endif

  case( 60_ip )
     !
     ! SVEGM: All State Dependent Variables (SDVs)
     !
     if ( INOTMASTER ) ger3p => svegm_sld(:)

  case( 61_ip )
     !
     ! CELEN: Characteristic element length
     !
     if ( INOTMASTER ) then
        call memgen(0_ip,nelem,0_ip)
        gesca(1:nelem) = celen_sld(1:nelem)
     end if

  case( 62_ip )
     !
     ! FIXBO: Fixity code on boundaries
     !
     if ( INOTMASTER ) then
        call memgen(0_ip,npoin,0_ip)
        if ( kfl_icodb > 0 ) then
           do iboun = 1,nboun
              pblty=ltypb(iboun)
              pnodb=nnode(pblty)
              do inodb=1,pnodb
                 ipoin=lnodb(inodb,iboun)
                 gesca(ipoin)= max(real(kfl_fixbo_sld(iboun),rp), real(gesca(ipoin),rp))
              end do
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
        else
           return
        end if
     end if

  case( 63_ip )
     !
     ! BOCOD: Boundary codes
     !
     if ( INOTMASTER ) then
        call memgen(0_ip,npoin,0_ip)
        if ( kfl_icodb > 0 ) then
           do iboun = 1,nboun
              pblty=ltypb(iboun)
              pnodb=nnode(pblty)
              do inodb=1,pnodb
                 ipoin=lnodb(inodb,iboun)
                 gesca(ipoin)= max(real(kfl_codbo(iboun),rp), real(gesca(ipoin),rp))
              end do
           end do
           call PAR_INTERFACE_NODE_EXCHANGE(gesca,'MAX','IN MY CODE')
        else
           return
        end if
     end if

  case( 64_ip )
     !
     ! ELNOR: Normal vector for HEX08 and QUA04 (Interface/contact elements)
     !
     if ( INOTMASTER ) then
        call memgen(0_ip,ndime,nelem)
        do ielem = 1, nelem
           pelty = ltype(ielem)
           pnode = nnode(abs(pelty))
           if ( pelty == HEX08 .or. pelty == QUA04 ) then
              call sld_csys_midsurface_element_normal(ielem, ndime, pnode, norvec)
              gevec(1:ndime,ielem) = norvec(1:ndime)
           else
              gevec(1:ndime,ielem) = 0.0_rp
           end if
        end do
     end if

  case( 65_ip )
     !
     ! FEXTE: External nodal forces
     !
     if ( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        fexte_aux(1:ndime,1:npoin) = reshape(fexte_sld,shape=(/ndime,npoin/))
        do ipoin = 1,npoin
           gevec(1:ndime,ipoin) = fexte_aux(1:ndime,ipoin)
        end do
     end if

  case ( 66_ip )
     ! 'WETNO'
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call sld_wetno(gesca(1:npoin))
     endif

  case ( 67_ip )
     !
     ! SBVNA: Solver bvnat
     !
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        if( solve_sol(1)%kfl_bvnat == 1_ip ) then
           gevec(1:ndime,1:npoin) = solve_sol(1)%bvnat(1:ndime,1:npoin)
        end if
     end if

  case ( 68_ip )
     !
     ! ELSET: Element sets
     !
     if( neset < 1 ) return
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,0_ip)
        gesca(1:nelem) = real(leset(1:nelem),rp)
     end if

  case ( 69_ip )
     !
     ! PARTI: MPI Partitions
     !
     if( INOTMASTER ) then
        call memgen(0_ip,nelem,0_ip)
        gesca(1:nelem) = real(kfl_paral,rp)
     end if

  case ( 70_ip)
     !
     ! elastic strains
     !
     if ( INOTMASTER ) then
  !     call memgen(zero, npoin, zero)
  !     call projec_elements_to_nodes(svegm_sld,gesca)
        if (kfote_sld == 1_ip) call sld_fotens
        call memgen(0.0_rp, 6_ip, npoin)
        gevec(1:nvoig_sld,1:npoin) = epsee_sld(1:nvoig_sld,1:npoin)
     end if

  case ( 71_ip )
     !
     ! ROTM1: Rotation matrix (1st vector)
     !
     if ( INOTMASTER ) then

        call memgen(0_ip,ndime,npoin)
        do ipoin = 1, npoin
           ibopo = lpoty(ipoin)
           if ( ibopo > 0 ) then
              gevec(1:ndime,ipoin) = jacrot_du_dq_sld(1:ndime,1,ipoin)
           else
              gevec(1:ndime,ipoin) = 0.0_rp
           end if
        end do

     end if

  case ( 72_ip )
     !
     ! ROTM2: Rotation vector (2nd vector)
     !
     if ( INOTMASTER ) then

        call memgen(0_ip,ndime,npoin)
        do ipoin = 1, npoin
           ibopo = lpoty(ipoin)
           if ( ibopo > 0 ) then
              gevec(1:ndime,ipoin) = jacrot_du_dq_sld(1:ndime,2,ipoin)
           else
              gevec(1:ndime,ipoin) = 0.0_rp
           end if
        end do

     end if

  case ( 73_ip )
     !
     ! ROTM3: Rotation vector (3rd vector)
     !
     if ( INOTMASTER .and. ndime == 3_ip ) then

        call memgen(0_ip,ndime,npoin)
        do ipoin = 1, npoin
           ibopo = lpoty(ipoin)
           if ( ibopo > 0 ) then
              gevec(1:ndime,ipoin) = jacrot_du_dq_sld(1:ndime,3,ipoin)
           else
              gevec(1:ndime,ipoin) = 0.0_rp
           end if
        end do

     end if

  case ( 74_ip )
     !
     ! FCONT: Contact force
     !
     if ( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        if( kfl_conta_sld /= 0_ip ) then
           gevec(1:ndime,1:npoin) = reshape(fcont_sld,shape=(/ndime,npoin/))
        else
           gevec(1:ndime,1:npoin) = 0.0_rp
        end if

     end if

  case ( 75_ip )
     !
     ! DAMAG: Damage variables
     ! sm152: d1, dG, dK, d6
     ! sm154: d1, d2, d6
     if( kfl_damag_sld == 0_ip ) return
     if( INOTMASTER ) then
        call memgen(0_ip,4_ip,nelem)
        do ielem = 1,nelem
           pgaus = ngaus(ltype(ielem))
           if(      lawst_sld(lmate(ielem)) == 152_ip ) then
              do isdva=1, 4
                 gevec(isdva,ielem) = sum(svegm_sld(ielem)%a(isdva+3,1:pgaus,ITER_K_STATE))/real(pgaus,rp)
              end do
           else if( lawst_sld(lmate(ielem)) == 154_ip ) then
              gevec(1,ielem) = sum(svegm_sld(ielem)%a( 5,1:pgaus,ITER_K_STATE))/real(pgaus,rp)
              gevec(2,ielem) = sum(svegm_sld(ielem)%a( 6,1:pgaus,ITER_K_STATE))/real(pgaus,rp)
              gevec(3,ielem) = sum(svegm_sld(ielem)%a(10,1:pgaus,ITER_K_STATE))/real(pgaus,rp)
              gevec(4,ielem) = sum(svegm_sld(ielem)%a(11,1:pgaus,ITER_K_STATE))/real(pgaus,rp)
           else
              gevec(:,ielem) = 0.0_rp
           end if

        end do

     else
        call memgen(0_ip,4_ip,1_ip)
     end if

  case (90_ip)

     ! MICNL
     if(INOTMASTER) then
        call memgen(zero, nelem, zero)
        gesca(1:nelem) = 0.0_rp
        do e = 1, nelem
           non_linear = fe2_is_non_linear(e)
           if (non_linear) then
                   gesca(e) = 1.0_rp
           else
                   gesca(e) = 0.0_rp
           endif
        end do
     endif

  case (91_ip)

     ! MICCO
     if(INOTMASTER) then
        call memgen(zero, nelem, zero)
        gesca(1:nelem) = 0.0_rp
        do e = 1, nelem
           cost = fe2_get_cost(e)
           gesca(e) = real(cost,rp)
        end do
     endif

  case (92_ip)

     ! MICCV
     if(INOTMASTER) then
        call memgen(zero, nelem, zero)
        gesca(1:nelem) = 1.0_rp
        do e = 1, nelem
           converged = fe2_has_converged(e)
           if (converged) then
                   gesca(e) = 1.0_rp
           else
                   gesca(e) = 0.0_rp
           endif
        end do
     endif

  case (93_ip)
     !
     ! Principal stretches
     !
     if (INOTMASTER) then
        kfote_sld= kfote_sld + 1
        if (kfote_sld == 1) then   ! fotens computes caust, green and lepsi, so do it only the first kfote
           call sld_fotens
        end if
     end if

     gesca => eprin_sld

  case (94_ip)
     ! Calcium in solidz
     !   WRITE(6,*) 'asssoc calc:', associated(calcium_ecc)
     !   WRITE(6,*) 'shape calc:', shape(calcium_ecc)
     !   WRITE(6,*) 'npoin:', npoin
     !   WRITE(6,*) 'kfl_exmsld_3Dcou_ecc:',  kfl_exmsld_3Dcou_ecc
     !   WRITE(6,*) 'INOTMASTER: ', INOTMASTER
     
     if(INOTEMPTY) then 
        call memgen(0_ip,npoin,0_ip)
        !if((( coupling('SOLIDZ','EXMEDI') >= 1 ) .or. ( coupling('EXMEDI','SOLIDZ') >= 1 ) .or. kfl_exmsld_3Dcou_ecc)) then
        if(has_exmsld_coupling() .or. kfl_exmsld_3Dcou_ecc) then
          do ipoin=1,npoin
              gesca(ipoin)= calcium_ecc(1,ipoin)
          enddo
        else
          gesca(:)= 0.0_rp
        endif
     endif


  case default

     return

  end select
  !
  ! Postprocess
  !
  call outvar(&
       ivari,&
       ittim,cutim,postp(1) % wopos(:,ivari),MESH_ID=imesh)

end subroutine sld_outvar
