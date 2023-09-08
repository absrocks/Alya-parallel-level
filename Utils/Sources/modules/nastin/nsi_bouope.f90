!-----------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    nsi_bouope.f90
!> @author  Guillaume Houzeaux
!> @brief   Matrix assembly: boundary contribution
!> @details Boundary operations
!>          - Element matrix calculation
!>          - Scatter in global matrix
!>          Used variables:
!>          LNNOB(1:NBOUN) ........... Number of boundary nodes for each IBOUN
!>          KFL_FIXBO(1:NBOUN) ....... Condition type (e.g.: 3 is wall law)
!>          LTYPB(1:NBOUN) ........... Boundary type
!>          LNODB(1:MNODB,1:NBOUN) ... Boundary connectvity
!>          LBOEL(MNODB,1:NBOUN) ..... Local boundary to local element numbering
!>
!>                          45        77
!>                           o---------o  Example: PNODE = 2
!>                           |         |  LNODS(1:PNODE,IELEM) = 77,45,6,12
!>                           |  IELEM  |  LNODB(1:PNODB,IBOUN) = 6,12
!>                           |         |
!>                           o--IBOUN--o----------O
!>                           6 \\\\\\\ 12 \\\\\\\\\
!>
!>                           1,2,3,4 = Local numbering
!>
!>                           2         1
!>                           o---------o  Example: PNODB = 2
!>                           |         |  LBOEL(1:PNODB,IBOUN) = 3,4
!>                           |  IELEM  |  LBOEL(PNODB+1,IBOUN) = IELEM
!>                           |         |  LNODB(1:PNODB,IBOUN) =
!>                           o--IBOUN--o----------O
!>                           3 \\\\\\\ 4 \\\\\\\\\
!>
!>         ITASK = 0 ... Default behavior
!>               = 1 ... Calculates traction on nodes for postprocessing
!> @}
!-----------------------------------------------------------------------
subroutine nsi_bouope(itask)
  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_nastin
  use mod_ker_proper
  use mod_nsi_assembly_global_system, only : nsi_assembly_fractional_step_boundary_scalar

  implicit none

  integer(ip), intent(in) :: itask

  real(rp)    :: elmat(nevat_nsi,nevat_nsi)                    ! Element matrices
  real(rp)    :: elmap(mnode,mnode)
  real(rp)    :: wmatr(nevat_nsi,nevat_nsi)
  real(rp)    :: elrhs(nevat_nsi)
  real(rp)    :: wrhsi(nevat_nsi)

  real(rp)    :: baloc(ndime,ndime)                            ! Gather
  real(rp)    :: bovel(ndime,mnodb)
  real(rp)    :: bovfi(ndime,mnodb)
  real(rp)    :: bocod(ndime,mnodb)
  real(rp)    :: bopre(mnodb)
  real(rp)    :: bomut(mnodb)
  real(rp)    :: botem(mnodb)
  real(rp)    :: elmas(mnodb)
  real(rp)    :: elcod(ndime,mnode)
  real(rp)    :: elvel(ndime,mnode)
  real(rp)    :: elpre(mnode)

  real(rp)    :: gbcar(ndime,mnode,mgaus)
  real(rp)    :: gpcar(ndime,mnode,mgaus)
  real(rp)    :: xjaci(9),xjacm(9),shapp(mnode)
  real(rp)    :: gbsur(mgaub),eucta,detjm                      ! Values at Gauss points
  real(rp)    :: gbden(mgaub)
  real(rp)    :: gbvis(mgaub)
  real(rp)    :: gbpor(mgaub)
  real(rp)    :: gbmut(mgaub)
  real(rp)    :: gbgvi(ndime,mgaub)
  real(rp)    :: grvis(ndime,mgaub)
  real(rp)    :: gbtem(mgaub)
  real(rp)    :: gbgve(ndime,ndime,mgaub)
  real(rp)    :: gbpre(mgaub)
  real(rp)    :: gbfle(mgaub)
  real(rp)    :: gpvol(mgaus)
  real(rp)    :: tract(3),chale(3),chave(3),dummr
  real(rp)    :: velex(3),veave(3)
  real(rp)    :: gpvis,ustar,tragl(9),hleng(3),kinen
  real(rp)    :: udotn,gbvel(3),roughness,dtaux
  integer(ip) :: ielem,ipoin,igaus,inode,idime,jnode,jdime     ! Indices and dimensions
  integer(ip) :: pnode,pgaus,iboun,igaub,inodb
  integer(ip) :: pelty,pmate,pblty,pnodb,pgaub
  integer(ip) :: pevat,ievat,jevat,porde,ndim1
  integer(ip) :: dummi,iflow
  !
  ! Adjoint problem
  !
  real(rp)    :: baloc_der(ndime,ndime,ndime,mnodb)
  real(rp)    :: eucta_der(ndime,mnodb)
  real(rp)    :: gbsur_der(mgaub,ndime,mnodb)
  real(rp)    :: tmpfp_der(3,ndime,mnodb),tmpfv_der(3,ndime,mnodb)
  real(rp)    :: setfp_derx(ndime,ndime,mnodb),setfp_derp(ndime,mnodb)
  real(rp)    :: eldcost_dx(ndime*mnode),elmat_aux(ndime*mnode,ndime*mnode)
  integer(ip) :: jnodb,pevat1
  logical(lg) :: if_boundaries
  !
  ! If Dynamic pressure cond. exist, compute average dynamic pressure
  !
  if( ittot_nsi >= itebc_nsi ) call nsi_bouave()
  !
  ! If stable outlfow condition is present
  !
  if( kfl_exist_fib20_nsi == 1 .and. itask /= 1_ip) call nsi_outflow_mass()


  if( IMASTER ) return
  !
  ! Initialize nodal traction and boundary mass
  !
  if( itask == 1_ip ) then
     if( associated(notra_nsi) ) notra_nsi = 0.0_rp
     if( associated(massb_nsi) ) massb_nsi = 0.0_rp
  end if
  !
  ! Check if loop should be carried out for general option
  !
  if_boundaries = .false.
  if( bemol_nsi > 0.0_rp ) if_boundaries = .true.              ! Convective term integrated by parts
  if( abs(fcons_nsi) > 0.0_rp                     .and. &      ! Adding this part looks a bit untidy but it can save time!!!
       ( kfl_stabi_nsi /= NSI_SPLIT_OSS           .and. &             
       & kfl_stabi_nsi /= NSI_GALERKIN            .and. &
       & kfl_stabi_nsi /= NSI_ALGEBRAIC_SPLIT_OSS       ) ) &
       if_boundaries = .true.
  if( kfl_regim_nsi == 3 ) if_boundaries = .true.              ! Low Mach
  !
  ! Loop over boundaries
  !
  !
  ! OpenMP declarations
  !
  !$OMP  PARALLEL DO SCHEDULE ( STATIC )                                           &
  !$OMP  DEFAULT      ( NONE )                                                     &
  !$OMP  FIRSTPRIVATE ( if_boundaries )                                            &
  !$OMP  PRIVATE      ( iboun, pblty, pnodb, ielem, pelty, pnode, pgaub, velex,    &
  !$OMP                 pgaus, porde, pevat, pmate, ievat, jevat, elmat, elrhs,    &
  !$OMP                 inode, ipoin, elcod, elvel, inodb, bovel, bovfi, bopre,    &
  !$OMP                 bomut, botem, tragl, hleng, gbcar, gbden, gbvis, veave,    &
  !$OMP                 gbgvi, gbmut, grvis, gbgve, gbpor, wmatr, wrhsi, tract,    &
  !$OMP                 bocod, baloc, eucta, kinen, ustar, dummr, shapp, xjacm,    &
  !$OMP                 xjaci, detjm, gpcar, gbpre, chave, chale, gpvol, gbfle,    &
  !$OMP                 jdime, jnodb, elpre, baloc_der, eucta_der, gbsur_der,      &
  !$OMP                 tmpfp_der,tmpfv_der, setfp_derx, setfp_derp, idime,        &
  !$OMP                 eldcost_dx,pevat1, elmat_aux,roughness,ndim1,              &
  !$OMP                 elmap, dummi, gbtem, gbsur, gpvis, gbvel, udotn, iflow,    &
  !$OMP                 elmas, igaub, jnode )                                      &
  !$OMP  SHARED       ( kfl_fixbo_nsi, bemol_nsi, kfl_regim_nsi,lelbo,             &
  !$OMP                 ltypb, nnode, lboel, ltype, ngaus, lorde, solve, itask,    &
  !$OMP                 lmate, nevat_nsi, veloc, press, coord, lnodb, nboun,       &
  !$OMP                 velom, kfl_coupl, hnatu, kfl_prope, kfl_cotur_nsi,         &
  !$OMP                 kfl_kemod_ker, ndimb, kfl_rough, rough_dom, dtinv_nsi,     &
  !$OMP                 rough, untur, kfl_ustar, kfl_waexl_ker, velel_ker, dtaux,  &
  !$OMP                 lexlo_ker, bvnat_nsi, bpess_nsi, gravi_nsi, grnor_nsi,     &
  !$OMP                 prthe, lowtr_nsi, gasco, cutim, lapla_nsi, nmate,          &
  !$OMP                 ivari_nsi, lnods, elmar, kfl_addpr_nsi, kfl_ellsh_nsi,     &
  !$OMP                 kfl_advec_nsi, kfl_matdi_nsi, nsi_schur_complement,        &
  !$OMP                 poauu_nsi, poaup_nsi, poapu_nsi, poapp_nsi, kfl_predi_nsi, &
  !$OMP                 kfl_adj_prob, kfl_cost_type, veloc_forw, dcost_dx_nsi,     &
  !$OMP                 fcons_nsi, outflow_mass , kfl_hydro_nsi, kfl_wlaav_ker,    &
  !$OMP                 NSI_FRACTIONAL_STEP, ndbgs_nsi, kfl_stabi_nsi, kfl_logva,  &
#ifndef NDIMEPAR
  !$OMP                 ndime,                                                     &
#endif
  !$OMP                 velav_ker, amatr, pmatr, rhsid, notra_nsi, massb_nsi,      &
  !$OMP                 btrac_nsi,kfl_noslw_ker   ) 

  boundaries: do iboun = 1,nboun

     if(    if_boundaries               .or.  &    ! Previously defined
          & kfl_fixbo_nsi(iboun)  ==  2 .or.  &    ! Pressure (in bvnat_nsi)
          & kfl_fixbo_nsi(iboun)  ==  3 .or.  &    ! Wall law
          & kfl_fixbo_nsi(iboun)  == 13 .or.  &    ! Wall law + open pressure
          & kfl_fixbo_nsi(iboun)  ==  5 .or.  &    ! Dynamic pressure
          & kfl_fixbo_nsi(iboun)  ==  6 .or.  &    ! Open flow
          & kfl_fixbo_nsi(iboun)  == 10 .or.  &    ! Dynamic pressure + open flow
          & kfl_fixbo_nsi(iboun)  == 11 .or.  &    ! Dynamic pressure + open flow
          & kfl_fixbo_nsi(iboun)  == 12 .or.  &    ! Nodal pressure (in bpess_nsi)
          & kfl_fixbo_nsi(iboun)  == 15 .or.  &    ! Outflow with pressure dependent on density (Low Mach regime)
          & kfl_fixbo_nsi(iboun)  == 17 .or.  &    ! Outflow with pressure dependent on density (Low Mach regime)
          & kfl_fixbo_nsi(iboun)  == 18 .or.  &    ! u.n in weak form
          & kfl_fixbo_nsi(iboun)  == 19 .or.  &    ! Impose traction
          & kfl_fixbo_nsi(iboun)  == 20 .or.  &    ! Stable outflow (Bazilevs and WK0)
          & kfl_fixbo_nsi(iboun)  == 21 .or.  &    ! Stable outflow (Bazilevs and WK1 and WK2)
          & kfl_fixbo_nsi(iboun)  == 22 .or.  &    ! Boundary traction is imposed from an auxiliary RANS simulation (Two-layer wall model)
!  beware 23 is used to identify boundaries where wall law will be applied with additional viscosity and no slip - it must not enter here          
          & (kfl_adj_prob         == 1  .and. &
          &  kfl_cost_type        == 5  .and. &
          &  kfl_fixbo_nsi(iboun) == 1) ) then     ! Adjoint terms: dF/dU and dF/dX
        !
        ! Element properties and dimensions
        !
        pblty = ltypb(iboun)
        pnodb = nnode(pblty)
        ielem = lelbo(iboun)
        pelty = ltype(ielem)

        if( pelty > 0 ) then

           pnode = nnode(pelty)
           pgaub = ngaus(pblty)
           pgaus = ngaus(pelty)
           porde = lorde(pelty)
           ndim1 = ndime + 1
           pevat = ndim1 * pnode
           pmate = 1

           if( nmate > 1 ) then
              pmate = lmate(ielem)
           end if

           if( pmate /= -1 ) then
              !
              ! Initialize
              !
              do ievat = 1,nevat_nsi
                 do jevat = 1,nevat_nsi
                    elmat(jevat,ievat) = 0.0_rp
                 end do
                 elrhs(ievat) = 0.0_rp
              end do
              elmas(1:pnodb) = 0.0_rp
              eldcost_dx     = 0.0_rp
              elmat_aux      = 0.0_rp
              !
              ! Gather operations: ELVEL, ELCOD, BOVEL
              !
              do inode = 1,pnode
                 ipoin = lnods(inode,ielem)
                 if (kfl_adj_prob == 0) then
                   elvel(1:ndime,inode) = veloc(1:ndime,ipoin,1)
                 else
                   elvel(1:ndime,inode) = veloc_forw(1:ndime,ipoin,1)
                 endif
                 elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                 elpre(inode) = press(ipoin,1)
              end do

              do inodb = 1,pnodb
                 ipoin = lnodb(inodb,iboun)
                 bovel(1:ndime,inodb) = veloc(1:ndime,ipoin,1)
                 if( kfl_coupl(ID_NASTIN,ID_ALEFOR) /= 0 ) then
                    bovfi(1:ndime,inodb) = velom(1:ndime,ipoin)
                    ! before I was using bvess_nsi(idime,ipoin,1)
                    ! but that introduced problems for cases without mesh movement where an inflow and
                    ! a wall law met because bvess_nsi/=0 there. The solution we decided
                    ! with guillaume is to use velom (and only in the cases with ale coupling).
                    ! Velom is already in the global system - no rotation needed.
                 else
                    bovfi(1:ndime,inodb) = 0.0_rp
                 end if
              end do
              call nsi_elmgap(&
                   pnodb,pmate,lnodb(1,iboun),bocod,bopre,&
                   bovel,bomut,botem)
              !
              ! Element length HLENG
              !
              call elmlen(&
                   ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
                   hnatu(pelty),hleng)
              !
              ! Properties: GPPOR, GPVIS and GPDEN
              !
              gbcar = 0.0_rp
              gbmut = 0.0_rp
              grvis = 0.0_rp
              gbgvi = 0.0_rp
              call ker_proper('DENSI','PGAUB',dummi,iboun,gbden)
              call ker_proper('VISCO','PGAUB',dummi,iboun,gbvis)
              call ker_proper('TURBU','PGAUB',dummi,iboun,gbmut)

              call nsi_turbul(&
                   -1_ip,0_ip,pnodb,pgaub,1_ip,pgaub,kfl_cotur_nsi,elmar(pblty)%shape,  &
                   gbcar,bovel,gbden,gbvis,gbmut,gbgvi,grvis,gbgve,ielem,kfl_kemod_ker)

              gauss_points: do igaub = 1,pgaub

                 do ievat = 1,nevat_nsi
                    do jevat = 1,nevat_nsi
                       wmatr(jevat,ievat) = 0.0_rp
                    end do
                 end do
                 wrhsi(1:nevat_nsi) = 0.0_rp
                 tract(1:3)         = 0.0_rp
                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&    ! Cartesian derivative
                      bocod,baloc,eucta)                                   ! and Jacobian
                 gbsur(igaub) = elmar(pblty)%weigp(igaub)*eucta
                 if( itask == 1_ip ) then
                    do inodb = 1,pnodb
                       elmas(inodb) = elmas(inodb) + gbsur(igaub) * elmar(pblty) % shape(inodb,igaub)
                    end do
                 end if
                 !
                 ! Adjoint problem
                 !
                 call chenor(pnode,baloc,bocod,elcod)
                 if( kfl_adj_prob == 1 ) then
                    call bouder_der(&
                         pnode,pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&
                         bocod,elcod,baloc_der,eucta_der)
                    do inodb = 1, pnodb
                       do idime = 1,ndime
                          gbsur_der(igaub,idime,inodb) = elmar(pblty) % weigp(igaub) * eucta_der(idime,inodb)
                       end do
                    end do
                 end if

#ifdef matiaslma
                 !
                 ! Boundary term coming from by parts integration in mass equation
                 !
                 !call nsi_boumas(pevat,pnodb,solve(ivari_nsi)%ndofn, &
                 !     lboel(1,iboun),elmar(pblty)%shape(1,igaub), &
                 !     gbden(igaub),baloc,wmatr)

                 call nsi_boumas(pevat,pnodb,ndim1, &
                      lboel(1,iboun),elmar(pblty)%shape(1,igaub), &
                      gbden(igaub),baloc,wmatr)

#endif

                 if(  kfl_fixbo_nsi(iboun) == 2 .or. &
                      kfl_fixbo_nsi(iboun) == 6 ) then
                    !
                    ! Pressure: sig.n = - p n
                    !
                    tract(1:ndime) = -bvnat_nsi(1,iboun,1) * baloc(1:ndime,ndime)

                 else if( kfl_fixbo_nsi(iboun) == 3 .or. kfl_fixbo_nsi(iboun) == 13 .or. kfl_fixbo_nsi(iboun) == 18 ) then
                    !
                    ! Wall law: sig.n = - rho* (U*^2) * (u_tan-u_fix_tan)/|u_tan-u_fix_tan|
                    !
                    if( kfl_rough > 0 ) then
                       roughness = 0.0_rp
                       do inodb = 1,pnodb
                          ipoin = lnodb(inodb,iboun)
                          roughness = roughness + rough(ipoin) * elmar(pblty) % shape(inodb,igaub)
                       end do
                    else
                       roughness = rough_dom
                    end if
                    kinen = 0.0_rp
                    if( kfl_ustar == 2 ) then
                       do inodb = 1,pnodb
                          ipoin = lnodb(inodb,iboun)
                          kinen = kinen + untur(1,ipoin,1) * elmar(pblty) % shape(inodb,igaub)
                       end do

                       if( kfl_logva == 1 ) kinen = exp(kinen)

                    end if

                    if ( kfl_waexl_ker == 1_ip ) then !if exchange location for wall law
                       velex(1:ndime) = velel_ker(1:ndime,lexlo_ker(igaub,iboun))
                    else
                       velex(1:ndime) = 0.0_rp
                    end if

                    if ( kfl_wlaav_ker == 1_ip ) then ! Time-averaged velocity for wall law
                       veave(1:ndime) = velav_ker(1:ndime,igaub,iboun)
                    else
                       veave(1:ndime) = 0.0_rp
                    end if

                    !call nsi_bouwal(&
                    !     1_ip,pevat,pnodb,solve(ivari_nsi) % ndofn,iboun,                &
                    !     lboel(:,iboun),elmar(pblty) % shape(1,igaub),bovel,bovfi,tract, &
                    !     gbvis(igaub),gbden(igaub),baloc,ustar,wmatr,roughness,kinen,    &
                    !     velex,veave,igaub,lelbo(iboun))
                    call nsi_bouwal(&
                         1_ip,pevat,pnodb,ndim1,iboun,                &
                         lboel(:,iboun),elmar(pblty) % shape(1,igaub),bovel,bovfi,tract, &
                         gbvis(igaub),gbden(igaub),baloc,ustar,wmatr,roughness,kinen,    &
                         velex,veave,igaub,lelbo(iboun))

                    if( kfl_fixbo_nsi(iboun) == 13 ) then
                       !
                       ! Weak imposition of u.n: assemble ( p n , v )_S on LHS
                       !
                       tract = 0.0_rp
                       wmatr = 0.0_rp
                       wrhsi = 0.0_rp

                       do igaus = 1,pgaus
                          call elmder(&
                               pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&      ! Cartesian derivative
                               elcod,gpcar(1,1,igaus),detjm,xjacm,xjaci)        ! and Jacobian
                       end do
                       call cartbo(&
                            2_ip,lboel(1,iboun),elmar(pblty)%shape(1,igaub),&
                            elmar(pelty)%shaga,gpcar,elmar(pelty)%shape,&
                            shapp,gbcar,pnodb,pnode,pgaus)
                       call nsi_zobi(&
                            pevat,pnodb,pnode,lboel(1,iboun),elmar(pblty)%shape(1,igaub),&
                            shapp,gbcar,gbvis(igaub),baloc(1,ndime),wmatr)
                    end if

                 else if( kfl_fixbo_nsi(iboun) == 22 ) then
                    !
                    ! Two-layer wall modelling: The traction is calculated from an auxiliary RANS simulation
                    !
                    tract = 0.0_rp
                    do inodb = 1,pnodb
                       ipoin = lnodb(inodb,iboun)
                       tract(1:ndime) = tract(1:ndime) &
                            + elmar(pblty)%shape(inodb,igaub) * btrac_nsi(1:ndime,ipoin)
                    end do

                 else if( kfl_fixbo_nsi(iboun) ==  5 .or.&
                      &   kfl_fixbo_nsi(iboun) == 10 ) then
                    !
                    ! Dynamic pressure: sig.n = - (-1/2*rho*u^2) n
                    !
                    tract(1:ndime) = -bvnat_nsi(1,iboun,1) * baloc(1:ndime,ndime)

                 end if

                 if(     kfl_fixbo_nsi(iboun) ==  6 .or.&
                      &  kfl_fixbo_nsi(iboun) == 10 ) then
                    !
                    ! Open boundary: assemble -2*mu*Sym(grad(u).n and possibly pI.n (check nsi_bouopb)
                    !
                    do igaus = 1,pgaus
                       call elmder(&
                            pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&      ! Cartesian derivative
                            elcod,gpcar(1,1,igaus),detjm,xjacm,xjaci)        ! and Jacobian
                    end do
                    call cartbo(&
                         2_ip,lboel(1,iboun),elmar(pblty) % shape(1,igaub),&
                         elmar(pelty) % shaga,gpcar,elmar(pelty) % shape,  &
                         shapp,gbcar,pnodb,pnode,pgaus)
                    gpvis = gbvis(igaub)
                    if( kfl_cotur_nsi == 1 ) then
                       gpvis = gpvis + gbmut(igaub)
                    end if
                    !call nsi_bouopb(&
                    !     lboel(1,iboun),elmar(pblty) % shape(1,igaub),gbcar,&
                    !     baloc(1,ndime),wmatr,pnode,pnodb,solve(ivari_nsi) % ndofn,&
                    !     pevat,gpvis,shapp)
                    call nsi_bouopb(&
                         lboel(1,iboun),elmar(pblty) % shape(1,igaub),gbcar,&
                         baloc(1,ndime),wmatr,pnode,pnodb,ndim1,&
                         pevat,gpvis,shapp)

                 else if( kfl_fixbo_nsi(iboun) == 11 .or. ( kfl_fixbo_nsi(iboun) == 3 .and. kfl_addpr_nsi == 1) ) then
                    !
                    ! Open boundary: assemble pI.n or R*T*rho*I.n - or wall law and correction 'do nothing'
                    !
                    do igaus=1,pgaus
                       call elmder(&
                            pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&      ! Cartesian derivative
                            elcod,gpcar(1,1,igaus),detjm,xjacm,xjaci)        ! and Jacobian
                    end do
                    call cartbo(&
                         2_ip,lboel(1,iboun),elmar(pblty)%shape(1,igaub),&
                         elmar(pelty)%shaga,gpcar,elmar(pelty)%shape,&
                         shapp,gbcar,pnodb,pnode,pgaus)
                    !call nsi_bouopp(&
                    !     lboel(1,iboun),lnodb(1,iboun),elmar(pblty)%shape(1,igaub),&
                    !     baloc(1,ndime),wmatr,pnode,pnodb,solve(ivari_nsi)%ndofn,&
                    !     pevat,shapp)
                    call nsi_bouopp(&
                         lboel(1,iboun),lnodb(1,iboun),elmar(pblty)%shape(1,igaub),&
                         baloc(1,ndime),wmatr,pnode,pnodb,ndim1,&
                         pevat,shapp)

                 else if( kfl_fixbo_nsi(iboun) == 12 ) then
                    !
                    ! Pressure: sig.n = - p n, with p is nodal
                    !
                    gbpre(igaub) = 0.0_rp
                    do inodb = 1,pnodb
                       ipoin = lnodb(inodb,iboun)
                       gbpre(igaub) = gbpre(igaub) + bpess_nsi(1,ipoin,1) * elmar(pblty)%shape(inodb,igaub)
                    end do
                    tract(1:ndime) = - gbpre(igaub) * baloc(1:ndime,ndime)

                 else if( kfl_fixbo_nsi(iboun) == 15 ) then
                    !
                    ! Pressure depends on density
                    ! sig.n = - rho g.r n
                    !
                    dummr=0.0_rp
                    do idime=1,ndime
                       dummr = dummr + bocod(idime,igaub)* gravi_nsi(idime)
                    enddo
                    tract(1:ndime) =  - (gbden(igaub)-prthe(1)/(lowtr_nsi * gasco) )  &
                         * abs(dummr) * grnor_nsi *  baloc(1:ndime,ndime) ! p = (rho- rho_0) g.r !-(lowpr_nsi/(lowtr_nsi * gasco))

                 else if( kfl_fixbo_nsi(iboun) == 17 ) then
                    call runend('NASTIN: Boundary code 17 not coded')
                    !
                    ! Open boundary: assemble 2*mu*eps(grad(u).n-div.u/3)-p I.n
                    !
                    do igaus=1,pgaus
                       call elmder(&
                            pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&      ! Cartesian derivative
                            elcod,gpcar(1,1,igaus),detjm,xjacm,xjaci)        ! and Jacobian
                    end do
                    call cartbo(&
                         2_ip,lboel(1,iboun),elmar(pblty)%shape(1,igaub),&
                         elmar(pelty)%shaga,gpcar,elmar(pelty)%shape,&
                         shapp,gbcar,pnodb,pnode,pgaus)
                    gpvis = gbvis(igaub)
                    if( kfl_cotur_nsi == 1 ) then
                       gpvis = gpvis + gbmut(igaub)
                    end if
                    !call nsi_bouopb(&
                    !     lboel(1,iboun),elmar(pblty)%shape(1,igaub),gbcar,&
                    !     baloc(1,ndime),wmatr,pnode,pnodb,solve(ivari_nsi)%ndofn,&
                    !     pevat,gpvis,shapp)
                    !call nsi_bouopp(&
                    !     lboel(1,iboun),elmar(pblty)%shape(1,igaub),gbcar,&
                    !     baloc(1,ndime),wmatr,pnode,pnodb,solve(ivari_nsi)%ndofn,&
                    !     pevat,shapp)

                    call nsi_bouopb(&
                         lboel(1,iboun),elmar(pblty)%shape(1,igaub),gbcar,&
                         baloc(1,ndime),wmatr,pnode,pnodb,ndim1,&
                         pevat,gpvis,shapp)
                    call nsi_bouopp(&
                         lboel(1,iboun),elmar(pblty)%shape(1,igaub),gbcar,&
                         baloc(1,ndime),wmatr,pnode,pnodb,ndim1,&
                         pevat,shapp)

                 else if( kfl_fixbo_nsi(iboun) == 18 ) then

                    !call nsi_boupen(&
                    !     pevat,pnodb,solve(ivari_nsi)%ndofn,lboel(1,iboun),&
                    !     elmar(pblty)%shape(1,igaub),bovel,bovfi,&
                    !     gbvis(igaub),gbden(igaub),hleng,baloc,wmatr)
                    call nsi_boupen(&
                         pevat,pnodb,ndim1,lboel(1,iboun),&
                         elmar(pblty)%shape(1,igaub),bovel,bovfi,&
                         gbvis(igaub),gbden(igaub),hleng,baloc,wmatr)

                 else if( bemol_nsi > 0.0_rp .and. kfl_fixbo_nsi(iboun) /= 13 .and. kfl_fixbo_nsi(iboun) /= 21 ) then
                    !
                    ! Integrate  ( bemol * rho * (u.n) u . v )_S on LHS
                    !
                    call nsi_boubem(&
                         pevat,pnodb,lboel(1,iboun),elmar(pblty)%shape(1,igaub),bovel,&
                         baloc(1,ndime),gbden(igaub),wmatr)

                 else if( abs(fcons_nsi) > 0.0_rp .and. kfl_fixbo_nsi(iboun) /= 13 .and. kfl_fixbo_nsi(iboun) /= 21 ) then
                    !
                    ! Boundary term  ( fcons * rho * (u.n) u . v )_S on LHS for the
                    ! conservative and skew-symmetric forms of the convective term
                    !
                    ! kfl_fixbo_nsi(iboun) /= 21 is a temporary solution for mass flow control
                    ! so that the boundary term is not applied on the auxiliary set
                    !
                    if(    kfl_stabi_nsi /= NSI_SPLIT_OSS           .and. &
                         & kfl_stabi_nsi /= NSI_GALERKIN            .and. &
                         & kfl_stabi_nsi /= NSI_ALGEBRAIC_SPLIT_OSS       ) then
                       call nsi_bouske(&
                            pevat,pnodb,lboel(1,iboun),elmar(pblty)%shape(1,igaub),bovel,&
                            baloc(1,ndime),gbden(igaub),wmatr)
                    end if
                    
                 else if( kfl_fixbo_nsi(iboun) == 19 ) then
                    !
                    ! Impose traction to the boundary
                    !
                    tract(1)   = bvnat_nsi(1,iboun,1)
                    tract(2:3) = 0.0_rp

                 else if( kfl_fixbo_nsi(iboun) == 20 ) then
                    !
                    ! Stable outflow condition with resistance, Bazilevs + WK0
                    !
                    !                           +-
                    ! p.n = sig.n = -p0 n - C * | u.n ds + rho*(u.n)_ u * beta
                    !                           -+S
                    !           +-
                    !           | u.n if u.n < 0 (inflow)
                    ! (u.n)_ = -|
                    !           | 0   otherwise  (outflow)
                    !           +-
                    !
                    ! p0   = BVNAT_NSI(1,IBOUN,1)
                    ! C    = BVNAT_NSI(2,IBOUN,1)
                    ! beta = BVNAT_NSI(3,IBOUN,1)
                    !


                    call nsi_bousto(pevat, pnodb,wmatr,tract, bvnat_nsi(:,iboun,1),&
                         bovel,lnodb(:,iboun),elmar(pblty)%shape(:,igaub), baloc, &
                         gbden(igaub), lboel(:,iboun),ndim1)


                    !iflow = int(bvnat_nsi(4,iboun,1),ip)
                    !gbvel = 0.0_rp
                    !do inodb = 1,pnodb
                    !   gbvel(1:ndime) = gbvel(1:ndime) + bovel(1:ndime,inodb) * elmar(pblty) % shape(inodb,igaub)
                    !end do
                    !udotn = dot_product(gbvel(1:ndime),baloc(1:ndime,ndime))


                 else if( kfl_fixbo_nsi(iboun) == 21 ) then
                    !
                    ! Stable outflow condition with resistance, Bazilevs + WK1 and WK2 (WK2 TO BE PROGRAMMED)
                    !
                    !
                    ! WK1:
                    !                                        +-
                    ! p_{n+1}.n = sig.n =  p_{n} n - C * dt* | u.n ds + rho*(u.n)_ u * beta - dt * p{n} / (C*R2)
                    !                                        -+S
                    !           +-
                    !           | u.n if u.n < 0 (inflow)
                    ! (u.n)_ = -|
                    !           | 0   otherwise  (outflow)
                    !           +-
                    !
                    ! C    = BVNAT_NSI(1,IBOUN,1)
                    ! R2   = BVNAT_NSI(2,IBOUN,1)
                    ! beta = BVNAT_NSI(3,IBOUN,1)
                    !
                    iflow = int(bvnat_nsi(4,iboun,1),ip)
                    gbvel = 0.0_rp
                    do inodb = 1,pnodb
                       gbvel(1:ndime) = gbvel(1:ndime) + bovel(1:ndime,inodb) * elmar(pblty) % shape(inodb,igaub)
                    end do
                    udotn = dot_product(gbvel(1:ndime),baloc(1:ndime,ndime))

                       !
                       ! Pressure: sig.n = - p n, with p is nodal + bazylev + WK1
                       !
                    gbpre(igaub) = 0.0_rp
                    do inodb = 1,pnodb
                       ipoin = lnodb(inodb,iboun)
                       gbpre(igaub) = gbpre(igaub) + bpess_nsi(1,ipoin,1) * elmar(pblty)%shape(inodb,igaub)
                    end do
                    dtaux = 1.0/dtinv_nsi
                    tract(1:ndime) =   (1_rp - dtaux /bvnat_nsi(1,iboun,1)/bvnat_nsi(2,iboun,1))      &
                         * gbpre(igaub) * baloc(1:ndime,ndime)                                        & ! p (1 + dt/C/R2) n
                         - dtaux * outflow_mass(iflow) * baloc(1:ndime,ndime) / bvnat_nsi(1,iboun,1)  & ! - dt * ( \int u.n ds ) n / C
                         + bvnat_nsi(3,iboun,1) * gbden(igaub) * gbvel(1:ndime) * min(udotn,0.0_rp)     ! + rho*(u.n)_ u
                    

                 else if( kfl_fixbo_nsi(iboun) == 1 .and. kfl_adj_prob == 1 ) then
                    !
                    ! Pressure forces derivatives: F = - p n w.r.t. unknown U and coordinates X
                    !
                    gbpre(igaub) = 0.0_rp
                    do inodb = 1,pnodb
                       gbpre(igaub) = gbpre(igaub) + elmar(pblty)%shape(inodb,igaub) * bopre(inodb)
                    end do
                    ! dF/dX
                    do idime = 1, ndime
                      do jdime = 1, ndime
                        do jnodb = 1,pnodb
                          setfp_derx(idime,jdime,jnodb) = -gbsur(igaub)*gbpre(igaub)*baloc_der(idime,ndime,jdime,jnodb) &
                                                        -gbsur_der(igaub,jdime,jnodb)*gbpre(igaub)*baloc(idime,ndime)
                        enddo
                      enddo
                    enddo
                    ! dF/dU
                    do idime = 1, ndime
                      do jnodb = 1,pnodb
                        setfp_derp(idime,jnodb) = -gbsur(igaub)*elmar(pblty)%shape(jnodb,igaub)*baloc(idime,ndime)
                      enddo
                    enddo

                 end if
                 !
                 ! Exact solution: GPRHS
                 !
                 call nsi_elmexa(                                                    &
                      -1_ip,pnodb,elmar(pblty)%shape(1,igaub),bocod,gbden(igaub),    &
                      gbvis(igaub),gbpor(igaub),gbgvi(1,igaub),cutim,baloc(1,ndime), &
                      tract,dummr,dummr)

                 !call nsi_bouass(&
                 !     pevat,solve(ivari_nsi)%ndofn,pnodb,lboel(1,iboun),&
                 !     elmar(pblty)%shape(1,igaub),&
                 !     gbsur(igaub),tract,wmatr,wrhsi,elmat,elrhs)
                 call nsi_bouass(&
                      pevat,ndim1,pnodb,lboel(1,iboun),&
                      elmar(pblty)%shape(1,igaub),&
                      gbsur(igaub),tract,wmatr,wrhsi,elmat,elrhs)

                 if( kfl_adj_prob == 1 ) then

                   !call nsi_bouass_der(&
                   !     pevat,solve(ivari_nsi)%ndofn,pnodb,pnode,lboel(1,iboun),&
                   !     setfp_derp,setfp_derx,elrhs,eldcost_dx)
                   call nsi_bouass_der(&
                        pevat,ndim1,pnodb,pnode,lboel(1,iboun),&
                        setfp_derp,setfp_derx,elrhs,eldcost_dx)

                 endif

              end do gauss_points
              !
              ! Calculating traction at the boundary nodes
              !
              if ( itask == 1_ip ) then
                 do inodb=1, pnodb
                    inode = lboel(inodb,iboun)
                    ipoin = lnods(inode,ielem)
                    !$OMP ATOMIC
                    massb_nsi(ipoin) = massb_nsi(ipoin) + elmas(inodb)
                    if ( kfl_noslw_ker == 0 ) then  ! for no slip wall I only need massb_nsi
                       do idime=1, ndime
                          !ievat = (inode-1) * solve(ivari_nsi)%ndofn + idime
                          ievat = (inode-1) * ndim1 + idime
                          !$OMP ATOMIC
                          notra_nsi(idime,ipoin) = notra_nsi(idime,ipoin) + elrhs(ievat)
                          do jnode =1, pnode
                             do jdime =1, ndime
                                !jevat = (jnode-1) * solve(ivari_nsi)%ndofn + jdime
                                jevat = (jnode-1) * ndim1 + jdime
                                !$OMP ATOMIC
                                notra_nsi(idime,ipoin) = notra_nsi(idime,ipoin) - &
                                     elmat(ievat,jevat) * elvel(jdime,jnode)
                             end do
                          end do
                       end do
                    end if
                 end do
              end if
              !
              ! Schur complement preconditioner: - ( tau*grad(p).n, q )_S
              !
              if( .false. .and. ( kfl_predi_nsi == 2 .or. kfl_predi_nsi == 3 ) ) then
                 call elmchl(&
                      tragl,hleng,elcod,elvel,chave,chale,pnode,&
                      porde,hnatu(pelty),kfl_advec_nsi,kfl_ellsh_nsi)
                 if( kfl_cotur_nsi == 1 ) then
                    do igaub = 1,pgaub
                       gbvis(igaub) = gbvis(igaub) + gbmut(igaub)
                    end do
                 end if
                 call elmcar(&
                      pnode,pgaus,0_ip,elmar(pelty)%weigp,elmar(pelty)%shape,&
                      elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
                      dummr,ielem)
                 call nsi_bousch(&
                      pnode,pnodb,pgaub,pgaus,lnods(1,ielem),lboel(1,iboun),&
                      gpcar,gbcar,elmar(pelty)%shaga,baloc(1,ndime),gbsur,&
                      gbden,gbvis,gbpor,elmar(pblty)%shape,bovel,chale,elmap)
                 call nsi_assmat(&
                      -1_ip,pnode,pnode,lnods(1,ielem),elmap,dummr,dummr,&
                      dummr,lapla_nsi)
              end if
              !
              ! Prescribe Dirichlet boundary conditions
              !
              if( kfl_matdi_nsi == 0 ) &
                   call nsi_elmdir(&
                   1_ip,1_ip,pnode,pevat,ndim1,lnods(1,ielem),&
                   elmat,elrhs)
              !
              ! Assembly
              !

              !if( NSI_MONOLITHIC ) then
              !   call nsi_assemble_monolithic(&
              !        pnode,pevat,lnods(1,ielem),elauu,elaup,elapp,elapu,&
              !        elrbu,elrbp,amatr,rhsid)
              !else
              !   call nsi_assemble_schur(&
              !        1_ip,pnode,pevat,ielem,lnods(1,ielem),elauu,elaup,elapp,elapu,&
              !        elrbu,elrbp,amatr(poauu_nsi),amatr(poaup_nsi),amatr(poapp_nsi),&
              !        amatr(poapu_nsi),rhsid,rhsid(ndbgs_nsi+1))
              !end if

              ! eldcost_dx ---> dcost_dx_nsi note that itask = 2
              if( kfl_adj_prob == 1 ) then

                pevat1 = pnode*ndime
                call nsi_assrhs(&
                     2_ip,ndim1,pnode,pevat1,ndime,&
                     solve(ivari_nsi)%kfl_algso,lnods(1,ielem),elvel,eldcost_dx,&
                     elmat_aux,dcost_dx_nsi)
              endif

              if( NSI_SCHUR_COMPLEMENT ) then
                 !
                 ! Schur complement
                 !
                 call nsi_assrhs(&
                      1_ip,ndim1,pnode,pevat,ndime,&
                      solve(1)%kfl_algso,lnods(1,ielem),elvel,elrhs,&
                      elmat,rhsid)
                 call nsi_assmat(&
                      1_ip,pnode,pevat,lnods(1,ielem),elmat,amatr(poauu_nsi),&
                      amatr(poaup_nsi),amatr(poapu_nsi),amatr(poapp_nsi))
                 
              else if( NSI_FRACTIONAL_STEP ) then
                 !
                 ! Fractional step
                 !
                 call nsi_assembly_fractional_step_boundary_scalar(&
                      pnode,pevat,ielem,lnods(:,ielem),elvel,elpre,elmat,elrhs,&
                      amatr(poaup_nsi:),amatr(poapu_nsi:),rhsid,&
                      rhsid(ndbgs_nsi+1:))
              else
                 !
                 ! Monolithic
                 !
                 !$OMP CRITICAL (nsi_bouope_lock1)
                 call nsi_assrhs(&
                      1_ip,ndim1,pnode,pevat,ndime,&
                      solve(ivari_nsi)%kfl_algso,lnods(1,ielem),elvel,elrhs,&
                      elmat,rhsid)
                 call assmat(&
                      ndim1,pnode,pevat,solve(ivari_nsi)%nunkn,&
                      solve(ivari_nsi)%kfl_algso,ielem,lnods(1,ielem),elmat,amatr)
                 !$OMP END CRITICAL (nsi_bouope_lock1)
              end if

           end if

        end if

     end if

  end do boundaries

end subroutine nsi_bouope

subroutine nsi_boubem(&
     pevat,pnodb,lboel,gbsha,bovel,baloc,gbden,wmatr)

  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_boubem
  ! NAME
  !    nsi_boubem
  ! DESCRIPTION
  !    Integrate int_\Gamma bemol * rho * (u.n) u . v ds
  !    on the left-hand side
  ! USES
  ! USED BY
  !    nsi_matrix
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_nastin, only     :  bemol_nsi
  implicit none
  integer(ip), intent(in)  :: pevat,pnodb
  integer(ip), intent(in)  :: lboel(pnodb)
  real(rp),    intent(in)  :: gbsha(pnodb)
  real(rp),    intent(in)  :: bovel(ndime,pnodb)
  real(rp),    intent(in)  :: baloc(ndime)
  real(rp),    intent(in)  :: gbden
  real(rp),    intent(out) :: wmatr(pevat,pevat)
  integer(ip)              :: inodb,jnodb,ievat,jevat,idime,ndofn
  real(rp)                 :: gbvel(3),udotn

  ndofn    = ndime + 1
  gbvel(1) = 0.0_rp
  gbvel(2) = 0.0_rp
  gbvel(3) = 0.0_rp
  do inodb = 1,pnodb
    do idime = 1,ndime
       gbvel(idime) = gbvel(idime) + bovel(idime,inodb) * gbsha(inodb)
    end do
  end do

  udotn = 0.0_rp
  do idime = 1,ndime
     udotn = udotn + baloc(idime) * gbvel(idime)
  end do
  udotn = bemol_nsi * gbden * udotn


  do inodb = 1,pnodb
     do jnodb = 1,pnodb
        do idime = 1,ndime
           ievat = (lboel(inodb)-1) * ndofn + idime
           jevat = (lboel(jnodb)-1) * ndofn + idime
           wmatr(jevat,ievat) = wmatr(jevat,ievat) + udotn * gbsha(inodb) * gbsha(jnodb)
        end do
     end do
  end do

end subroutine nsi_boubem


subroutine nsi_bouske(&
     pevat,pnodb,lboel,gbsha,bovel,baloc,gbden,wmatr)

  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_bouske
  ! NAME
  !    nsi_bouske
  ! DESCRIPTION
  !    Integrate int_\Gamma fcons * rho * (u.n) u . v ds
  !    on the left-hand side for the conservative and skew-symmetric
  !    forms of the convective term
  ! USES
  ! USED BY
  !    nsi_matrix
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_nastin, only     :  fcons_nsi
  implicit none
  integer(ip), intent(in)  :: pevat,pnodb
  integer(ip), intent(in)  :: lboel(pnodb)
  real(rp),    intent(in)  :: gbsha(pnodb)
  real(rp),    intent(in)  :: bovel(ndime,pnodb)
  real(rp),    intent(in)  :: baloc(ndime)
  real(rp),    intent(in)  :: gbden
  real(rp),    intent(out) :: wmatr(pevat,pevat)
  integer(ip)              :: inodb,jnodb,ievat,jevat,idime,ndofn
  real(rp)                 :: gbvel(3),udotn

  ndofn    = ndime + 1
  gbvel(1) = 0.0_rp
  gbvel(2) = 0.0_rp
  gbvel(3) = 0.0_rp
  do inodb = 1,pnodb
    do idime = 1,ndime
       gbvel(idime) = gbvel(idime) + bovel(idime,inodb) * gbsha(inodb)
    end do
  end do

  udotn = 0.0_rp
  do idime = 1,ndime
     udotn = udotn + baloc(idime) * gbvel(idime)
  end do
  udotn = abs(fcons_nsi) * gbden * udotn

  do inodb = 1,pnodb
     do jnodb = 1,pnodb
        do idime = 1,ndime
           ievat = (lboel(inodb)-1) * ndofn + idime
           jevat = (lboel(jnodb)-1) * ndofn + idime
           wmatr(jevat,ievat) = wmatr(jevat,ievat) + udotn * gbsha(inodb) * gbsha(jnodb)
        end do
     end do
  end do

end subroutine nsi_bouske

subroutine nsi_zobi(&
     pevat,pnodb,pnode,lboel,gbsha,shapp,gbcar,gbvis,baloc,wmatr)

  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_boubem
  ! NAME
  !    nsi_boubem
  ! DESCRIPTION
  !    Integrate int_\Gamma ( p n - 2*mu*div(u) n ).v ds on LHS
  ! USES
  ! USED BY
  !    nsi_matrix
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_nastin, only       :  fvins_nsi
  implicit none
  integer(ip), intent(in)    :: pevat,pnodb,pnode
  integer(ip), intent(in)    :: lboel(pnodb)
  real(rp),    intent(in)    :: gbsha(pnodb)
  real(rp),    intent(in)    :: shapp(pnode)
  real(rp),    intent(in)    :: gbcar(ndime,pnode)
  real(rp),    intent(in)    :: gbvis
  real(rp),    intent(in)    :: baloc(ndime)
  real(rp),    intent(inout) :: wmatr(pevat,pevat)
  integer(ip)                :: inodb,jnodb,ievat,jevat,idime,ndofn
  integer(ip)                :: jdime,jnode
  real(rp)                   :: xfact

  ndofn = ndime + 1

  if( 1 == 2 ) then

     do inodb = 1,pnodb
        do jnode = 1,pnode
           do idime = 1,ndime
              ievat = ( lboel(inodb) - 1 ) * ndofn + idime
              jevat =   jnode * ndofn
              wmatr(ievat,jevat) = wmatr(ievat,jevat) + gbsha(inodb) * shapp(jnode) * baloc(ndime)
           end do
        end do
     end do

  else
     !
     ! (pn,v)
     !
     do inodb = 1,pnodb
        ievat = ( lboel(inodb) - 1 ) * ndofn
        do idime = 1,ndime
           ievat =  ievat + 1
           do jnodb = 1,pnodb
              jevat = lboel(jnodb) * ndofn
              wmatr(ievat,jevat) = wmatr(ievat,jevat) + gbsha(inodb) * gbsha(jnodb) * baloc(ndime)
           end do
        end do
     end do

     return

     !
     ! Laplacian form:    -mu*( div(u) , v )
     ! Divergence form: -2 mu*( div(u) , v )
     !
     if( fvins_nsi > 0.9_rp ) then
        xfact = 2.0_rp * gbvis
     else
        xfact = 1.0_rp * gbvis
     end if
     do inodb = 1,pnodb
        do idime = 1,ndime
           ievat = ( lboel(inodb) - 1 ) * ndofn + idime
           do jnode = 1,pnode
              do jdime = 1,ndime
                 jevat = ( jnode - 1 ) * ndofn + jdime
                 wmatr(ievat,jevat) = wmatr(ievat,jevat) - xfact * gbcar(jdime,jnode) * baloc(ndime) * gbsha(inodb)
              end do
           end do
        end do
     end do
  end if
! do inodb = 1,pnodb
!        do idime = 1,ndime
!           ievat = ( lboel(inodb) - 1 ) * ndofn + idime
!           do jnode = 1,pnode
!              do jdime = 1,ndime
!                 jevat = ( jnode - 1 ) * ndofn + jdime
!                 wmatr(ievat,jevat) = wmatr(ievat,jevat) - gbvis * gbcar(idime,jnode) * baloc(jdime,ndime) * gbsha(inodb)
!              end do
!           end do
!        end do
!     end do

end subroutine nsi_zobi


subroutine nsi_pipi(npopo,nbvar,kfl_symme,ia,ja,an,wa1)
  !-----------------------------------------------------------------------
  !
  ! Compute the diagonal
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  INOTMASTER,NPOIN_TYPE
  implicit none
  integer(ip), intent(in)    :: npopo,nbvar,kfl_symme
  integer(ip), intent(in)    :: ia(*),ja(*)
  real(rp),    intent(in)    :: an(nbvar,nbvar,*)
  real(rp),    intent(inout) :: wa1(*)
  integer(ip)                :: ii,jj,kk,ll

  if( INOTMASTER ) then
     ii = 24753
     jj = ia(ii)
     ll = -1
     do while (jj< ia(ii+1) .and. ll ==-1)
        if(ja(jj)==ii) then
           ll = jj
        end if
        jj = jj+1
     end do
     if(ll/=-1) then
        jj = (ii-1) * nbvar
        do kk= 1, nbvar
           print*,'a=',ii,kk,an(kk,kk,ll)
        end do
     end if
     ii = 25975
     jj = ia(ii)
     ll = -1
     do while (jj< ia(ii+1) .and. ll ==-1)
        if(ja(jj)==ii) then
           ll = jj
        end if
        jj = jj+1
     end do
     if(ll/=-1) then
        jj = (ii-1) * nbvar
        do kk= 1, nbvar
           print*,'a=',ii,kk,an(kk,kk,ll)
        end do
     end if

  end if

end subroutine nsi_pipi

