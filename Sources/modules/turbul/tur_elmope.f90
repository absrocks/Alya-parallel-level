subroutine tur_elmope(itask)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmop2
  ! NAME
  !   tur_elmop2
  ! DESCRIPTION
  !    Compute elemental matrix and RHS for the Spalart-Allmaras model.
  ! USES
  !    tur_elmgat
  ! USED BY
  !    tur_matrix
  !***
  !-----------------------------------------------------------------------
  use def_elmtyp
  use def_master
  use def_domain
  use def_turbul
  use mod_ker_proper 
  use def_kermod 
  use mod_ADR, only : ADR_element_assembly
  use mod_ADR, only : ADR_bubble_assembly
  use mod_ADR, only : ADR_projections_and_sgs_assembly
  use mod_ADR, only : ADR_add_sgs_or_bubble
  use mod_ADR, only : ELEMENT_ASSEMBLY             ! 1
  use mod_ADR, only : PROJECTIONS_AND_SGS_ASSEMBLY ! 4
  use mod_ADR, only : BUBBLE_ASSEMBLY              ! 5
  use mod_ADR, only : mreac_adr

  implicit none

  integer(ip), intent(in) :: itask
  real(rp)    :: elmat(mnode,mnode)                       ! Element matrices
  real(rp)    :: elrhs(2*mnode)                           ! For properties
  real(rp)    :: elrpr(2*mnode)

  real(rp)    :: eledd(mnode,2)                           ! Gather 
  real(rp)    :: elcod(ndime,mnode)                       ! Coordinates x
  real(rp)    :: elvel(ndime,mnode)                       ! Velocity u
  real(rp)    :: eltem(mnode)                             ! Temperature T
  real(rp)    :: elwal(mnode)                             ! Distance to wall y
  real(rp)    :: eltur(nturb_tur,mnode,3)                 ! Turb. variables
  real(rp)    :: elunk(mnode,2)
  real(rp)    :: elust(mnode)                             ! U*
  real(rp)    :: elrg2(mnode)                             ! (d^2 ui) / ( dxm dxn )
  real(rp)    :: elsqk(mnode)                             ! sqrt(k)
  real(rp)    :: elgrp(ndime,mnode)                       ! grad(phi)
  real(rp)    :: elpro(mnode)                             ! Projection
  real(rp)    :: elprr(mnode)                             ! Projection of reaction (split oss)
  real(rp)    :: elfle(mnode)                             ! Level set
  real(rp)    :: elmsh(ndime,mnode)                       ! u mesh
  real(rp)    :: elprd(mnode)                             ! Smoothed production term

  real(rp)    :: tragl(ndime,ndime),chave(ndime,2)        ! Element length
  real(rp)    :: chale(2),hleng(ndime) 
  real(rp)    :: gpstt(mgaus)                             ! Laplacian

  integer(ip) :: ielem,igaus,ipoin,dummi                  ! Indices and dimensions
  integer(ip) :: pelty,pnode,pgaus,porde,pevat
  integer(ip) :: plapl,pmate,ptopo,kelem

  real(rp)    :: dtcri                                    ! Critical time step  
  real(rp)    :: fddes(mgaus)                             ! DDES blending factor
  real(rp)    :: gddes(mgaus)                             ! DDES blending function
  real(rp)    :: sstf1(mgaus)                             ! SST blending function
  real(rp)    :: sasso(mgaus)                             ! SAS source term
  real(rp)    :: gpvol(mgaus)                             ! |J|*w,|J|
  real(rp)    :: gprea(mgaus)                             ! r
  real(rp)    :: gpvel(ndime,mgaus)                       ! a
  real(rp)    :: gppro(mgaus)                             ! Weighted residual L2-projection
  real(rp)    :: gpprr(mgaus)                             ! Weighted residual L2-projection
  real(rp)    :: gpdiv(mgaus)                             ! Divergence of convection
  real(rp)    :: gpdif(mgaus)                             ! k
  real(rp)    :: gpgrd(ndime,mgaus)                       ! k, grad(k)
  real(rp)    :: gprhs(mgaus)                             ! f
  real(rp)    :: gprec(mgaus)                             ! reactive term
  real(rp)    :: gpcon(mgaus)                             ! convective term
  real(rp)    :: gpres(mgaus)                             ! residual term
  real(rp)    :: gpgrv(ndime,ndime,mgaus)                 ! grad(a)
  real(rp)    :: gpden(mgaus),gpvis(mgaus)                ! rho,mu
  real(rp)    :: gptur(nturb_tur,3,mgaus)                 ! Turb. variables
  real(rp)    :: gpsha(mnode,mgaus)                       ! N
  real(rp)    :: gpder(ndime,mnode,mgaus)                 ! dNk/dsj
  real(rp)    :: gpcar(ndime,mnode,mgaus)                 ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)                 ! dNk/dxidxj
  real(rp)    :: gpcan(mgaus)                             ! canopy flow ( = cd*LAD*|vel|)
  real(rp)    :: sreac(mgaus)                             ! reaction term in perturbation of test function

  real(rp)    :: gpunk(mgaus,ncomp_tur)                  ! Turb. variables
  real(rp)    :: gprea1(mgaus,mreac_adr)                  ! Reaction
  real(rp)    :: gptvi(mgaus)
  !
  ! Initialization
  !
  gpstt     =  1.0_rp
  gpdiv     =  0.0_rp
  gppro     =  0.0_rp
  gpprr     =  0.0_rp
  gpcan     =  0.0_rp
  gpres     =  0.0_rp
  gprec     =  0.0_rp
  gpcon     =  0.0_rp
  gptur     =  0.0_rp
  gpunk     =  0.0_rp
  gprea1    =  0.0_rp
  dtmax_tur = -1.0_rp
  !
  ! Initialize postprocessing values
  !
  if( kfl_ddesm_tur >= 1 .and.  iunkn_tur == 1 ) then
     do ipoin = 1,npoin
        fddes_tur(ipoin) = 0.0_rp
        gddes_tur(ipoin) = 0.0_rp
     end do
  end if
  if( TUR_SST_K_OMEGA ) then
     do ipoin = 1,npoin
        sstf1_tur(ipoin) = 0.0_rp
     end do
     if ( kfl_sasim_tur == 1 ) then
        do ipoin = 1,npoin
           sasso_tur(ipoin) = 0.0_rp
        end do
     end if
  end if
  !
  ! Loop over elements
  !
  elements: do ielem = 1,nelem

     if( lelch(ielem) /= ELHOL ) then
        !
        ! Initialize
        !
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        porde = lorde(pelty)
        plapl = llapl(pelty)
        ptopo = ltopo(pelty)
        pevat = pnode
        !
        ! Check if element is a solid
        !
        pmate = 1
        if( nmate > 1 ) pmate = lmate(ielem) 
        !
        ! Gather operations
        !
        call tur_elmgat(&
             pnode,lnods(1,ielem),eltur,elvel,elcod,elwal,eledd,&
             elust,eltem,elrg2,elsqk,elgrp,elpro,elunk,elfle,elmsh,&
             elprd, elprr)
        !
        ! HLENG and TRAGL at center of gravity
        !
        call elmlen(& 
             ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),&
             hleng)
        !
        ! Compute the characteristic length CHALE
        !
        call elmchl(&
             tragl,hleng,elcod,elvel,chave,chale,pnode,&
             porde,hnatu(pelty),kfl_advec_tur,kfl_ellen_tur)
        !
        ! Local time step
        !
        if( kfl_timco == 2 .and. solve(iunkn_tur) % kfl_algso /= -2 ) then 
           call tur_elmtss(&
                pelty,pnode,ielem,elvel,eledd,dtcri,eltur,elfle,&
                lnods(:,ielem),chale,hleng)
           ! Maximum time step between that given by the global safety factor saflo_nsi, 
           ! and local safet_nsi
           dtinv_tur = min(1.0_rp / (dtcri*safet_tur), 1.0_rp/(dtcri_tur*saflo_tur))
           if( kfl_stead_tur == 1 ) dtinv_tur = 0.0_rp
           if( kfl_timei_tur == 0 ) dtinv_tur = 0.0_rp
           ! Stores maximum time step
           dtmax_tur = max(dtmax_tur, 1.0_rp/dtinv_tur)           
        end if
        !
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
        !
        gphes = 0.0_rp
        call elmca2(&
             pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
             elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
             gpder,gpcar,gphes,ielem)
        !call elmcar(&
        !     pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
        !     elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
        !     gphes,ielem)
        ! 
        !  Properties: GPDEN and GPVIS
        !
        call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,gpsha,gpcar)
        call ker_proper('VISCO','PGAUS',dummi,ielem,gpvis,pnode,pgaus,gpsha,gpcar)
        call ker_proper('POROS','PGAUS',dummi,ielem,gpcan,pnode,pgaus,gpsha,gpcar) ! Porosity for canopy flow
        call ker_proper('TURBU','PGAUS',dummi,ielem,gptvi,pnode,pgaus,gpsha,gpcar)  
        call ker_proper('GRTUR','PGAUS',dummi,ielem,gpgrd,pnode,pgaus,gpsha,gpcar) 
    
        if( itask >= 10 .and. itask <= 20 ) then
           !
           ! Properties
           !
           call asspro(&
                itask,pnode,2_ip,pgaus,lnods(:,ielem),lelch(ielem),gpden,&
                gpvis,gpvol,gpsha,elrpr,rhsid)

        else
           !
           ! Equation coefficients: GPDIF, GPREA, GPRHS, GPGRD and GPVEL(NDIME)
           !
           do igaus = 1,pgaus      
              call tur_elmco2(&
                   pnode,gpvis(igaus),gpden(igaus),gpsha(1,igaus),&
                   gpcar(1,1,igaus),elvel,elmsh,eltur,eledd,elwal,elust,eltem,&
                   elrg2,elsqk,elgrp,elpro,elprd,gpvol(igaus),chale,gpdif(igaus),&
                   gpgrd(1,igaus),gprea(igaus),gpvel(1,igaus),gprhs(igaus),&
                   gpgrv(1,1,igaus),gptur(1,1,igaus),gppro(igaus),hleng,fddes(igaus),&
                   gddes(igaus),sstf1(igaus),sasso(igaus), gpcan(igaus),sreac(igaus), &
                   gpprr(igaus),elprr,pmate, gptvi(igaus))
           end do

           if ( kfl_ddesm_tur >= 1 .and.  iunkn_tur == 1) then    ! DDES model flag
              call tur_pprdes(lnods(:,ielem),pnode,fddes,gddes,pgaus,gpvol,gpsha)
           end if
           if ( TUR_SST_K_OMEGA ) then
              call tur_pprsst(1_ip,lnods(:,ielem),pnode,sstf1,pgaus,gpvol,gpsha)
           end if
           if ( kfl_sasim_tur == 1 ) then
              call tur_pprsst(2_ip,lnods(:,ielem),pnode,sasso,pgaus,gpvol,gpsha)
           end if

           do igaus = 1,mgaus
              gpunk(igaus,1:ncomp_tur) = gptur(iunkn_tur,1:ncomp_tur,igaus)
              gprea1(igaus,1)           = gprea(igaus)
           end do

           if( itask == ELEMENT_ASSEMBLY ) then
              ! 
              ! Assemble equation
              ! 
              call ADR_element_assembly(&
                   ielem,pnode,pgaus,elcod,gpsha,gpcar,elmar(pelty) % deriv,gphes,gpvol,chale,&
                   elmar(pelty) % shape_bub,elmar(pelty) % deriv_bub,ADR_tur(iunkn_tur),&
                   cutim,gpden,gpvel,gpdif,gpgrd,gprea,gprhs,&
                   gpunk,elunk,elmat,elrhs)              
              if( solve(iunkn_tur) % kfl_iffix == 0 ) &
                   call tur_elmdir(&
                   1_ip,pnode,pevat,lnods(:,ielem),elmat,elrhs)
              call assrhs(&
                   solve(1)%ndofn,pnode,lnods(:,ielem),elrhs,rhsid)
              call assmat(&
                   1_ip,pnode,pnode,npoin,solve(iunkn_tur)%kfl_algso,&
                   ielem,lnods(:,ielem),elmat,amatr)         

           else if( itask == PROJECTIONS_AND_SGS_ASSEMBLY ) then
              !
              ! Assemble residual projection
              !
              call ADR_projections_and_sgs_assembly(&
                   ielem,pnode,pgaus,elcod,gpsha,gpcar,gphes,gpvol,chale,ADR_tur(iunkn_tur),&
                   cutim,gpden,gpvel,gpdif,gpgrd,gprea1,&
                   gprhs,gpunk,elunk)

           else if( itask == BUBBLE_ASSEMBLY ) then
              !
              ! Update bubble
              !
              call ADR_bubble_assembly(&
                   ielem,pnode,pgaus,elcod,gpsha,gpcar,elmar(pelty) % deriv,gphes,gpvol,chale,&
                   elmar(pelty) % shape_bub,elmar(pelty) % deriv_bub,ADR_tur(iunkn_tur),&
                   cutim,gpden,gpvel,gpdif,gpgrd,gprea1,gprhs,&
                   gpunk,elunk,elmat,elrhs)

           end if
        end if

     end if

  end do elements

end subroutine tur_elmope
!-----------------------------------------------------------------------
! NOTES
! 
! Governing equation:
! L(T) = rho*u/dt + rho*cp*a.grad(u) - div[k*grad(T)] + r*u = f
! The equation is stabilized using the ASGS model:
!  +-            +-                        +-             
!  | L(u)*v dw - | tau*L'(v)*(L(u)-f) dw = | f*v dw 
! -+            -+                        -+
! 
! where L'(v) = -a.grad(v) - k*Lapl(v) - grad(k).grad(v) + r*v
! Let gpadj = -tau*L'(v)
! The weak form is:
!                                        Gal. ASGS
!                                         |    |
!  +-                                     |    |
!  | ( rho*u/dt + rho*a.grad(u) + r*u )*( v + gpadj ) dw
! -+ 
!     +-                            +-
!   + | -div[k*grad(u)]*gpadj dw +  | k*grad(u).grad(v) dw 
!    -+                            -+
!     +-                              
!   = | ( f + rho*u^(n-1)/dt )*( v + gpadj ) dw 
!    -+                   
!
!***
!-----------------------------------------------------------------------

