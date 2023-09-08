!------------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    nsi_elmope.f90
!> @author  Guillaume Houzeaux
!> @brief   Matrix assembly: boundary contribution
!> @details Elemental operations
!!          - Element matrix calculation
!!          - Scatter in global matrix
!!          - Compute SGS
!> @} 
!------------------------------------------------------------------------
subroutine nsi_elmope(itask)
  use def_kintyp,            only : ip,rp
  use def_master,            only : amatr,rhsid,cutim,kfl_timco,lumma
  use def_master,            only : tesgs,current_zone,vesgs
  use def_kermod,            only : kfl_kemod_ker,kfl_prope,kfl_duatss
  use def_domain,            only : ltype,nnode,nelem
  use def_domain,            only : ngaus,llapl,lorde,ltopo
  use def_domain,            only : ngaus,ndime,lnods,nmate
  use def_domain,            only : lelch,elmar,mnode,mgaus
  use def_domain,            only : lmate,hnatu,ntens
  use def_elmtyp,            only : ELEXT
  use def_nastin,            only : ncomp_nsi,dtmax_nsi
  use def_nastin,            only : kfl_stabi_nsi
  use def_nastin,            only : tamin_nsi,tamax_nsi
  use def_nastin,            only : kfl_sgste_nsi,kfl_advec_nsi
  use def_nastin,            only : dtinv_nsi,safet_nsi  
  use def_nastin,            only : dtcri_nsi,saflo_nsi
  use def_nastin,            only : dtsgs_nsi,kfl_stead_nsi
  use def_nastin,            only : kfl_timei_nsi
  use def_nastin,            only : kfl_cotur_nsi,kfl_grvir_nsi
  use def_nastin,            only : prope_nsi,gradv_nsi
  use def_nastin,            only : lapla_nsi,kfl_sgscp_nsi
  use def_nastin,            only : kfl_force_nsi,kfl_shock_nsi
  use def_nastin,            only : kfl_matdi_nsi,poauu_nsi
  use def_nastin,            only : poaup_nsi,poapp_nsi
  use def_nastin,            only : poapu_nsi,ndbgs_nsi
  use def_nastin,            only : kfl_predi_nsi,kfl_ellen_nsi
  use def_nastin,            only : kfl_ellsh_nsi,lforc_material_nsi
  use def_nastin,            only : xforc_material_nsi
  use def_nastin,            only : resis_nsi,itsta_nsi
  use def_nastin,            only : rmsgs_nsi,resgs_nsi
  use mod_ker_proper,        only : ker_proper
  use mod_nsi_subgrid_scale, only : nsi_subgrid_scale_gather
  use mod_nsi_subgrid_scale, only : nsi_subgrid_scale_residual_and_update
  implicit none

  integer(ip), intent(in) :: itask
  !
  ! Element matrices and vectors (stiffness and preconditioner)
  ! 
  real(rp)    :: elauu(mnode*ndime,mnode*ndime)
  real(rp)    :: elaup(mnode*ndime,mnode)
  real(rp)    :: elapp(mnode,mnode)
  real(rp)    :: elapu(mnode,mnode*ndime)
  real(rp)    :: elrbu(ndime,mnode)
  real(rp)    :: elrbp(mnode)
  real(rp)    :: elrhs(6*mnode)
  real(rp)    :: elmap(mnode,mnode)
  !
  ! Gather 
  !
  real(rp)    :: elvel(ndime,mnode,ncomp_nsi)          ! u
  real(rp)    :: elpre(mnode,ncomp_nsi-1)              ! p
  real(rp)    :: elfle(mnode)                          ! phi
  real(rp)    :: elcod(ndime,mnode)                    ! x
  real(rp)    :: elvep(ndime,mnode)                    ! Pi(momentum)
  real(rp)    :: elprp(mnode)                          ! Pi(div(u))
  real(rp)    :: elgrp(ndime,mnode)                    ! Pi(grad(p))
  real(rp)    :: elmut(mnode)                          ! mut
  real(rp)    :: eltem(mnode,ncomp_nsi)                ! T
  real(rp)    :: elwmean(mnode,ncomp_nsi)              ! W mean
  real(rp)    :: elmsh(ndime,mnode)                    ! u mesh
  real(rp)    :: elnor(ndime,mnode)                    ! Normal to the Level Set interface 
  real(rp)    :: elcur(mnode)                          ! Level Set interface curvature
  real(rp)    :: ellum(mnode)                          ! Lumped mass matrix
  !
  ! Indices and dimensions
  !
  integer(ip) :: ielem,igaus,idime,inode,jdime,ipoin
  integer(ip) :: pnode,pgaus,pevat,kelem,dummi
  integer(ip) :: pelty,plapl,porde,pmate,ptopo
  !
  ! Gauss point values
  !
  real(rp)    :: gpsha(mnode,mgaus)                    ! N
  real(rp)    :: gpder(ndime,mnode,mgaus)              ! dN/dsi                
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dN/dxi
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! d2N/dxidxj
  real(rp)    :: gplap(mnode,mgaus)                    ! Lapl(N)
  real(rp)    :: gpvol(mgaus)                          ! w*|J|, |J|
  real(rp)    :: gpvis(mgaus)                          ! Viscosity 
  real(rp)    :: gpgvi(ndime,mgaus)                    ! Viscosity gradients
  real(rp)    :: gpmut(mgaus)                          ! mut
  real(rp)    :: grvis(ndime,mgaus)                    ! grad(mut)
  real(rp)    :: gpgde(ndime,mgaus)                    ! grad(den)
  real(rp)    :: gppor(mgaus)                          ! Porosity
  real(rp)    :: gpden(mgaus)                          ! Density
  real(rp)    :: gpfle(mgaus)                          ! Level set function
  real(rp)    :: gpst1(mgaus)                          ! tau1
  real(rp)    :: gpst2(mgaus)                          ! tau2
  real(rp)    :: gpsp1(mgaus)                          ! tau1'
  real(rp)    :: gpsp2(mgaus)                          ! tau2'
  real(rp)    :: gptt1(mgaus)                          ! tau1'/tau1
  real(rp)    :: gptt2(mgaus)                          ! tau2'/tau2
  real(rp)    :: gpadv(ndime,mgaus)                    ! u+u'
  real(rp)    :: gprhs(ndime,mgaus)                    ! RHS
  real(rp)    :: gprhs_sgs(ndime,mgaus)                ! RHS due to subscales
  real(rp)    :: gprhc(mgaus)                          ! RHS for the continuity equation (Low Mach)
  real(rp)    :: gprh2(mgaus)                          ! RHS for the residual of continuity equation (Low Mach)
  real(rp)    :: gpsgs(ndime,mgaus,2)                  ! u'
  real(rp)    :: gpsgt(mgaus)                          ! T'
  real(rp)    :: gppre(mgaus,ncomp_nsi-1)              ! p
  real(rp)    :: gpvel(ndime,mgaus,ncomp_nsi-1)        ! u
  real(rp)    :: gpgve(ndime,ndime,mgaus)              ! grad(u)
  real(rp)    :: gpgpr(ndime,mgaus,2)                  ! grad(p)
  real(rp)    :: gptem(mgaus,ncomp_nsi)                ! T
  real(rp)    :: gpsgi(ndime,mgaus)                    ! SGS (work array)
  real(rp)    :: gpvep(ndime,mgaus)                    ! -tau1'*R(u) or tau1*rho*(a.grad)u
  real(rp)    :: gpprp(mgaus)                          ! tau2*div(u)
  real(rp)    :: gpgrp(ndime,mgaus)                    ! tau1'*( grad(p) - rho*f )
  real(rp)    :: gpfli(mgaus)                          ! phy_hyd
  real(rp)    :: gphyd(mgaus)                          ! rho_hyd
  real(rp)    :: gpmsh(ndime,mgaus)                    ! u_mesh
  real(rp)    :: gpnor(ndime,mgaus)                    ! LS normal
  real(rp)    :: gpcur(mgaus)                          ! LS curvature
  !
  ! Element characteristics (to compute tau1 and tau2)
  !
  real(rp)    :: tragl(ndime,ndime)
  real(rp)    :: chave(ndime,2)
  real(rp)    :: chale(2),dtmax
  real(rp)    :: hleng(ndime)
  real(rp)    :: dummr,dtcri
  !
  ! Perturbation and residuals
  !
  real(rp)    :: rmomu(mnode,mgaus)                    ! Residual velocity in momentum
  real(rp)    :: rmom2(ndime,ndime,mnode,mgaus)        ! Residual velocity in momentum
  real(rp)    :: rcont(ndime,mnode,mgaus)              ! Residual velocity in continuity
  real(rp)    :: wgrgr(mnode,mnode,mgaus)              ! grad(Ni).grad(Nj)
  real(rp)    :: wgrvi(mnode,mgaus)                    ! grad(mu).grad(Ni)
  real(rp)    :: agrau(mnode,mgaus)                    ! a.grad(Ni)
  real(rp)    :: p1vec(ndime,ndime,mnode,mgaus)        ! Test funct. velocity in momentum
  real(rp)    :: p1ve2(ndime,ndime,mnode,mgaus)        ! Test funct. velocity in momentum
  real(rp)    :: p2vec(ndime,mnode,mgaus)              ! Test funct. velocity in continuity
  real(rp)    :: p2sca(mnode,mgaus)                    ! Test function pressure in continuity
  !
  ! Event and counters
  !
#ifdef EVENT
  call mpitrace_user_function(1)
#endif
  !
  ! Initialization
  !
  do igaus = 1,mgaus   
     gpden(igaus) = 0.0_rp
     gpvis(igaus) = 0.0_rp
     gppor(igaus) = 0.0_rp
     gprhc(igaus) = 0.0_rp
     gprh2(igaus) = 0.0_rp
     do idime = 1,ndime
        gpgvi(idime,igaus)   = 0.0_rp
        gpgde(idime,igaus)   = 0.0_rp
        gpsgs(idime,igaus,1) = 0.0_rp
        gpsgs(idime,igaus,2) = 0.0_rp
        gpmsh(idime,igaus)   = 0.0_rp 
     end do     
     do inode = 1,mnode
        do idime = 1,ndime
           do jdime = 1,ndime
              rmom2(jdime,idime,inode,igaus) = 0.0_rp
              p1ve2(jdime,idime,inode,igaus) = 0.0_rp
           end do
        end do
     end do
  end do
  dtmax = -1.0_rp
  !
  ! Loop over elements
  !
  elements: do ielem = 1,nelem
     pelty = ltype(ielem)
     !
     ! Element properties and dimensions
     !
    !write(*,*)'elements loop',kelem
    if( pelty > 0 ) then
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        plapl = llapl(pelty) 
        porde = lorde(pelty)
        ptopo = ltopo(pelty)
        pevat = ndime * pnode
        if( kfl_stabi_nsi == 2 ) plapl = 0
        !
        ! Check if element is a solid
        !
        if( nmate > 1 ) then
           pmate = lmate(ielem)
        else
           pmate = 1
        end if
        !
        ! Initializations of the subgrid scales
        !
        call nsi_subgrid_scale_gather(ndime,pgaus,ielem,vesgs,gpsgs)
        if( kfl_sgste_nsi == 1 ) then
           do igaus = 1,pgaus 
              gpsgt(igaus) = tesgs(ielem)%a(1,igaus,1)
           end do
        else
           do igaus = 1,pgaus 
              gpsgt(igaus) = 0.0_rp
           end do
        end if
        !
        ! Gather
        !
        call nsi_elmga3(&
             pnode,lnods(1,ielem),elcod,elpre,elvel,elfle,&
             elvep,elprp,elgrp,elmut,eltem,elmsh,elnor,elcur,elwmean)
        !
        ! HLENG and TRAGL at center of gravity
        !
        call elmlen(&
             ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
             hnatu(pelty),hleng)
        !
        ! Compute the characteristic length: CHALE
        ! 
        call elmchl(&
             tragl,hleng,elcod,elvel,chave,chale,pnode,&
             porde,hnatu(pelty),kfl_advec_nsi,kfl_ellen_nsi)
        !
        ! Local time step: DTINV_NSI
        ! 
        if( kfl_timco == 2 ) then 
           call nsi_elmtss(&
                pelty,pmate,pnode,lnods(1,ielem),ielem,elcod,elvel,&
                gpcar,chale,hleng,dtcri)
           !
           ! Maximum time step between that given by the global safety factor saflo_nsi, 
           ! and local safet_nsi
           !
           dtinv_nsi = min(1.0_rp / (dtcri*safet_nsi), 1.0_rp/(dtcri_nsi*saflo_nsi))
           dtsgs_nsi = min(1.0_rp / (dtcri*safet_nsi), 1.0_rp/(dtcri_nsi*saflo_nsi))
           if( kfl_stead_nsi == 1 ) dtinv_nsi = 0.0_rp
           if( kfl_timei_nsi == 0 ) dtinv_nsi = 0.0_rp
           !
           ! Stores maximum time step
           !
           dtmax = max(dtmax, 1.0_rp/dtinv_nsi)           
        end if
        !
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
        !
        call elmca2(&
             pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
             elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpsha,&
             gpder,gpcar,gphes,ielem)
        !
        ! Properties
        !
        call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,gpsha,gpcar)
        call ker_proper('GRDEN','PGAUS',dummi,ielem,gpgde,pnode,pgaus,gpsha,gpcar)
        call ker_proper('VISCO','PGAUS',dummi,ielem,gpvis,pnode,pgaus,gpsha,gpcar)
        !  BIFLU ignores laminar viscosity gradient - BIFL2 includes it
        call ker_proper('GRVIS','PGAUS',dummi,ielem,gpgvi,pnode,pgaus,gpsha,gpcar)
        call ker_proper('POROS','PGAUS',dummi,ielem,gppor,pnode,pgaus,gpsha,gpcar)
        call ker_proper('TURBU','PGAUS',dummi,ielem,gpmut,pnode,pgaus,gpsha,gpcar)  
        call ker_proper('GRTUR','PGAUS',dummi,ielem,grvis,pnode,pgaus,gpsha,gpcar) 

        call nsi_turbul(&
             itask,1_ip,pnode,pgaus,1_ip,pgaus,kfl_cotur_nsi,&
             gpsha,gpcar,elvel,gpden,gpvis,gpmut,&
             gpgvi,grvis,gpgve, ielem, kfl_kemod_ker)
        
        if (kfl_grvir_nsi==0) then ! zero viscosity gradient
           do igaus = 1,mgaus   
              do idime = 1,ndime
                 gpgvi(idime,igaus)   = 0.0_rp
              end do
           end do
        end if

        if( itask >= 10 .and. itask < 20 ) then
           call asspro(&
                itask,pnode,2_ip,pgaus,lnods(1,ielem),lelch(ielem),gpden,gpvis,&
                gpvol,gpsha,elrhs,prope_nsi)

        else if( itask == 20 ) then
           !
           ! Assemble gradients
           !
           call assgra(&
                pnode,pgaus,lnods(1,ielem),gpden,gpvis,gpvol,gpsha,&
                gpcar,elvel,elrhs,gradv_nsi)

        else if( itask == 6 ) then

           !-------------------------------------------------------------
           !
           ! Assemble pressure equation only
           !
           !-------------------------------------------------------------

           call nsi_elmpri(&
                pnode,pgaus,lnods(1,ielem),gpcar,gpvol,gpden,&
                gpvis,gppor,gpsha,chale,elmap)
           call nsi_assmat(&
                -1_ip,pnode,pnode,lnods(1,ielem),elmap,dummr,dummr,&
                dummr,lapla_nsi)                

        else
           !
           ! Residual and RHS
           !
           call nsi_elmres(                                            &
                pnode,pgaus,plapl,gpsha,gpcar,gphes,gpgvi,gpden,gpvis, &
                gppor,gptem,gpsgs,elvel,elpre,elvep,elprp,elgrp, &
                eltem,elmsh,elcod,elnor,elcur,elwmean,hleng,chale,gpvel,gpgpr, &
                rmomu,rmom2,rcont,gprhs,gprhc,gplap,gpadv,gpvep,gpprp, &
                gpgrp,gphyd,gpmsh,gpgve,gpnor,gpcur,gpfle,ielem, &
                gprh2,gppre,gprhs_sgs,dtinv_nsi,gpgde) 
           !
           ! Compute SGS
           !
           if( kfl_sgscp_nsi == 1 ) then
              call nsi_updsgs(                                            &
                   pgaus,pnode,ndime,ielem,chale,elvel,gpadv,gpvis,gpden, &
                   rmomu,rmom2,gprhs,gpgpr,gpvel,gpcar,gpsp1,gpsgs,gpsgi, &
                   gpgve,gpvep,gpgrp,gpst1,gprhs_sgs,dtsgs_nsi,resis_nsi, &
                   itsta_nsi,rmsgs_nsi,resgs_nsi,gppor)    
              call nsi_subgrid_scale_residual_and_update(                 &
                   ndime,pgaus,ielem,gpsgs,vesgs,resgs_nsi)
             
           end if
           !
           ! Exact solution: GPRHS
           !
           call nsi_elmexa(                                            &
                pgaus,pnode,gpsha,elcod,gpden,gpvis,gppor,gpgvi,cutim, &
                dummr,gprhs,gprhc,gprh2)
           !
           ! External force: GPRHS
           !
           if( kfl_force_nsi == 1 ) then
              call nsi_elmexf(                                  &
                   ndime,pgaus,lforc_material_nsi(pmate),gpden, &
                   xforc_material_nsi(1,pmate),gprhs, gpvel)
           end if

           if( itask == 4 ) then

              !-------------------------------------------------------------
              !
              ! Compute Subgrid scale
              !
              !-------------------------------------------------------------

              call nsi_updsgs(                                            &
                   pgaus,pnode,ndime,ielem,chale,elvel,gpadv,gpvis,gpden, &
                   rmomu,rmom2,gprhs,gpgpr,gpvel,gpcar,gpsp1,gpsgs,gpsgi, &
                   gpgve,gpvep,gpgrp,gpst1,gprhs_sgs,dtsgs_nsi,resis_nsi, &
                   itsta_nsi,rmsgs_nsi,resgs_nsi,gppor)
              call nsi_subgrid_scale_residual_and_update(                 &
                   ndime,pgaus,ielem,gpsgs,vesgs,resgs_nsi)
              call nsi_elmsgs(                                            &
                   pgaus,pnode,chale,hleng,gpadv,gpvis,gpden,gpcar,gpst1, &
                   gpst2,gpsp1,gpsp2,gptt1,gptt2,rmomu,gppor,dtsgs_nsi,   &
                   tamin_nsi,tamax_nsi, elvel,gprhs,gpgpr,rmom2,gpgve)
              call nsi_elmort(                                            &
                   ielem,pgaus,pnode,ndime,elvel,elpre,rmomu,rmom2,gprhs, &
                   gpgpr,gpsha,gpvol,gpden,gpadv,gpcar,gpsp1,gpsp2,gpst1  )

           else if( itask == 5 ) then
              !
              ! Limiter
              !
              call nsi_elmsgs(                                               &
                   pgaus,pnode,chale,hleng,gpadv,gpvis,gpden,gpcar,gpst1,    &
                   gpst2,gpsp1,gpsp2,gptt1,gptt2,rmomu, gppor,dtsgs_nsi,     &
                   tamin_nsi,tamax_nsi, elvel,gprhs,gpgpr,rmom2,gpgve)
              call nsi_asslim(                                               &
                   pnode,pgaus,lnods(1,ielem),gpden,gpsp1,gpadv,gpvep,gpvol, &
                   elvel,gpsha,gpcar,wgrvi,elrhs,rhsid)

           else if( itask /= 4 ) then

              !-------------------------------------------------------------
              !
              ! Assemble equations
              !
              !-------------------------------------------------------------
              !
              ! Stabilization parameters
              !
              call nsi_elmsgs(                                              &
                   pgaus,pnode,chale,hleng,gpadv,gpvis,gpden,gpcar,gpst1,   &
                   gpst2,gpsp1,gpsp2,gptt1,gptt2,rmomu, gppor,dtsgs_nsi,    &
                   tamin_nsi,tamax_nsi, elvel,gprhs,gpgpr,rmom2,gpgve)
              !
              ! Assembly
              !
              if( kfl_stabi_nsi == 2 ) then
                 call nsi_elmma4(                                            &
                      pnode,pgaus,pevat,gpden,gpvis,gppor,gpsp1,gpsp2,gpvol, &
                      gpsha,gpcar,gpadv,gpvep,gpprp,gpgrp,gprhs,gpvel,gpsgs, &
                      wgrgr,agrau,elvel,elauu,elaup,elapp,elapu,elrbu,elrbp, &
                      dtinv_nsi,dtsgs_nsi)
              else
                 !if( pgaus == 1 .and. pelty == TET04 ) then
                 !   call nsi_elmp13(                                            &
                 !        pnode,gpden,gpvis,gppor,gpgvi,gpsp1,gptt1,gpsp2,gptt2, &
                 !        gpvol,gpsha,gpcar,gpadv,gpvep,gprhs,rmomu,rcont,p1vec, &
                 !        p2vec,p2sca,wgrgr,wgrvi,elauu,elaup,elapp,elapu,elrbu, &
                 !        elrbp,rmom2,p1ve2,gpst1,gpsgs)
                 !else
                 call nsi_elmmat(                                      &
                      pnode,pgaus,pevat,gpden,gpvis,gppor,gpgvi,gpsp1, &
                      gptt1,gpsp2,gptt2,gpvol,gpsha,gpcar,gplap,gphes, &
                      gpadv,gpvep,gpprp,gpgrp,gprhs,gprhc,rmomu,rcont, &
                      p1vec,p2vec,p2sca,wgrgr,wgrvi,elauu,elaup,elapp, &
                      elapu,elrbu,elrbp,rmom2,p1ve2,gpst1,gpsgs,gpgve, &
                      gprh2,gppre,gprhs_sgs,elvel,ellum,dtinv_nsi)
                 !end if
              end if
              !
              ! Shock capturing
              !
              if( kfl_shock_nsi /= 0 ) then 
                 call nsi_elmshc(                                      &
                      pnode,pgaus,ptopo,pevat,ndime,gpden,gpvel,gprhs, &
                      gpsp1,gpsgs,gpvol,elvel,gpcar,chale,rmomu,rmom2, &
                      elauu)
              end if
              !
              ! Extension elements (Dodeme) 
              !
              if( lelch(ielem) == ELEXT ) then
                 call nsi_elmext(&
                      1_ip,pnode,elauu,elaup,elapu,elapp,elmap,elrbu,elrbp)
              end if

              !if( lelch(ielem) == 32 .and. kfl_colev_nsi /= 0 ) then
              !   call nsi_caca(&
              !        pnode,lnods(1,ielem),elauu,elaup,elapu,elapp,elmap,elrbu,elrbp)
              !end if
              !
              ! Prescribe Dirichlet boundary conditions
              !
              if( kfl_matdi_nsi == 0 ) &
                   call nsi_elmdi3(&
                   pnode,pevat,lnods(1,ielem),&
                   elauu,elaup,elapp,elapu,elrbu,elrbp)
              !
              ! Assembly: AMATR and RHSID
              !
              call nsi_assemble_schur(&            
                   1_ip,pnode,pevat,lnods(1,ielem),elauu,elaup,elapp,elapu,&
                   elrbu,elrbp,amatr(poauu_nsi),amatr(poaup_nsi),amatr(poapp_nsi),&
                   amatr(poapu_nsi),rhsid,rhsid(ndbgs_nsi+1))

              if (kfl_duatss==1) then ! dual time step preconditioner
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    lumma(ipoin) = lumma(ipoin) + ellum(inode)
                 end do
              end if
              !
              ! Schur complement preconditioner
              !
              if( kfl_predi_nsi == 2 .or. kfl_predi_nsi == 3 .or. kfl_predi_nsi == 4 ) then

                 call elmchl(&
                      tragl,hleng,elcod,elvel,chave,chale,pnode,&
                      porde,hnatu(pelty),kfl_advec_nsi,kfl_ellsh_nsi)
                 call nsi_elmsch(&
                      pnode,pgaus,lnods(1,ielem),gpcar,gpvol,gpden,&
                      gpvis,gppor,gpsha,elvel,chale,gpsp1,&
                      elmap,dtinv_nsi)
                 !
                 ! Extension elements (Dodeme)
                 !
                 if( lelch(ielem) == ELEXT ) then
                    call nsi_elmext(&
                         2_ip,pnode,elauu,elaup,elapu,elapp,elmap,elrbu,elrbp)
                 end if
                 call nsi_assmat(&
                      -1_ip,pnode,pnode,lnods(1,ielem),elmap,dummr,dummr,&
                      dummr,lapla_nsi)
              end if
           end if
        end if

     end if
  end do elements
  !
  ! Global variables
  ! 
  dtmax_nsi = dtmax
  !
  ! Event and counters
  !
#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine nsi_elmope

subroutine nsi_caca(&
     pnode,lnods,elauu,elaup,elapu,elapp,elmap,elrbu,elrbp)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmext
  ! NAME 
  !    nsi_elmext
  ! DESCRIPTION
  !    Modify element matrix when extension elements are used
  !    Only equation of the first node should be assembled as
  !    it corresponds to its extension test function
  !
  ! USES
  ! USED BY
  !    nsi_elmext
  !***
  !----------------------------------------------------------------------- 
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime
  use def_master, only       :  fleve
  use def_nastin, only       :  prout_nsi
  implicit none
  integer(ip), intent(in)    :: pnode
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(out)   :: elauu(pnode*ndime,pnode*ndime)
  real(rp),    intent(out)   :: elaup(pnode*ndime,pnode)
  real(rp),    intent(out)   :: elapu(pnode,pnode*ndime)
  real(rp),    intent(out)   :: elapp(pnode,pnode)
  real(rp),    intent(in)    :: elmap(pnode,pnode)
  real(rp),    intent(out)   :: elrbu(pnode*ndime)
  real(rp),    intent(out)   :: elrbp(pnode)
  integer(ip)                :: inode,jnode,idofn,jdofn
  integer(ip)                :: idime,jdime,ipoin,jpoin
  !
  ! Auu, Aup, Apu, App, bu, bp
  !
  do inode = 1,pnode
     ipoin = lnods(inode)
     if( prout_nsi(ipoin) == 1 ) then
        idofn = (inode-1)*ndime
        do idime = 1,ndime
           idofn = idofn + 1
           do jdofn = 1,ndime*pnode
              elauu(idofn,jdofn) = 0.0_rp
           end do
           do jnode = 1,pnode
              elaup(idofn,jnode) = 0.0_rp
           end do
           elrbu(idofn) = 0.0_rp
        end do
        do jdofn = 1,ndime*pnode
           elapu(inode,jdofn) = 0.0_rp
        end do
        do jnode = 1,pnode
           elapp(inode,jnode) = 0.0_rp
        end do
        elrbp(inode) = 0.0_rp
     end if
  end do
  
  return
  do inode = 1,pnode
     ipoin = lnods(inode)
     idofn = (inode-1)*ndime
     do idime = 1,ndime
        idofn = idofn + 1
        do jnode = 1,pnode
           jpoin = lnods(jnode)
           if( prout_nsi(jpoin) == 1 ) then
              jdofn = (jnode-1)*ndime
              do jdime = 1,ndime
                 jdofn = jdofn + 1
                 elauu(idofn,jdofn) = 0.0_rp
              end do
              elaup(idofn,jnode) = 0.0_rp
           end if
        end do
     end do
  end do

  do inode = 1,pnode
     ipoin = lnods(inode)
     idofn = (inode-1)*ndime
     do idime = 1,ndime
        idofn = idofn + 1
        do jnode = 1,pnode           
           jpoin = lnods(jnode)
           if( prout_nsi(jpoin) == 1 ) then
              jdofn = (jnode-1)*ndime
              do jdime = 1,ndime
                 jdofn = jdofn + 1
                 elapu(inode,jdofn) = 0.0_rp
              end do
              elapp(inode,jnode) = 0.0_rp
           end if
        end do
     end do
  end do

end subroutine nsi_caca
