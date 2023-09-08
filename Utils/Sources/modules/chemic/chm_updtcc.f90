subroutine chm_updtcc(dtmin)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_updtsc
  ! NAME 
  !    chm_updtsc
  ! DESCRIPTION
  !    This routine computes the critical time step size
  ! USED BY
  !    chm_updtss
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_ker_proper 
  use def_kermod
  use mod_ADR,            only : ADR_critical_time_step
  use mod_ADR,            only : FROM_CRITICAL
  implicit none 
  real(rp),   intent(inout) :: dtmin
  real(rp)                  :: dtcri(2)
  integer(ip) :: ielem,igaus,idime,iclas               ! Indices and dimensions
  integer(ip) :: pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo,dummi

  real(rp)    :: elcon(mnode,nspec_chm,ncomp_chm)      ! <=> conce
  real(rp)    :: elco2(mnode,nspec_chm)                ! <=> Previous iteration
  real(rp)    :: elcod(ndime,mnode)                    ! <=> coord
  real(rp)    :: elden(mnode,2)                        ! <=> global density
  real(rp)    :: eltem(mnode)                          ! <=> tempe
  real(rp)    :: elmut(mnode)                          ! Turbulence viscosity
  real(rp)    :: elDik(mnode,nspec_chm)                ! Species diffusion coefficient
  real(rp)    :: elmas(mnode,nspec_chm)                ! Mass source terms
  real(rp)    :: elmol(mnode)
  real(rp)    :: elvel(ndime,mnode)
  real(rp)    :: elsen(mnode)                          ! Flame sensor for DTFLES
  real(rp)    :: gpvol(mgaus)                          ! |J|*w 
  real(rp)    :: gpcon(mgaus,nspec_chm,ncomp_chm)      ! <=> conce
  real(rp)    :: gprea(mgaus)                          ! r
  real(rp)    :: gpvel(ndime,mgaus)                    ! u
  real(rp)    :: gpvec(ndime,mgaus)                    ! uc
  real(rp)    :: gpgve(mgaus)                          ! div.uc
  real(rp)    :: gpadv(ndime,mgaus)                    ! 
  real(rp)    :: gpdif(mgaus)                          ! D_k
  real(rp)    :: gpgrd(ndime,mgaus)                    ! grad(k) = grad(D_k)
  real(rp)    :: gprhs(mgaus)                          ! f (all terms)
  real(rp)    :: gpden(mgaus),gpgde(ndime,mgaus)       ! fake rho for elmadr
  real(rp)    :: gpDik(mgaus,nspec_chm)                ! Species difussion coefficient
  real(rp)    :: gpgDk(ndime,mgaus,nspec_chm)          ! Gradient of difussion coefficient
  real(rp)    :: gpdiv(mgaus)                          ! Divergence of convection
  real(rp)    :: gpmas(mgaus,nspec_chm)                ! Mass realease of each reaction
  real(rp)    :: gpmol(mgaus)                          ! Average molar mass
  real(rp)    :: gpgmo(ndime,mgaus)                    ! Average molar mass Gradient
  real(rp)    :: gphmo(mgaus)                          ! Average molar mass Laplacian
  real(rp)    :: gpfar(mgaus)                          ! Fuel / air ratio
  real(rp)    :: gpspe(mgaus)                          ! Tabulated flame speed
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj
  real(rp)    :: gplap(mgaus,mnode)                    ! Laplacian
  real(rp)    :: gpgac(ndime,mgaus,nspec_chm)          ! Gradient(activity) / activity
  real(rp)    :: gplac(mgaus,nspec_chm)                ! Laplacian(activity)/activity
  real(rp)    :: gpsgs(mgaus)
  real(rp)    :: gpthi(mgaus)                          ! Tabulated flame thickness
  real(rp)    :: sgsef(mgaus)                          ! Subgrid scale wrinkling factor for DTFLES
  real(rp)    :: gpfac(mgaus)                          ! Dynamic thickening factor F for DTFLES
  real(rp)    :: gphco(mgaus)                          ! heat conductivity
  real(rp)    :: gpsph(mgaus)                          ! specific heat
  real(rp)    :: ellev(mnode)
  real(rp)    :: gplev(mgaus)
  real(rp)    :: gpgle(mgaus,nspec_chm)

  real(rp)    :: dummr(mgaus*ndime*mnode)
  real(rp)    :: dummw(mgaus,nspec_chm)
  real(rp)    :: chale(3),chave(3),hleng(3),tragl(9)
  integer(ip) :: inode


  if( INOTMASTER ) then
     do ielem = 1,nelem
        !
        ! Initialization
        !
        gpden = 0.0_rp
        gpdiv = 0.0_rp
        gprhs = 0.0_rp
        gpdif = 0.0_rp
        gprea = 0.0_rp
        gpsgs = 0.0_rp
        gpmol = 0.0_rp
        gpadv = 0.0_rp
        gpvel = 0.0_rp
        gpgrd = 0.0_rp
        gphco = 0.0_rp
        gpsph = 0.0_rp
        !
        ! Element dimensions
        !
        pelty = ltype(ielem)
        pnode = nnode(pelty)
        pgaus = ngaus(pelty)
        plapl = llapl(pelty) 
        porde = lorde(pelty)
        ptopo = ltopo(pelty)
        !
        ! Gather all
        !
        call ker_proper('DENSI','PNODE',dummi,ielem,elden(:,1),pnode)
        call ker_proper('TURBU','PNODE',dummi,ielem,elmut,pnode)
        !
 
        call chm_elmgac(&
             2_ip,ielem,pnode,lnods(1,ielem),elden,elcod,elcon,elco2,elvel,eltem,elDik,&
             elmas,elmol,elmut,dummr,dummr,dummr,elsen)

        !
        ! CHALE, HLENG and TRAGL 
        !
        if( kfl_taust_chm /= 0 .or. kfl_shock_chm /= 0 ) then
           call elmlen(&
                ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),&
                hleng)
           call elmchl(&
                tragl,hleng,elcod,dummr,chave,chale,pnode,&
                porde,hnatu(pelty),kfl_advec_chm,kfl_ellen_chm)
        else
           plapl = 0
        end if
        !
        ! Cartesian derivatives, Hessian matrix and volume: GPCAR, GPHES, PGVOL
        !
        call elmcar(&
             pnode,pgaus,plapl,elmar(pelty)%weigp,elmar(pelty)%shape,&
             elmar(pelty)%deriv,elmar(pelty)%heslo,elcod,gpvol,gpcar,&
             gphes,ielem)
        !
        ! Compute laplacian
        !
        if (plapl /= 0 ) then
           do igaus = 1,pgaus
              do inode = 1,pnode
                 gplap(igaus,inode)=0.0_rp
                 do idime=1,ndime
                    gplap(igaus,inode) = gplap(igaus,inode) + gphes(idime,inode,igaus)
                 enddo
              enddo
           enddo
        else
           do igaus = 1,pgaus
              do inode = 1,pnode
                 gplap(igaus,inode)=0.0_rp
              enddo
           enddo
        endif
        !
        ! Send stuff to gauss points
        !
        call chm_elmpre(&
             pnode,pgaus,elcon,elden,elvel,elDik,elmas,elmut,dummr,dummr,dummr,elsen,elmar(pelty)%shape,gpcar,gplap, &
             gpcon,gpvel,gpDik,gpgDk,gpmas,elmol,gpmol,gpgmo,gphmo,gpfar,gpvol,gpdiv,dummr,dummw,gpspe,&
             gpthi,sgsef,gpfac,dummr,dummr,dummr) 

        !
        ! Density from kermod
        !
        call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty) % shape,gpcar)
        call ker_proper('GRDEN','PGAUS',dummi,ielem,gpgde,pnode,pgaus,elmar(pelty) % shape,gpcar)
        call ker_proper('CONDU','PGAUS',dummi,ielem,gphco,pnode,pgaus,elmar(pelty) % shape,gpcar)
        call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty) % shape,gpcar)
        !
        ! Next check coupling with levels and if sucessful then complete the activities
        !
        if (kfl_activ_chm == 1) then
           if( kfl_coupl(ID_CHEMIC,ID_LEVELS) == 1 ) then
              call ker_proper('DUMMY','PNODE',dummi,ielem,ellev,pnode,pgaus,elmar(pelty) % shape,gpcar)
              call ker_proper('DUMMY','PGAUS',dummi,ielem,gplev,pnode,pgaus,elmar(pelty) % shape,gpcar)
              call ker_proper('GRDUM','PGAUS',dummi,ielem,gpgle,pnode,pgaus,elmar(pelty) % shape,gpcar)
              call chm_levels(pnode,pgaus,gplap,ellev,gplev,gpgle,gpgac,gplac)
           else
              call chm_activi( &
                   pnode,pgaus,gpcar,gplap,eltem,elcon,gpgac,gplac)
           endif
        endif
        !
        ! Correction velocity
        !
        call chm_elmvel(&
             pnode,pgaus,elco2,gpDik,gpgDk,elmar(pelty)%shape,gpcar,gplap,&
             gpvec,gpgve,gpmol,gpgmo,gphmo)

        if (kfl_activ_chm == 1) call chm_levvel(&  ! Correct with activities
             pnode,pgaus,elco2,gpDik,gpgDk,elmar(pelty)%shape,gpcar,gplap,&
             gpvec,gpgve,gpmol,gpgmo,gphmo,gpgac,gplac)
        !
        !
        do iclas = 1,nspec_chm
           call chm_elmprc(&
                iclas,pgaus,gpcon(:,:,1),gpden,gpgde,gpDik,gpgDk,gpmas,gpvel,gpvec,gpgve,&
                gpadv,gpdif,gpgrd,gprea,gprhs,gpmol,gpgmo,gphmo,gpgac,gplac)

           if (kfl_model_chm == 5) gpdif(1:pgaus) = gphco(1:pgaus) / gpsph(1:pgaus)

           call ADR_critical_time_step(ADR_chm(iclas),gpden,gpvel,gpdif,gprea,dtcri,chale(1),chale(2))
           !!DMM-ADR  call chm_adr_critical_time(pgaus,chale,gpadv,gpdif,gprea,gpmas(:,iclas),gpden,gpcon(:,iclas),dtcri)

           dtmin= min(dtmin,dtcri(1))
        end do

     end do


  end if
  !
  ! Look for minimum over subdomains
  !
  call pararr('MIN',0_ip,1_ip,dtmin)

end subroutine chm_updtcc





subroutine  chm_adr_critical_time(pgaus,chale,gp_advection,gp_diffusion,gp_reaction,gp_sourceterm,gp_density,gp_concentration,crit_time)
  use def_kintyp, only : ip,rp
  use def_domain, only : ndime
  use def_chemic, only : staco_chm,kfl_stagg_chm
  use mod_tauadr, only :  tauadr
  implicit none
  integer(ip),intent(in) :: pgaus
  real(rp),intent(in)    :: chale(3)
  real(rp),intent(in)    :: gp_advection(ndime,pgaus)
  real(rp),intent(in)    :: gp_diffusion(pgaus)
  real(rp),intent(in)    :: gp_reaction(pgaus)
  real(rp),intent(in)    :: gp_sourceterm(pgaus)
  real(rp),intent(in)    :: gp_density(pgaus)
  real(rp),intent(in)    :: gp_concentration(pgaus)
  real(rp),intent(out)   :: crit_time
  real(rp) :: adv, dif,rea,chtim,advdime
  integer(ip) :: idime, igaus
  !
  ! We keep only the maximum value among the Gauss points
  !
  adv=0.0
  dif=0.0
  rea=0.0
  do igaus=1,pgaus
     advdime = 0.0_rp
     do idime = 1,ndime
        advdime = advdime+gp_advection(idime,igaus)**2
     enddo
     adv = max(adv,sqrt(advdime))
     dif = max(dif,abs(gp_diffusion(igaus)))
     rea = max(rea,abs(gp_reaction(igaus)))
  enddo
  ! ADR critical time
  call tauadr(&
       1_ip,staco_chm,adv,dif,rea,&
       chale(1),chale(2),crit_time)
  !
  ! Reaction critical time
  if (kfl_stagg_chm .eq. 0) then 
     call chm_timchm(pgaus,gp_sourceterm,gp_density,gp_concentration,chtim)
  else
     chtim = 1.e6
  endif
  
  crit_time = min(crit_time,chtim)

end subroutine chm_adr_critical_time
