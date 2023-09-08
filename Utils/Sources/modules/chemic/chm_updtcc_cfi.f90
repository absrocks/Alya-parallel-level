subroutine chm_updtcc_cfi(dtmin)
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

  real(rp)    :: elcon(mnode,nspec_chm,ADR_chm(1)%ntime)  ! <=> conce
  real(rp)    :: elco2(mnode,nspec_chm)                   ! <=> Previous iteration
  real(rp)    :: elcod(ndime,mnode)                       ! <=> coord
  real(rp)    :: elden(mnode,2)                           ! <=> global density
  real(rp)    :: eltem(mnode)                             ! <=> tempe
  real(rp)    :: elmut(mnode)                             ! Turbulence viscosity
  real(rp)    :: elDik(mnode,nspec_chm)                   ! Species diffusion coefficient
  real(rp)    :: elmas(mnode,nspec_chm)                   ! Mass source terms
  real(rp)    :: elmol(mnode)
  real(rp)    :: elvel(ndime,mnode)
  real(rp)    :: elsen(mnode)                             ! Flame sensor for DTFLES
  real(rp)    :: elrrt(mnode)                             ! {w_c * c} for CFI-LES
  real(rp)    :: elkey(mnode)                             ! Turbulent kinetic energy from TURBUL (RANS)
  real(rp)    :: eleps(mnode)                             ! Dissipation rate from TURBUL (RANS)
  real(rp)    :: gpvol(mgaus)                             ! |J|*w 
  real(rp)    :: gpcon(mgaus,nspec_chm,ncomp_chm)         ! <=> conce
  real(rp)    :: gprea(mgaus)                             ! r
  real(rp)    :: gpvel(ndime,mgaus)                       ! u
  real(rp)    :: gpvec(ndime,mgaus)                       ! uc
  real(rp)    :: gpgve(mgaus)                             ! div.uc
  real(rp)    :: gpadv(ndime,mgaus)                       ! 
  real(rp)    :: gpdif(mgaus)                             ! D_k
  real(rp)    :: gpgrd(ndime,mgaus)                       ! grad(k) = grad(D_k)
  real(rp)    :: gprhs(mgaus)                             ! f (all terms)
  real(rp)    :: gpden(mgaus),gpgde(ndime,mgaus)          ! fake rho for elmadr
  real(rp)    :: gpDik(mgaus,nspec_chm)                   ! Species difussion coefficient
  real(rp)    :: gpgDk(ndime,mgaus,nspec_chm)             ! Gradient of difussion coefficient
  real(rp)    :: gpdiv(mgaus)                             ! Divergence of convection
  real(rp)    :: gpmas(mgaus,nspec_chm)                   ! Mass realease of each reaction
  real(rp)    :: gpmol(mgaus)                             ! Average molar mass
  real(rp)    :: gpgmo(ndime,mgaus)                       ! Average molar mass Gradient
  real(rp)    :: gphmo(mgaus)                             ! Average molar mass Laplacian
  real(rp)    :: gpfar(mgaus)                             ! Fuel / air ratio
  real(rp)    :: gpspe(mgaus)                             ! Tabulated flame speed
  real(rp)    :: gpcar(ndime,mnode,mgaus)                 ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)                 ! dNk/dxidxj
  real(rp)    :: gplap(mgaus,mnode)                       ! Laplacian
  real(rp)    :: gpgac(ndime,mgaus,nspec_chm)             ! Gradient(activity) / activity
  real(rp)    :: gplac(mgaus,nspec_chm)                   ! Laplacian(activity)/activity
  real(rp)    :: gpsgs(mgaus)
  real(rp)    :: gpthi(mgaus)                             ! Tabulated flame thickness
  real(rp)    :: sgsef(mgaus)                             ! Subgrid scale wrinkling factor for DTFLES
  real(rp)    :: gpfac(mgaus)                             ! Dynamic thickening factor F for DTFLES
  real(rp)    :: gphco(mgaus)                             ! heat conductivity
  real(rp)    :: gpsph(mgaus)                             ! specific heat
  real(rp)    :: ellev(mnode)
  real(rp)    :: gplev(mgaus)
  real(rp)    :: gpgle(mgaus,nspec_chm)
  real(rp)    :: gpprd(mgaus,nspec_chm)                   ! Production term of c and f equations in the CFI model
  real(rp)    :: gprrt(mgaus)
  real(rp)    :: gpdis(mgaus)
  real(rp)    :: gptau(mgaus)
  real(rp)    :: gppro(mgaus)
  real(rp)    :: gptur(mgaus)

  real(rp)    :: dummr(mgaus*ndime*mnode)
  real(rp)    :: dummw(mgaus,nspec_chm)
  real(rp)    :: chale(3),chave(3),hleng(3),tragl(9)
  integer(ip) :: inode


  if( INOTMASTER ) then
     do ielem = 1,nelem
        !
        ! Initialization
        !
        gptau = 0.0_rp
        gpdiv = 0.0_rp
        gprhs = 0.0_rp
        gpdif = 0.0_rp
        gprea = 0.0_rp
        gpsgs = 0.0_rp
        gptur = 0.0_rp
        gppro = 0.0_rp
        gpmol = 0.0_rp
        gpgrd = 0.0_rp
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
        !
        call chm_elmgac_cfi(&
          ielem,pnode,lnods(1,ielem),elden,elcod,elcon(1:pnode,:,:),elvel,eltem,elDik,&
          elmas,elkey,eleps,elrrt)
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

        call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty) % shape,gpcar)
        call ker_proper('CONDU','PGAUS',dummi,ielem,gphco,pnode,pgaus,elmar(pelty) % shape,gpcar)
        call ker_proper('SPHEA','PGAUS',dummi,ielem,gpsph,pnode,pgaus,elmar(pelty) % shape,gpcar)
        call ker_proper('TURBU','PGAUS',dummi,ielem,gptur,pnode,pgaus,elmar(pelty) % shape,gpcar)

        !
        ! Send quantities to gauss points
        !
        call chm_elmpre_cfi(&
             pnode,pgaus,elcon(1:pnode,:,:),elden,elvel,elDik,elmas,elkey,eleps,elrrt,elmar(pelty)%shape,gpcar,gplap, &
             gpcon(1:pgaus,:,:),gpvel,gpDik,gpgDk,gpmas,gpmol,gpgmo,gphmo,gpvol,gpdiv,gpdis,gpprd(1:pgaus,:),&
             gprrt,gptur,gpden,gphco,gpsph)
        !
        !
        do iclas = 1,nspec_chm
           call chm_elmprc_cfi(&
                iclas,pgaus,gpcon(1:pgaus,:,:),gpden,gpgDk,gpmas,gphco,gpsph,gptur,gpdis,gpprd(1:pgaus,:), &
                gprrt,gpdif,gpgrd,gprea,gprhs)

           call ADR_critical_time_step(ADR_chm(iclas),gpden,gpvel,gpdif,gprea,dtcri,chale(1),chale(2))

           dtmin= min(dtmin,dtcri(1))
        end do

     end do

  end if
  !
  ! Look for minimum over subdomains
  !
  call pararr('MIN',0_ip,1_ip,dtmin)

end subroutine chm_updtcc_cfi
