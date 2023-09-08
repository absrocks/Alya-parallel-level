subroutine chm_elmcfi(order)
  !------------------------------------------------------------------------
  ! NAME 
  !    chm_elmcfi
  ! DESCRIPTION
  !    Elemental operations for the CFI combustion model
  ! USES
  ! USED BY
  !    chm_matrix
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_ker_proper 
  use def_kermod
  use mod_matrix
  use def_solver
 
  use mod_ADR,    only : ADR_element_assembly
  use mod_ADR,    only : ADR_bubble_assembly
  use mod_ADR,    only : ADR_projections_and_sgs_assembly
  use mod_ADR,    only : ADR_add_sgs_or_bubble
  use mod_ADR,    only : ELEMENT_ASSEMBLY                 ! 1
  use mod_ADR,    only : PROJECTIONS_AND_SGS_ASSEMBLY     ! 4
  use mod_ADR,    only : BUBBLE_ASSEMBLY                  ! 5
  use mod_ADR,    only : mreac_adr

  implicit none

  integer(ip), intent(in) :: order                        ! =1 defaul or =2 compute SGS only

  real(rp)    :: elmat(mnode,mnode) 
  real(rp)    :: elrhs(mnode)
  integer(ip) :: ielem,igaus,idime,iclas                  ! Indices and dimensions
  integer(ip) :: izmat,izrhs,pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo
  integer(ip) :: dummi
  integer(ip) :: inode
  integer(ip) :: iboun,pnodb,pblty,inodb,ipoin

  real(rp)    :: elcon(mnode,nspec_chm,ADR_chm(1) % ntime)! <=> conce
  real(rp)    :: elco1(mnode,ncomp_chm)                   ! <=> conce: work array
  real(rp)    :: elcod(ndime,mnode)                       ! <=> coord
  real(rp)    :: elden(mnode,2)                           ! <=> global density
  real(rp)    :: eltem(mnode)                             ! <=> tempe
  real(rp)    :: elrrt(mnode)                             ! {w_c * c} for CFI-LES
  real(rp)    :: elkey(mnode)                             ! Turbulent kinetic energy from TURBUL (RANS)
  real(rp)    :: eleps(mnode)                             ! Dissipation rate from TURBUL (RANS)
  real(rp)    :: elDik(mnode,nspec_chm)                   ! Species diffusion coefficient
  real(rp)    :: elmas(mnode,nspec_chm)                   ! Mass source terms
  real(rp)    :: elvel(ndime,mnode)
  real(rp)    :: gpvol(mgaus)                             ! |J|*w 
  real(rp)    :: gphco(mgaus)                             ! heat conductivity
  real(rp)    :: gpsph(mgaus)                             ! specific heat 
  real(rp)    :: gpcon(mgaus,nspec_chm,ncomp_chm)         ! <=> conce
  real(rp)    :: gprea(mgaus)                             ! r
  real(rp)    :: gpvel(ndime,mgaus)                       ! u
  real(rp)    :: gpvec(ndime,mgaus)                       ! uc
  real(rp)    :: gpgve(mgaus)                             ! div.uc
  real(rp)    :: gpdif(mgaus)                             ! D_k
  real(rp)    :: gpgrd(ndime,mgaus)                       ! grad(k) = grad(D_k)
  real(rp)    :: gprhs(mgaus)                             ! f (all terms)
  real(rp)    :: gpden(mgaus)                             ! fake rho for elmadr
  real(rp)    :: gpDik(mgaus,nspec_chm)                   ! Species difussion coefficient
  real(rp)    :: gptur(mgaus)                             ! turbulent viscosity
  real(rp)    :: gpgDk(ndime,mgaus,nspec_chm)             ! Gradient of difussion coefficient
  real(rp)    :: gppro(mgaus)                             ! Weighted residual L2-projection
  real(rp)    :: gpdiv(mgaus)                             ! Divergence of convection
  real(rp)    :: gpmas(mgaus,nspec_chm)                   ! Mass realease of each reaction
  real(rp)    :: gpmol(mgaus)                             ! Average molar mass
  real(rp)    :: gpgmo(ndime,mgaus)                       ! Average molar mass Gradient
  real(rp)    :: gphmo(mgaus)                             ! Average molar mass Laplacian
  real(rp)    :: gpcar(ndime,mnode,mgaus)                 ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)                 ! dNk/dxidxj
  real(rp)    :: gplap(mgaus,mnode)                       ! Laplacian
  real(rp)    :: gpgac(ndime,mgaus,nspec_chm)             ! Gradient(activity) / activity
  real(rp)    :: gplac(mgaus,nspec_chm)                   ! Laplacian(activity)/activity
  real(rp)    :: gpdis(mgaus)                             ! Dissipation rate 
  real(rp)    :: gpprd(mgaus,nspec_chm)                   ! Production term 
  real(rp)    :: gprrt(mgaus)                             ! Transport of reaction rate fluctuations
  real(rp)    :: ellev(mnode)
  real(rp)    :: gplev(mgaus)
  real(rp)    :: gpgle(mgaus,nspec_chm)

  real(rp)    :: dummr(mgaus*ndime)
  real(rp)    :: gpsgs(mgaus)
  real(rp)    :: gptau(mgaus)
  real(rp)    :: chale(3),chave(3),hleng(3),tragl(9)

  real(rp)    :: dtmin,dtcri

  !
  ! Loop over elements
  !
  dtmin=1.0e6_rp
  !
  elements: do ielem = 1,nelem
     !
     ! Initialization: unused variables
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
        gphes = 0.0_rp
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
     ! Add sgs or bubble
     !
     do iclas = iclai_chm,iclaf_chm
        call ADR_add_sgs_or_bubble(&
             ielem,pgaus,elmar(pelty) % shape_bub,ADR_chm(iclas),gpcon(1:pgaus,iclas:iclas,:))
     end do
     !
     ! Assemble matrix
     !  
     izmat = 1
     izrhs = 1
     do iclas = iclai_chm,iclaf_chm

        call chm_elmprc_cfi(&
                iclas,pgaus,gpcon(1:pgaus,:,:),gpden,gpgDk,gpmas,gphco,gpsph,gptur,gpdis,gpprd(1:pgaus,:), &
                gprrt,gpdif,gpgrd,gprea,gprhs)

        call chm_elmwor(&
             pnode,iclas,elcon(1:pnode,:,:),elco1) ! elco1(:,:) <= elcon(:,iclas,:)
        !
        ! Reset check
        !
        if (kfl_reset == 0) then
           call chm_adr_critical_time(pgaus,chale,gpvel,gpdif,gprea,gpmas(1:pgaus,iclas),gpden,gpcon(1:pgaus,iclas:iclas,1),dtcri)
           if (dtcri /= 0.0 ) then
              dtmin= min(dtmin,dtcri)
           endif
        endif
        ! 
        ! End reset check
        !

        if( order == ELEMENT_ASSEMBLY ) then
           ! 
           ! Assemble equation
           ! 
           call ADR_element_assembly(&
                ielem,pnode,pgaus,elcod,elmar(pelty)%shape,gpcar,elmar(pelty)%deriv,gphes,gpvol,chale,&
                elmar(pelty) % shape_bub,elmar(pelty) % deriv_bub,ADR_chm(iclas),&
                cutim,gpden,gpvel,gpdif,gpgrd,gprea,gprhs,&
                gpcon(1:pgaus,iclas:iclas,:),elcon(1:pnode,iclas:iclas,:),elmat,elrhs)

           call chm_elmdir(&
                iclas,pnode,lnods(1,ielem),elmat,elrhs)
           !
           ! Kernel assembly routines allow us to assemble each class in a single matrix,
           ! or all classes in a single monolothic matrix but one matrix at the time
           !
           if (kfl_coupl_chm==0) then
              if( solve(1) % kfl_algso == SOL_SOLVER_RICHARDSON ) then 
                 call matrix_assexp(solve(1)% ndofn, 1_ip , pnode,    npoin,lnods(1:pnode,ielem),elrhs,elmat,elco1,rhsid)
              else
                 call matrix_assrhs(solve(1)% ndofn, 1_ip , pnode,    npoin,lnods(1:pnode,ielem),elrhs,rhsid)
                 call matrix_asscsr(solve(1)% ndofn, 1_ip,r_sol,c_sol,pnode,lnods(1:pnode,ielem),elmat,amatr)
              end if
           else
              call matrix_assrhs(solve(1)% ndofn, 1_ip , pnode,    npoin,lnods(1:pnode,ielem),elrhs,rhsid,iclas)
              call matrix_asscsr(solve(1)% ndofn, 1_ip,r_sol,c_sol,pnode,lnods(1:pnode,ielem),elmat,amatr,iclas)
           endif

           izrhs = izrhs + npoin !solve(1)%nzrhs
           izmat = izmat + solve(1)%nzmat/nspec_chm**2

        else if ( order == PROJECTIONS_AND_SGS_ASSEMBLY ) then

        else if ( order == BUBBLE_ASSEMBLY ) then

        end if

     end do
     
  end do elements

  if (kfl_reset==0) then
     if (dtmin /= 0.0) then
        dtcri_chm = min(dtcri_chm, dtmin)
     endif
  endif

end subroutine chm_elmcfi

