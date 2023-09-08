subroutine chm_elmcom()
  !------------------------------------------------------------------------
  ! NAME 
  !    chm_elmcom
  ! DESCRIPTION
  !    Elemental operations for standard combustion models
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

  implicit none
  real(rp)    :: elmat(mnode,mnode) 
  real(rp)    :: elrhs(mnode)
  integer(ip) :: ielem,igaus,idime,iclas               ! Indices and dimensions
  integer(ip) :: izmat,izrhs,pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo
  integer(ip) :: dummi
  integer(ip) :: inode

  real(rp)    :: elcon(mnode,nspec_chm,ncomp_chm)      ! <=> conce
  real(rp)    :: elco2(mnode,nspec_chm)                ! <=> Previous iteration
  real(rp)    :: elco1(mnode,ncomp_chm)                ! <=> conce: work array
  real(rp)    :: elcod(ndime,mnode)                    ! <=> coord
  real(rp)    :: elden(mnode,2)                        ! <=> global density
  real(rp)    :: eltem(mnode)                          ! <=> tempe
  real(rp)    :: elmut(mnode)                          ! Turbulence viscosity
  real(rp)    :: elsen(mnode)                          ! Flame sensor for DTFLES
  real(rp)    :: elDik(mnode,nspec_chm)                ! Species diffusion coefficient
  real(rp)    :: elmas(mnode,nspec_chm)                ! Mass source terms
  real(rp)    :: elmol(mnode)                          ! Avergae molar mass
  real(rp)    :: elvel(ndime,mnode)
  real(rp)    :: gpvol(mgaus)                          ! |J|*w 
  real(rp)    :: gpcon(mgaus,nspec_chm)                ! <=> conce
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
  real(rp)    :: gppro(mgaus)                          ! Weighted residual L2-projection
  real(rp)    :: gpdiv(mgaus)                          ! Divergence of convection
  real(rp)    :: gpmas(mgaus,nspec_chm)                ! Mass realease of each reaction
  real(rp)    :: gpmol(mgaus)                          ! Average molar mass
  real(rp)    :: gpgmo(ndime,mgaus)                    ! Average molar mass Gradient
  real(rp)    :: gphmo(mgaus)                          ! Average molar mass Laplacian
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj
  real(rp)    :: gplap(mgaus,mnode)                    ! Laplacian
  real(rp)    :: gpgac(ndime,mgaus,nspec_chm)          ! Gradient(activity) / activity
  real(rp)    :: gplac(mgaus,nspec_chm)                ! Laplacian(activity)/activity
  real(rp)    :: gpfar(mgaus)                          ! Fuel / air ratio
  real(rp)    :: gpspe(mgaus)                          ! Tabulated flame speed
  real(rp)    :: gpthi(mgaus)                          ! Tabulated flame thickness
  real(rp)    :: sgsef(mgaus)                          ! Subgrid scale wrinkling factor for TFLES
  real(rp)    :: gpfac(mgaus)                          ! Dynamic thickening factor F for DTFLES
  real(rp)    :: gpsen(mgaus)                          ! Flame sensor DTFLES
  real(rp)    :: gptur(mgaus)                          ! Flame sensor DTFLES

  real(rp)    :: ellev(mnode)
  real(rp)    :: gplev(mgaus)
  real(rp)    :: gpgle(mgaus,nspec_chm)

  real(rp)    :: dummr(mgaus*ndime*mnode)
  real(rp)    :: dummw(mgaus,nspec_chm)
  real(rp)    :: gpsgs(mgaus)
  real(rp)    :: gptau(mgaus)
  real(rp)    :: chale(3),chave(3),hleng(3),tragl(9)

  real(rp)    :: dtmin,dtcri
  !
  ! Loop over elements
  !
  dtmin=1.e6
  !

  elements: do ielem = 1,nelem
     !
     ! Initialization: unused variables
     !
     do igaus = 1,mgaus
        gptau(igaus) = 0.0_rp
        gpdiv(igaus) = 0.0_rp
        gprhs(igaus) = 0.0_rp
        gpdif(igaus) = 0.0_rp
        gprea(igaus) = 0.0_rp
        gpsgs(igaus) = 0.0_rp
        gppro(igaus) = 0.0_rp
        gpmol(igaus) = 0.0_rp
        do idime = 1,ndime
           gpgrd(idime,igaus) = 0.0_rp
        end do
     end do
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
     call ker_proper('DENSI','PNODE',dummi,ielem,elden(:,1),pnode,pgaus)
     call ker_proper('TURBU','PNODE',dummi,ielem,elmut,pnode,pgaus)
     !
     call chm_elmgac_species(&
          3_ip,ielem,pnode,lnods(1,ielem),elden,elcod,elcon,elco2,elvel,eltem,elDik,&
          elmas,elmol,elmut,elsen)
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
     ! Send quantities to gauss points
     !
     call chm_elmpre_species(&
          pnode,pgaus,elcon(:,:,1),elden,elvel,elDik,elmas,elmut,elsen,elmar(pelty)%shape,gpcar,gplap, &
          gpcon,gpvel,gpDik,gpgDk,gpmas,elmol,gpmol,gpgmo,gphmo,gpfar,gpvol,gpdiv,gpspe,&
          gpthi,sgsef,gpfac,gpsen,gptur) 

     !
     ! Density from kermod
     !
     call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty) % shape,gpcar)
     call ker_proper('GRDEN','PGAUS',dummi,ielem,gpgde,pnode,pgaus,elmar(pelty) % shape,gpcar)
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
     ! pre-compute correction velocity to ensure total mass conservation
     !
     if (kfl_gauss_chm == 0_ip) then
        call chm_elmvel(&
             pnode,pgaus,elco2,gpDik,gpgDk,elmar(pelty)%shape,gpcar,gplap,&
             gpvec,gpgve,gpmol,gpgmo,gphmo)
        if (kfl_activ_chm == 1) call chm_levvel(&  ! Correct with activities
             pnode,pgaus,elco2,gpDik,gpgDk,elmar(pelty)%shape,gpcar,gplap,&
             gpvec,gpgve,gpmol,gpgmo,gphmo,gpgac,gplac)
     else
        call chm_elmvel(&
             pnode,pgaus,elcon,gpDik,gpgDk,elmar(pelty)%shape,gpcar,gplap,&
             gpvec,gpgve,gpmol,gpgmo,gphmo)
        if (kfl_activ_chm == 1) call chm_levvel(&
             pnode,pgaus,elcon,gpDik,gpgDk,elmar(pelty)%shape,gpcar,gplap,&
             gpvec,gpgve,gpmol,gpgmo,gphmo,gpgac,gplac)
     endif
     !
     ! Assemble matrix
     !  
     izmat = 1
     izrhs = 1
     do iclas = iclai_chm,iclaf_chm
        call chm_elmprc(&
                iclas,pgaus,gpcon,gpden,gpgde,gpDik,gpgDk,gpmas,gpvel,gpvec,gpgve,&
                gpadv,gpdif,gpgrd,gprea,gprhs,gpmol,gpgmo,gphmo,gpgac,gplac)
        call chm_elmwor(&
                pnode,iclas,elcon,elco1) ! elco1(:,:) <= elcon(:,iclas,:)
        !
        ! Reset check
        !
        if (kfl_reset == 0) then
           call chm_adr_critical_time(pgaus,chale,gpadv,gpdif,gprea,gpmas(:,iclas),gpden,gpcon(:,iclas),dtcri)
           if (dtcri /= 0.0 ) then
              dtmin= min(dtmin,dtcri)
           endif
        endif

        !
        ! Addition of extra artificial diffusion for undershoots/overshoots
        !
!        if ( kfl_shock_chm > 0.0_rp) then
!           rhok = 0.0_rp  !  averaged value of element mass fraction
!
!           do igaus = 1,mgaus
!              rhok = rhok + gpcon(igaus,iclas)
!           enddo
!
!           rhok = rhok / dble(mgaus)
!
!           if (rhok > 1.0_rp .or. rhok < 0.0_rp) then
!              shock_chm = min(3.0_rp, shock_chm*2.0_rp) 
!           endif
!
!        endif

        ! 
        ! End reset check
        !
        call elmadr(& ! Elemental assembly
             1_ip,ielem,pnode,pgaus,ptopo,plapl,pelty,1_ip,1_ip,lnods(1,ielem),kfl_shock_chm,&
             kfl_taust_chm,kfl_stabi_chm,0_ip,0_ip,kfl_limit_chm,staco_chm,gpdif,&
             gprea,gpden,gpadv,gppro,gpvol,gpgrd,gprhs,elmar(pelty)%shape,gpcar,gphes,&
             elco1,elcod,chale,gptau,dtinv_chm,shock_chm,bemol_chm,gpdiv,gplap,gpsgs,&
             elmat,elrhs)
        call chm_elmdir(&
             iclas,pnode,lnods(1,ielem),elmat,elrhs)
        !
        ! Kernel assembly routines allow us to assemble each class in a single matrix,
        ! or all classes in a single monolothic matrix but one matrix at the time
        !

        if (kfl_coupl_chm==0) then
           !
           ! SHOULD UPDATE TO
           ! matrix_assemble_element_RHS
           ! matrix_assemble_element_matrix_to_CSR
           !
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

!!$        call assrhs(&
!!$             solve(1)%ndofn,pnode,lnods(1,ielem),elrhs,rhsid(izrhs))
!!$        call assmat(&
!!$             solve(1)%ndofn,pnode,pnode,npoin,solve(1)%kfl_algso,&
!!$             ielem,lnods(1,ielem),elmat,amatr(izmat))
        izrhs = izrhs + npoin !solve(1)%nzrhs
        izmat = izmat + solve(1)%nzmat/nspec_chm**2
     end do

  end do elements

  if (kfl_reset==0) then
     if (dtmin /= 0.0) then
        dtcri_chm = min(dtcri_chm, dtmin)
     endif
  endif

end subroutine chm_elmcom

