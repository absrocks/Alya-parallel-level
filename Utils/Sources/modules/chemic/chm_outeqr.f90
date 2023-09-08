subroutine chm_outeqr()
  !------------------------------------------------------------------------
  ! NAME 
  !    chm_outeqr
  ! DESCRIPTION
  !    Projection of gauss points gpfar to equiv_chm
  ! USES
  ! USED BY
  !    chm_outvar
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
  integer(ip) :: ielem,igaus,ipoin               ! Indices and dimensions
  integer(ip) :: pelty,pnode
  integer(ip) :: pgaus,plapl,porde,ptopo
  integer(ip) :: inode,dummi

  real(rp)    :: elcon(mnode,nspec_chm,ncomp_chm)      ! <=> conce
  real(rp)    :: elco2(mnode,nspec_chm)                ! <=> Previous iteration
  real(rp)    :: elcod(ndime,mnode)                    ! <=> coord
  real(rp)    :: elden(mnode,2)                        ! <=> global density
  real(rp)    :: eltem(mnode)                          ! <=> tempe
  real(rp)    :: elDik(mnode,nspec_chm)                ! Species diffusion coefficient
  real(rp)    :: elmas(mnode,nspec_chm)                ! Mass source terms
  real(rp)    :: elmol(mnode)                          ! Avergae molar mass
  real(rp)    :: elvel(ndime,mnode)
  real(rp)    :: elsen(mnode)                          ! Flame sensor for DTFLES
  real(rp)    :: elmut(mnode)
  real(rp)    :: gpvol(mgaus)                          ! |J|*w 
  real(rp)    :: gpcon(mgaus,nspec_chm)                ! <=> conce
  real(rp)    :: gpvel(ndime,mgaus)                    ! u
  real(rp)    :: gpDik(mgaus,nspec_chm)                ! Species difussion coefficient
  real(rp)    :: gpgDk(ndime,mgaus,nspec_chm)          ! Gradient of difussion coefficient
  real(rp)    :: gpmas(mgaus,nspec_chm)                ! Mass realease of each reaction
  real(rp)    :: gpmol(mgaus)                          ! Average molar mass
  real(rp)    :: gpgmo(ndime,mgaus)                    ! Average molar mass Gradient
  real(rp)    :: gphmo(mgaus)                          ! Average molar mass Laplacian
  real(rp)    :: gpcar(ndime,mnode,mgaus)              ! dNk/dxj
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj
  real(rp)    :: gplap(mgaus,mnode)                    ! Laplacian
  real(rp)    :: gpfar(mgaus)                          ! Fuel /air ratio
  real(rp)    :: gpspe(mgaus)                          ! Tabulated flame speed
  real(rp)    :: gpthi(mgaus)                          ! Tabulated flame thickness
  real(rp)    :: sgsef(mgaus)                          ! Subgrid scale wrinkling factor
  real(rp)    :: gpfac(mgaus)                          ! Dynamic thickening factor F for DTFLES
  real(rp)    :: gpsen(mgaus)                          ! Flame sensor DTFLES

  real(rp)    :: dummr(mgaus*ndime*mnode)
  real(rp)    :: dummw(mgaus,nspec_chm)
  real(rp)    :: chale(3),chave(3),hleng(3),tragl(9)

  !
  ! Loop over elements
  !
  elements: do ielem = 1,nelem
     !
     ! Initialization: unused variables
     !
     do igaus = 1,mgaus
        gpmol(igaus) = 0.0_rp
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
     call ker_proper('DENSI','PNODE',dummi,ielem,elden(:,1),pnode)
     call ker_proper('TURBU','PNODE',dummi,ielem,elmut,pnode)
     !

     call chm_elmgac(&
          1_ip,ielem,pnode,lnods(1,ielem),elden,elcod,elcon,elco2,elvel,eltem,elDik,&
          elmas,elmol,dummr,dummr,dummr,dummr,elsen)
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
     ! Send quantities to gauss points
     !
     call chm_elmpre(&
          pnode,pgaus,elcon(:,:,1),elden,elvel,elDik,elmas,dummr,dummr,dummr,dummr,elsen,elmar(pelty)%shape,gpcar,gplap, &
          gpcon,gpvel,gpDik,gpgDk,gpmas,elmol,gpmol,gpgmo,gphmo,gpfar,gpvol,dummr,dummr,dummw,gpspe,&
          gpthi,sgsef,gpfac,dummr,gpsen,dummr) 
     !
     ! Projection 
     ! 
     do inode = 1,pnode
        ipoin = lnods(inode,ielem)
        do igaus = 1,pgaus
           equiv_chm(ipoin) = equiv_chm(ipoin) + elmar(pelty)%shape(inode,igaus) * gpfar(igaus) * gpvol(igaus)            
           flspe_chm(ipoin) = flspe_chm(ipoin) + elmar(pelty)%shape(inode,igaus) * gpspe(igaus) * gpvol(igaus)
           flthi_chm(ipoin) = flthi_chm(ipoin) + elmar(pelty)%shape(inode,igaus) * gpthi(igaus) * gpvol(igaus)
           flsgs_chm(ipoin) = flsgs_chm(ipoin) + elmar(pelty)%shape(inode,igaus) * sgsef(igaus) * gpvol(igaus)
           flfac_chm(ipoin) = flfac_chm(ipoin) + elmar(pelty)%shape(inode,igaus) * gpfac(igaus) * gpvol(igaus)
        end do
     end do
     
          
  end do elements


end subroutine chm_outeqr

