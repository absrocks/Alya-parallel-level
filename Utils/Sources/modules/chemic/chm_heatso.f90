subroutine chm_heatso() 
  !-----------------------------------------------------------------------
  !****f* partis/chm_heatso
  ! NAME 
  !    chm_heatso
  ! DESCRIPTION
  !    Compute heat release (source term that ends up in Temper and Nastal)
  ! USES
  ! USED BY
  !   
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_master
  use def_domain
  use def_chemic
  use mod_ker_proper 
  use def_kermod

  implicit none
  integer(ip) :: ielem,pnode
  integer(ip) :: dummi
  real(rp)    :: gpcar(ndime,mnode,mgaus)
  real(rp)    :: gphes(ntens,mnode,mgaus)              ! dNk/dxidxj
  real(rp)    :: gpvol(mgaus)                          ! |J|*w 
  real(rp)    :: elcod(ndime,mnode)                    ! <=> coord
  real(rp)    :: elden(mnode,2)                        ! <=> global density
  real(rp)    :: elcon(mnode,nspec_chm,ncomp_chm)      ! <=> conce
  real(rp)    :: elco2(mnode,nspec_chm)                ! <=> Previous iteration
  real(rp)    :: eltem(mnode)                          ! <=> tempe
  real(rp)    :: elDik(mnode,nspec_chm)                ! Species diffusion coefficient
  real(rp)    :: elmas(mnode,nspec_chm)                ! Mass source terms
  real(rp)    :: elmut(mnode)                          ! Turbulence viscosity
  real(rp)    :: elmol(mnode)                          ! Avergae molar mass
  real(rp)    :: elsen(mnode)                          ! Flame sensor
  real(rp)    :: gpfar(mgaus)                          ! Fuel / air ratio
  real(rp)    :: gpdiv(mgaus)                          ! Divergence of convection
  real(rp)    :: gpspe(mgaus)                          ! Tabulated flame speed
  real(rp)    :: gpthi(mgaus)                          ! Tabulated flame thickness
  real(rp)    :: elvel(ndime,mnode)
  real(rp)    :: gpcon(mgaus,nspec_chm)
  real(rp)    :: gpgco(ndime,mgaus,nspec_chm)                ! <=> conce
  real(rp)    :: gp_lap_conce(mgaus,nspec_chm)  
  real(rp)    :: gpden(mgaus),gpgde(ndime,mgaus)       ! fake rho for elmadr
  real(rp)    :: gpmas(mgaus,nspec_chm) ! Mass source term
  real(rp)    :: gpDik(mgaus,nspec_chm)                ! Species difussion coefficient
  real(rp)    :: gpgDk(ndime,mgaus,nspec_chm)          ! Gradient of difussion coefficient
  real(rp)    :: gpvel(ndime,mgaus)                    ! u
  integer(ip) :: igaus,inode,idime
  real(rp)    :: enthalpy(mgaus,nspec_chm)    ! Enthalpy per species
  real(rp)    :: grad_enthalpy(ndime,mgaus,nspec_chm)    ! Enthalpy per species
  real(rp)    :: Cp(mgaus,nspec_chm)
  real(rp)    :: gpgrt(ndime,mgaus)
  real(rp)    :: dummr(mgaus*mnode*ndime)
  real(rp)    :: dummw(mgaus,nspec_chm)
  real(rp)    :: gplap(mgaus,mnode)                    ! Laplacian
  integer(ip) :: pelty
  integer(ip) :: pgaus,plapl,porde,ptopo  

!!!! FER   Awful awful routine I must take out some gathering operations and optimize

  if (INOTMASTER) then
     if (kfl_coupl(ID_TEMPER,ID_CHEMIC) >= 1 .or. kfl_coupl(ID_NASTAL,ID_CHEMIC) >= 1 ) then
        do ielem=1,nelem
           !
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           pgaus = ngaus(pelty)
           plapl = llapl(pelty) 
           porde = lorde(pelty)
           ptopo = ltopo(pelty)
           !
           gpgrt    = 0.0_rp
           cp       = 0.0_rp
           gpmas    = 0.0_rp
           gpdik    = 0.0_rp
           enthalpy = 0.0_rp
           gpgco    = 0.0_rp
           gpdiv    = 0.0_rp
           !
           !
           call ker_proper('DENSI','PNODE',dummi,ielem,elden(:,1),pnode)
           elvel    = 0.0_rp
           call chm_elmgac(&
                1_ip,ielem,pnode,lnods(1,ielem),elden,elcod,elcon,elco2,elvel,eltem,elDik,&
                elmas,elmol,elmut,dummr,dummr,dummr,elsen)
           !
           ! Geometric information
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
              gplap=0.0_rp
           endif
           !
           ! Send quantities to gauss points
           !
           call ker_proper('DENSI','PGAUS',dummi,ielem,gpden,pnode,pgaus,elmar(pelty) % shape,gpcar)
           call ker_proper('DENSI','GRDEN',dummi,ielem,gpgde,pnode,pgaus,elmar(pelty) % shape,gpcar)

           call chm_elmpre(&
                pnode,pgaus,elcon(:,:,1),elden,elvel,elDik,elmas,elmut,dummr,dummr,dummr,elsen,elmar(pelty)%shape,gpcar,gplap, &
                gpcon,gpvel,gpDik,gpgDk,gpmas,elmol,dummr,dummr,dummr,gpfar,gpvol,gpdiv,dummr,dummw,gpspe,&
                gpthi,tfles_sgseff(ielem)%a,tfles_factor(ielem)%a,dummr,tfles_sensor(ielem)%a,dummr) 

           call chm_heapre(&
                ielem,pnode,pgaus,lnods(1,ielem),elcon,eltem,elmar(pelty)%shape,gpcar,gplap, &
                enthalpy,cp,gpgco,gpgrt, grad_enthalpy, gp_lap_conce)

           ! FINALLY
           ! We put the heat source into the global variable used by other modules
           ! It is done inside a subroutine because of the pnode vs mnode declaration
           ! in different routines...
           if (kfl_coupl(ID_CHEMIC, ID_TEMPER)==1_ip) then
              call chm_temper(&
                   pgaus,enthalpy,gpmas,gpden,gpgrt,cp,gpdik,gpcon,gpgco,&
                   div_enthalpy_transport(ielem)%a,chemical_heat(ielem)%a)
           else if (kfl_coupl(ID_CHEMIC, ID_NASTAL)==1_ip) then
              call chm_nastal(&
                   pgaus, enthalpy, grad_enthalpy, gpmas, gpdik, gpgDk, &
                   gpcon, gpgco,gp_lap_conce, enthalpy_transport(ielem)%a,div_enthalpy_transport(ielem)%a, chemical_heat(ielem)%a)
           endif

        enddo
     endif
  endif
  
end subroutine chm_heatso
