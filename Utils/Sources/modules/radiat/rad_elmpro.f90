subroutine rad_elmpro(&
     ielem,pmate,pnode,pgaus,igaui,igauf,&
     gpsha,gpcar,gpdif,gpabs,gpbbr,gptem,gpgrd)   
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_elmpro
  ! NAME 
  !    rad_elmpro
  ! DESCRIPTION
  !    Compute properties:
  !    1. GPDIF: Total difusion constant Gamma
  !    2. GPABS: Absorption coeff, goes into Reaction term
  !    3. GPTEM: Temperature
  !    4. GPGRD: Diffusion coefficient derivative
  ! USES
  ! USED BY
  !    rad_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp           
  use def_domain, only     :  lnods,mnode,ndime,coord
  use def_radiat, only     :  tempe_rad, conce_rad, nspec_rad, kfl_parti_rad, &
       &                      scatt_rad,absor_rad,aniso_rad,kfl_atest_rad,expar_rad,steph_rad
  !use def_master, only     :  mlagr,lagrtyp
  use def_kermod
  use mod_ker_proper 

  implicit none
  integer(ip), intent(in)  :: ielem,pmate,pnode,pgaus,igaui,igauf
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),gpcar(ndime,mnode,pgaus)
  real(rp),    intent(out) :: gpdif(pgaus),gpabs(pgaus),gpbbr(pgaus),gptem(pgaus),gpgrd(ndime,pgaus)
  real(rp)                 :: ELMscattering(pnode), ELMAbsorption(pnode)  
  real(rp)                 :: ELMAnisotropyCoefficient(pnode)
  integer(ip)              :: igaus,inode,ipoin,idime,ispec,ilagr,dummi
  real(rp)                 :: ELMGamma(pnode), eltem(pnode),elden(pnode)
  real(rp),external        :: rad_gamma
  !Zero out quantities before start
  do igaus=1,pgaus
     gpdif(igaus)=0.0_rp
     gpabs(igaus)=0.0_rp
     gpbbr(igaus)=0.0_rp
     gptem(igaus)=0.0_rp
     do idime=1,ndime
        gpgrd(idime,igaus)=0.0_rp
     end do
  end do
  !
  ! Absorption and scattering coefficients
  !!F   This will need to take into account different materials
  !
  ! We get density first
  call ker_proper('DENSI','PNODE',dummi,ielem,elden,pnode,pgaus,gpsha,gpcar)
  !
  do inode=1,pnode
     ipoin=lnods(inode,ielem) 
     eltem(inode)=tempe_rad(ipoin)
     ELMscattering(inode)=0.0_rp
     ELMAbsorption(inode)=0.0_rp
     ELMAnisotropyCoefficient(inode)=0.0_rp
     do ispec=1,nspec_rad
        ELMscattering(inode) = ELMscattering(inode)+scatt_rad(ispec) *conce_rad(ipoin,ispec)*elden(inode)
        ELMAbsorption(inode) = ELMAbsorption(inode)+absor_rad(ispec) *conce_rad(ipoin,ispec)*elden(inode)
        ELMAnisotropyCoefficient(inode) = ELMAnisotropyCoefficient(inode) + scatt_rad(ispec) * aniso_rad(ispec)
     enddo
     if (kfl_parti_rad>0) then
        
     endif
     ELMGamma(inode) = rad_gamma(ELMAbsorption(inode),ELMscattering(inode),ELMAnisotropyCoefficient(inode))
  end do
  !
  !
  do inode= 1,pnode
     do igaus=igaui,igauf
        gpdif(igaus)=gpdif(igaus) &   ! Diffusion coefficient is GDIF=Gamma
             + gpsha(inode,igaus) * ELMGamma(inode)
        gpabs(igaus) =  gpabs(igaus) & ! Reaction Term: GPABS= a(x)
             + gpsha(inode,igaus) * ELMAbsorption(inode) 
        gptem(igaus)=gptem(igaus) &    ! XXX Temperature to the fourth power??
             & + gpsha(inode,igaus) * eltem(inode)
        do idime=1,ndime   ! Diffusion coefficient derivative GPGRD = grad(Gamma)
           gpgrd(idime,igaus)=gpgrd(idime,igaus)&
                +gpcar(idime,inode,igaus)*ELMGamma(inode)
        end do
     end do
  enddo
 
  if (kfl_parti_rad>0) then
!!$     do ilagr=1,mlagr
!!$        ! Check if particle belongs to this element
!!$        if (lagrtyp(ilagr) % kfl_exist == -1 .and. lagrtyp(ilagr) % ielem == ielem) then
!!$           gpbbr(igaus) = gpbbr(igaus)   &    ! Source term: GPBBR is the Black Body radiation term coming from particles
!!$                & + gpsha(inode,igaus) * &
!!$                & 4.0_rp* steph_rad  * ELMAbsorption(inode) * eltem(inode)**4         
!!$        endif
!!$     enddo
  endif
end subroutine rad_elmpro


