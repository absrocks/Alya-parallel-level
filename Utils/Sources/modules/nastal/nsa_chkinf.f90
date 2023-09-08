!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_chkinf.f90
!> @author  Mariano Vazquez
!> @date    07/11/2013
!> @brief   Check wether a boundary node is inflow or outflow 
!> @details Check wether a boundary node is inflow or outflow 
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_chkinf(kinfl,xmach,ipoin,xadve,xvelo,xpres,xdens)
  use def_master
  use def_domain
  use def_nastal
  use def_kermod
  use mod_ker_proper

  implicit none
  integer(ip)   :: kinfl
  integer(ip)   :: ipoin
  integer(ip)   :: idime,ibopo,dummi,kpoin
  real(rp)      :: xmach
  real(rp)      :: xdens
  real(rp)      :: xpres
  real(rp)      :: xvelo(ndime),xadve(ndime)
  real(rp)      :: uproj,uveno(ndime),dummy(ndime,ndime),adgam,sound,vmodu,vadve

  kinfl= 0

  if (kfl_prope /= 0 ) then
     !           wmean = mowei_nsa
     !call ker_proper('SPHEA','NPOIN',&
     !     dummi,dummi,shecp_nsa,dummi,dummi,dummy(:,1),dummy(:,1))
     call ker_proper('SPHEA','NPOIN',&
          dummi,dummi,shecp_nsa)
     adgam = shecp_nsa(ipoin) / (shecp_nsa(ipoin) - runiv_nsa/ wmean(ipoin,1))
  else
     adgam = cpcoe_nsa / (cpcoe_nsa - runiv_nsa/ mowei_nsa)
  endif
  
  ! Compute mach number, to see if sub or supersonic
  
  sound = sqrt(adgam * xpres / xdens)
  vmodu = xvelo(1)*xvelo(1) + xvelo(2)*xvelo(2)
  if (ndime == 3) vmodu = vmodu + xvelo(3)*xvelo(3)
  vmodu = sqrt(vmodu)
  xmach = vmodu/sound

  vadve = xadve(1)*xadve(1) + xadve(2)*xadve(2)
  if (ndime == 3) vadve = vadve + xadve(3)*xadve(3)
  vadve = sqrt(vadve)

  uveno=0.0_rp
  if (vadve > 1.0e-10) then
     do idime= 1,ndime
        uveno(idime)= xadve(idime)/vadve
     end do
  end if

  ibopo= lpoty(ipoin)
  
  uproj= 0.0_rp
  do idime= 1,ndime
     uproj= uproj+uveno(idime)*exnor(idime,1,ibopo)
  end do
  if (uproj <= 0.0_rp) kinfl= 1  ! inflow detected

end subroutine nsa_chkinf
