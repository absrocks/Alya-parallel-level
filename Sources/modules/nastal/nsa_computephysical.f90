!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_computephysical.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Compute physical variables, always from the conservative set at icomp=ITER_K
!> @details Compute physical variables, always from the conservative set at icomp=ITER_K
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_computephysical(imode)
  use      def_master
  use      def_domain
  use      def_parame
  use      def_kermod
  use      mod_ker_proper
  use      def_nastal

  implicit none
  integer(ip)       :: imode !> Initial = 0 or Running = 1
  integer(ip)       :: ipoin,kpoin,idime,kfixi,idofn,idummy,icomp
  real(rp)          :: xrano(3),xhecv,xhecp, adgam, vmodu, &
       tempe_old, tempe_new, densi_old, densi_new, press_old, press_new, &
       veloc_new(3),umome_new(3),energ_new,rdummy,xconc

  ! Compute physical variables, always from the conservative set at icomp=ITER_K
  icomp= ITER_K

  rdummy= 0.0_rp
  xconc = 0.0_rp

  if (INOTMASTER) then                                            

     if (kfl_prope /= 0 ) then
        call ker_proper('SPHEA','NPOIN',idummy,idummy,shecp_nsa)
        call ker_proper('VISCO','NPOIN',idummy,idummy,visco(:,icomp))
        !call ker_proper('SPHEA','NPOIN',idummy,idummy,shecp_nsa,idummy,idummy,xrano,xrano)
        !call ker_proper('VISCO','NPOIN',idummy,idummy,visco(:,icomp),idummy,idummy,xrano,xrano)
        if (kfl_coupl(ID_NASTAL,ID_CHEMIC) == 0 ) then
           wmean = mowei_nsa
        endif
     else
        wmean     = mowei_nsa
        shecp_nsa = cpcoe_nsa
     endif

     do ipoin = 1,npoin

        densi_new   = densi(ipoin,icomp) 
        energ_new   = energ(ipoin,icomp)
        umome_new(1:ndime) = umome(1:ndime,ipoin,icomp) 
        xhecp = shecp_nsa(ipoin)                    ! C_p
        xhecv = xhecp - runiv_nsa / wmean(ipoin,1)  ! C_v
        adgam = xhecp / xhecv                       ! gamma = C_p / C_v
        veloc_new(1:ndime) = umome_new(1:ndime) / densi_new
        vmodu = veloc_new(1) * veloc_new(1) + veloc_new(2) * veloc_new(2)
        if (ndime == 3) vmodu = vmodu + veloc_new(3) * veloc_new(3)
        
        tempe_old   = tempe(ipoin,icomp) 
        press_old   = press(ipoin,icomp) 

        !check if there is any energy prescription
        kfixi= kfl_fixno_nsa(ndime+2,ipoin)
        if (kfixi <= 1) then 
           ! no energy prescription or energy is fixed: compute T from E
           tempe_new   = (energ_new   /  densi_new  - 0.5_rp * vmodu) / xhecv
           if (kfixi ==1) tempe_new   = bvess_nsa(ndime+2,ipoin,1)
           ! compute p from rho and T
           call nsa_stalaw(2,imode,densi_new,press_new,tempe_new,&
                rdummy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))           

        else if (kfixi == 3) then
           ! temperature is fixed: compute T from E 
           ! due to the BC prescription in nsa_elmatrixsetboundary, T retains here its bvess
!!           tempe_new   = (energ_new   /  densi_new  - 0.5_rp * vmodu) / xhecv
           tempe_new   = tstag_nsa
           ! compute p from rho and T
           call nsa_stalaw(2,imode,densi_new,press_new,tempe_new,&
                rdummy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))           

        else if (kfixi == 2) then
           ! pressure is fixed: compute p from E 
           ! due to the BC prescription in nsa_elmatrixsetboundary, p retains here its bvess
!!           press_new   = (energ_new - 0.5_rp * densi_new * vmodu) * (adgam - 1.0_rp)
           press_new   = bvess_nsa(ndime+2,ipoin,1)
           ! compute T from rho and p
           call nsa_stalaw(3,imode,densi_new,press_new,tempe_new,&
                rdummy,xconc,wmean(ipoin,1),shecp_nsa(ipoin))           

        end if
        
        if (tempe_new < 0.0_rp) then
           write(6,*) 'tempe podrida=', ipoin, tempe_old, tempe_new              
        end if

        press(ipoin,icomp) = press_new
        tempe(ipoin,icomp) = tempe_new
        veloc(1:ndime,ipoin,icomp) = veloc_new(1:ndime)


!        if (coord(1,ipoin) < -4.99999) then
!           write (6,*) ipoin,press(ipoin,icomp), tempe(ipoin,icomp)
!        end if

     end do

  end if


end subroutine nsa_computephysical
