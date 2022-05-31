!------------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_gauvalxy.f90
!> @date    01/08/2012
!> @author  Mariano Vazquez
!> @brief   Computes all the variables on the Gauss points
!> @details Computes all the variables on the Gauss points.\n
!!          The unknowns are interpolated from its value on the nodes using the test functions.\n
!!          The other variables are directely computed on the Gauss points from the value of the unknowns.\n
!!          The mean value of the Jacobian (conme) and diffusion (difme) matrix in element ielem\n
!!          is also computed at the end of this subroutine.!!\n
!> @}
!------------------------------------------------------------------------
subroutine nsa_gauvalxy(&
     ielem,igaus,pnode,pgaus,weigh,&
     elcon,eldif,elunk,elsub,elort,eldtt,elmsh,&
     xunkn,xdtix,elpre,&
     xconv,dconv,xdiff,ddiff,&
     xconv_der,&
     gunkn,dflux_conv,gsube, &
     xshap,cartd,hesma,&
     hessi,xsube,xortp,xresi, &
     sound,xpres,xtemp,xvisc,xdith,xvelo,xvmsh,&
     gpres,gtemp, &
     gvisc,gvelo,dvelo,velmo,xlade,xldve,&
     xjaci,conme,difme,dvolu,xnutu,xvofo,htrad,dhtra,&
     xmowe,heats,xtide,xlopr_conservative,hleng)
  use      def_master
  use      def_domain
  use      def_nastal
  use      def_kermod
  use      mod_ker_proper

  implicit none

  integer(ip), intent(in):: pnode,igaus,pgaus,ielem
  integer(ip) :: idime,kdime,idofn,inode,jdofn,jdime,ipoin,kdofn,ldofn,mdofn,odofn,pdofn,ievat

  real(rp),    intent(in):: dvolu,weigh,elcon(ndofn_nsa,ndofn_nsa,ndime,mnode),eldif(ndofn_nsa,ndofn_nsa,ndime,ndime,mnode), &
       elunk(ndofn_nsa,mnode,ncomp_nsa),elsub(ndofn_nsa,mnode),elort(ndofn_nsa,mnode),elmsh(ndime,mnode), &
       eldtt(ndofn_nsa,mnode,2),elpre(mnode),hessi(ntens,mnode),xshap(mnode),cartd(ndime,mnode), &
       xjaci(ndime,ndime),htrad(ndime),dhtra,heats,hleng(ndime)
  real(rp),   intent(out):: xnutu,xvofo(ndofn_nsa,ndofn_nsa),xunkn(ndofn_nsa,mgaus,3),xdtix(ndofn_nsa,mgaus,2), &
       xconv(ndofn_nsa,ndofn_nsa,ndime),xdiff(ndofn_nsa,ndofn_nsa,ndime,ndime), &
       xconv_der(ndofn_nsa,ndofn_nsa,ndofn_nsa,ndime),&
       dconv(ndofn_nsa,ndofn_nsa),ddiff(ndofn_nsa,ndofn_nsa,2,ndime),hesma(ndime,ndime,mnode), &
       gunkn(ndofn_nsa,ndime),dflux_conv(ndofn_nsa),gsube(ndofn_nsa,ndime),xsube(ndofn_nsa,mgaus,3),xortp(ndofn_nsa), &
       xresi(ndofn_nsa),xvelo(ndime),xvmsh(ndime),gpres(ndime),gtemp(ndime),gvisc(ndime),gvelo(ndime,ndime), &
       xldve(ndime),xdith,xtemp,xvisc,velmo,xpres,dvelo,sound,xlade,xmowe, &
       xtide(ndofn_nsa) ! Preconditioner * ( Unknown(last_iteration) - Unknown(last_time_step) ) / physical_time_step
  real(rp),  intent(inout):: conme(ndofn_nsa,ndofn_nsa,ndime),difme(ndofn_nsa,ndofn_nsa,ndime,ndime)

  real(rp)                :: prope_tmp(1),gtunk(ndofn_nsa,ndime),xshai,rgasc,rgacv, &
       gpreo(ndime),xpreo,dvite,velsq,enepe,dicod,visci,rdumy, &
       hunkn(ndofn_nsa,ndime,ndime),xtunk(ndofn_nsa),xhecp,xhecv,adgam,dwall, &
       sgsdi,xmile,seci4, &
       xconv_test_1(ndofn_nsa,ndofn_nsa,ndime), & ! To test!!!!!
       xconv_test_2(ndofn_nsa,ndofn_nsa,ndime), & ! To test!!!!!
       xconv_ori(ndofn_nsa,ndofn_nsa,ndime), & ! Original convective Jacobian without preconditioning
       xconv_lopre_diff(ndofn_nsa,ndofn_nsa,ndime), & ! dP * K for CM local preconditioner
       xdiff_ori(ndofn_nsa,ndofn_nsa,ndime,ndime), & ! Original diffusion matrix without preconditioning
       xttra(ndofn_nsa,ndofn_nsa), & ! Symmetrizing transformation of the N-S equations
       xttra_inv(ndofn_nsa,ndofn_nsa), & ! Inverse of the Symmetrizing transformation of the N-S equations
       xqtra(ndofn_nsa,ndofn_nsa), & ! Streamline transformation of the N-S equations
       xqtra_inv(ndofn_nsa,ndofn_nsa), & ! Inverse of the Streamline transformation of the N-S equations
       xftra(ndofn_nsa,ndofn_nsa), & ! Symmetrizing streamline transformation of the N-S equations
       xftra_inv(ndofn_nsa,ndofn_nsa), & ! Inverse of the Symmetrizing streamline transformation of the N-S equations
       xlopr(ndofn_nsa,ndofn_nsa), & ! Local preconditioner in the symmetrizing variables and streamline coordinates
       xlopr_conservative(ndofn_nsa,ndofn_nsa), & ! Local preconditioner in the conservative variables and cartesian coordinates
       velmo_xy, & ! sqrt(u_1^2 + u_2^2)
       xmach, & ! Mach number at Guass points
       xmach_ref, & ! Reference Mach number at Guass points
       xbeta, &  ! local preconditioning parameter
       xbeta_visc, &  ! local preconditioning parameter
       xtaup, & ! local preconditioning parameter
       cfl_sound, & ! local preconditioning parameter
       reyno_cell, reyno_cell_inv, & ! cell Reynolds number and its inverse
       gmach(ndime), gsoun(ndime), gmach_ref(ndime), &
       gvemo(ndime), gbeta(ndime), gbeta_visc(ndime), greyn_cell_inv(ndime), &
       glopr(ndofn_nsa,ndofn_nsa,ndime), & ! gradient of CM preconditioner
       termA,termB,termC,termD,termA_der,termB_der,termC_der,termD_der, veloj, momei, energy_pressure, &
       termE,termF,termG,termH,termI,term1,term2,term3,term4,term5,zensa_lopre,hmini,hmaxi,factor_skew

  xlade       = 0.0_rp
  xpreo       = 0.0_rp       ! old pressure and pressure gradient: only 
  xnutu       = 0.0_rp
  rdumy       = 0.0_rp
  dwall       = 0.0_rp
  xmowe       = 0.0_rp
  zensa_lopre = 0.1_rp
  !  zensa_lopre= 0.01_rp

  factor_skew = 1.0_rp
  if (kfl_skews_nsa == 1) factor_skew = 0.5_rp

  hmini = hleng(ndime)  ! hleng(ndime) is the smallest
  hmaxi = hleng(1)      ! hleng(1) is the largest 

  do idofn= 1,ndofn_nsa
     dflux_conv(idofn)    = 0.0_rp
     xdtix(idofn,igaus,1) = 0.0_rp
     xdtix(idofn,igaus,2) = 0.0_rp
     xtide(idofn) = 0.0_rp
     do jdofn=1,ndofn_nsa
        xttra(idofn,jdofn) = 0.0_rp
        xttra_inv(idofn,jdofn) = 0.0_rp
        xqtra(idofn,jdofn) = 0.0_rp
        xqtra_inv(idofn,jdofn) = 0.0_rp
        xlopr(idofn,jdofn) = 0.0_rp
        xlopr_conservative(idofn,jdofn) = 0.0_rp

        xftra(idofn,jdofn) = 0.0_rp
        xftra_inv(idofn,jdofn) = 0.0_rp

        do idime=1,ndime
           xconv(idofn,jdofn,idime) = 0.0_rp           
           xconv_ori(idofn,jdofn,idime) = 0.0_rp           
           xconv_test_1(idofn,jdofn,idime) = 0.0_rp           
           xconv_test_2(idofn,jdofn,idime) = 0.0_rp           
           xconv_lopre_diff(idofn,jdofn,idime) = 0.0_rp
           ddiff(idofn,jdofn,1,idime) = 0.0_rp           
           ddiff(idofn,jdofn,2,idime) = 0.0_rp           
           do jdime=1,ndime
              xdiff(idofn,jdofn,idime,jdime) = 0.0_rp           
           end do
        end do
        dconv(idofn,jdofn) = 0.0_rp
        xvofo(idofn,jdofn)    = 0.0_rp
     end do

     xsube(idofn,igaus,2) = 0.0_rp
     xunkn(idofn,igaus,1) = 0.0_rp
     xunkn(idofn,igaus,2) = 0.0_rp
     xunkn(idofn,igaus,3) = 0.0_rp
     xresi(idofn) = 0.0_rp
     xortp(idofn) = 0.0_rp
  end do

  do idime=1,ndime
     xldve(idime) = 0.0_rp
     gpres(idime) = 0.0_rp
     gpreo(idime) = 0.0_rp !     used as limiter when current pressure goes below zero
     xvmsh(idime) = 0.0_rp
     do jdime=1,ndime
        gvelo(idime,jdime) = 0.0_rp
     end do
     do idofn=1,ndofn_nsa
        gunkn(idofn,idime) = 0.0_rp
        gsube(idofn,idime) = 0.0_rp
        do jdofn=1,ndofn_nsa
           glopr(idofn,jdofn,idime) = 0.0_rp
        end do
        do jdime=1,ndime
           hunkn(idofn,idime,jdime) = 0.0_rp
        end do
     end do
     xsube(idime,igaus,2)= umosg_nsa(idime,ielem,igaus,1)
     xsube(idime,igaus,3)= umosg_nsa(idime,ielem,igaus,2)
  end do

  xsube(ndime+1,igaus,2)= densg_nsa(ielem,igaus,1)
  xsube(ndime+2,igaus,2)= enesg_nsa(ielem,igaus,1)
  xsube(ndime+1,igaus,3)= densg_nsa(ielem,igaus,2)
  xsube(ndime+2,igaus,3)= enesg_nsa(ielem,igaus,2)  

  do inode= 1,pnode
     xshai= xshap(inode)
     xpreo= xpreo + xshai * elpre(inode)
     do idofn= 1,ndofn_nsa
        !!        xsube(idofn,igaus,2) = xsube(idofn,igaus,2) + xshai * elsub(idofn,inode)
        xunkn(idofn,igaus,ITER_K) = xunkn(idofn,igaus,ITER_K) &
             + xshai*elunk(idofn,inode,ITER_K)
        xunkn(idofn,igaus,ITER_AUX) = xunkn(idofn,igaus,ITER_AUX) &
             + xshai*elunk(idofn,inode,ITER_AUX)
        xunkn(idofn,igaus,TIME_N) = xunkn(idofn,igaus,TIME_N) + xshai*elunk(idofn,inode,TIME_N)
        xdtix(idofn,igaus,1) = xdtix(idofn,igaus,1) + xshai*eldtt(idofn,inode,1)     
        xdtix(idofn,igaus,2) = xdtix(idofn,igaus,2) + xshai*eldtt(idofn,inode,2)
        xortp(idofn) = xortp(idofn) + xshai * elort(idofn,inode)

        ievat = (inode-1) * ndofn_nsa + idofn

        do jdofn=1,ndofn_nsa
           do idime=1,ndime
              dconv(idofn,jdofn) = dconv(idofn,jdofn) &
                   + cartd(idime,inode)*elcon(idofn,jdofn,idime,inode)
              do jdime=1,ndime
                 ddiff(idofn,jdofn,1,idime) = ddiff(idofn,jdofn,1,idime) &
                      + cartd(jdime,inode) * eldif(idofn,jdofn,jdime,idime,inode)
                 ddiff(idofn,jdofn,2,idime) = ddiff(idofn,jdofn,2,idime) &
                      + cartd(jdime,inode) * eldif(idofn,jdofn,idime,jdime,inode)
              end do
           end do
        end do
     end do

     if (kfl_skews_nsa == 1) then           
        do idime=1,ndime
           do jdime=1,ndime           
              veloj= elunk(jdime,inode,ITER_K) / elunk(ndime+1,inode,ITER_K) 
              momei= elunk(idime,inode,ITER_K)
              dflux_conv(idime) = dflux_conv(idime) + cartd(jdime,inode) * (veloj * momei)
           end do
           dflux_conv(idime) = dflux_conv(idime) + elpre(inode) * cartd(idime,inode)
           
           dflux_conv(ndime+1) = dflux_conv(ndime+1) + cartd(idime,inode) * elunk(idime,inode,ITER_K)
           veloj= elunk(idime,inode,ITER_K) / elunk(ndime+1,inode,ITER_K)  
           energy_pressure = elunk(ndime+2,inode,ITER_K) + elpre(inode)
           dflux_conv(ndime+2) = dflux_conv(ndime+2) + cartd(idime,inode) * veloj * energy_pressure
        end do
     end if
     
     do idime=1,ndime

        gpreo(idime) = gpreo(idime) - elpre(inode) * cartd(idime,inode)

        do jdime=1,ndime
           hesma(idime,jdime,inode) = hessi(nindx_nsa(idime,jdime),inode)
        end do
        do idofn = 1,ndofn_nsa
           gunkn(idofn,idime) = gunkn(idofn,idime) &
                + cartd(idime,inode)*elunk(idofn,inode,ITER_K)
           gsube(idofn,idime) = gsube(idofn,idime) + cartd(idime,inode)*elsub(idofn,inode)
           do jdime=1,ndime
              hunkn(idofn,idime,jdime) = hunkn(idofn,idime,jdime) + hesma(idime,jdime,inode) * elunk(idofn,inode,ITER_K)
           end do
        end do
        xlade        = xlade        + hessi(idime,inode) * elunk(ndime+1,inode,ITER_K) 
        xldve(idime) = xldve(idime) + hessi(idime,inode) * elunk(ndime+1,inode,ITER_K)
        xvmsh(idime) = xvmsh(idime) + xshai*elmsh(idime,inode)                      ! only different than zero when coupled to alefor
     end do

     

     !
     ! Turbulent viscosity at gauss point when coupled with TURBUL
     !
     if (kfl_cotur_nsa /= 0 ) then 
        ipoin = lnods(inode,ielem)
        dwall = dwall + xshai * walld(ipoin)    ! Interpolation of the wall distance at gauss point
        if (kfl_cotur_nsa == 1) then
           xnutu = xnutu + xshai * turmu(ipoin)
        endif
     endif
     !
     ! Molecular weight at gauss point when coupled with CHEMIC
     !
     if (kfl_coupl(ID_NASTAL,ID_CHEMIC) >= 1 ) then
        ipoin = lnods(inode,ielem)
        xmowe = xmowe + xshai * wmean(ipoin,1)
     endif

  end do

  !
  ! Mesh velocity
  !     
  do inode = 1,pnode
     xshai= xshap(inode)
     ipoin = lnods(inode,ielem)
     do idime = 1,ndime
     end do
  end do

  do idofn=1,ndofn_nsa
     xtunk(idofn) = xunkn(idofn,igaus,ITER_K)
     do idime=1,ndime           
        gtunk(idofn,idime) = gunkn(idofn,idime)
     end do
  end do

  if (kfl_track_nsa == 1) then
     do idofn=1,ndofn_nsa
        xtunk(idofn) = xunkn(idofn,igaus,ITER_K) + xsube(idofn,igaus,2)
        do idime=1,ndime           
           gtunk(idofn,idime) = gunkn(idofn,idime) + gsube(idofn,idime)
        end do
     end do
  end if

  velsq= 0.0_rp
  xpres= 0.0_rp
  velmo= 0.0_rp
  do idime=1,ndime
     xvelo(idime) =  xtunk(idime) / xtunk(ndime+1)
     velsq = velsq + xvelo(idime) * xvelo(idime)
     xpres = xpres + xtunk(idime) * xtunk(idime)
     do jdime= 1,ndime
        gvelo(idime,jdime) = (gtunk(idime,jdime) &
             - xtunk(idime) * gtunk(ndime+1,jdime) / xtunk(ndime+1)) &
             / xtunk(ndime+1)
     end do
  end do

  if (kfl_prope /= 0 ) then
     call ker_proper('SPHEA','IGAUS',1_ip,ielem,prope_tmp,pnode,pgaus,xshap,cartd) 
     xhecp = prope_tmp(1)
     if (kfl_coupl(ID_NASTAL,ID_CHEMIC) == 0 ) then
        xmowe = mowei_nsa
     endif
     rgasc = runiv_nsa / xmowe  
     xhecv = xhecp - rgasc  !Cv is computed from R & Cp
  else
     xhecp = cpcoe_nsa 
     xmowe = mowei_nsa
     rgasc = runiv_nsa / xmowe
     xhecv = xhecp - rgasc  !Cv is computed from R & Cp
  endif

  do idime=1,ndime
     do jdime=1,ndime
        gpres(idime) = gpres(idime) + gvelo(jdime,idime) * xtunk(jdime) &
             + xvelo(jdime) * gtunk(jdime,idime)
     end do
     gpres(idime) = rgasc * (gtunk(ndime+2,idime) &
          - 0.5_rp * gpres(idime)) / xhecv 
  end do

  velmo = sqrt(velsq)
  xpres = rgasc * (xtunk(ndime+2) - 0.5_rp * xpres &
       / xtunk(ndime+1))  /   xhecv

  if (xpres .lt. zensa) then
     xpres = xpreo
     do idime=1,ndime
        gpres(idime)= gpreo(idime)
     end do
  end if

  xtemp = xpres / xtunk(ndime+1) / rgasc

  dvite = 0.0_rp
  xvisc = 0.0_rp
  xdith = 0.0_rp
  sgsdi = 0.0_rp

  if (kfl_visco_nsa > 0) then 
     !
     ! SGS viscous dissipation for LES
     !
     if (kfl_cotur_nsa < 0) then     
        seci4 = 0.0_rp
        do idime = 1,ndime                     ! 2 S_ij : S_ij
           do jdime = 1,ndime         
              seci4 = seci4 + gvelo(idime,jdime) * (gvelo(idime,jdime) + gvelo(jdime,idime))
           end do
        end do
        xmile =  dvolu**0.3333333_rp
        !
        ! SGS_DISSIPATION = C_eps * (K^sgs)**3/2 / V**1/3, K^sgs = sqrt(3/4) * nut * |S|
        !
        sgsdi = 0.916_rp * 0.866025_rp * (sqrt(seci4)**1.5_rp) / xmile 
     endif
     !
     ! Properties from the kernel: viscosity mu, c_p, K
     !
     if (kfl_prope /= 0 ) then
        call ker_proper('VISCO','IGAUS',1_ip,ielem,prope_tmp,pnode,pgaus,xshap,cartd) 
        xvisc = prope_tmp(1)
        call ker_proper('CONDU','IGAUS',1_ip,ielem,prope_tmp,pnode,pgaus,xshap,cartd) 
        xdith = prope_tmp(1)
        !
        ! Computation viscosity derivative dmu/dT
        !
        if (lawvi_nsa == 1) then               ! power law
           dvite = vispa_nsa(1)*vispa_nsa(2)*xtemp**(vispa_nsa(2)-1.0_rp)        
        else if (lawvi_nsa==2) then            ! sutherland law
           dvite = xvisc*( 1.5_rp/xtemp - 1.0_rp/(xtemp+vispa_nsa(2)))
        end if

        if (kfl_cotur_nsa /= 0) then           ! If turbulence model ON 
           if (kfl_cotur_nsa < 0) then        ! Turbulent viscosity 
              call nsa_turbul(dwall,velmo,xvisc,gvelo,dvolu,xnutu,xunkn(1,igaus,ITER_K))
              sgsdi = sgsdi*xnutu
              xnutu = xnutu*xtunk(ndime+1)
           endif
           xvisc = xvisc + xnutu              ! If turbulence model ON, xnutu /= 0, otherwise = 0 
           xdith = xdith + xnutu  * xhecp / prand_nsa
        endif

     else
        call nsa_lawvis(-1,1,xvisc,xtemp,dvite)! Dynamic viscosity 

        if (kfl_cotur_nsa < 0) then            ! LES Turbulent viscosity 
           call nsa_turbul(dwall,velmo,xvisc,gvelo,dvolu,xnutu,xunkn(1,igaus,ITER_K))
           sgsdi = sgsdi*xnutu
           xnutu = xnutu*xtunk(ndime+1)
        endif
        xvisc = xvisc + xnutu                  ! If turbulence model ON, xnutu /= 0, otherwise = 0 
        xdith = xvisc * xhecp / prand_nsa
     endif
  endif

  dvelo = 0.0_rp

  do idime=1,ndime
     gtemp(idime) = (gpres(idime) - xpres * gtunk(ndime+1,idime) &
          / xtunk(ndime+1)) &
          / xtunk(ndime+1) / rgasc
     gvisc(idime) = dvite * gtemp(idime)
     dvelo = dvelo + gvelo(idime,idime)
  end do
  enepe = xtunk(ndime+2) / xtunk(ndime+1)

  visci = xvisc / xtunk(ndime+1)
  dicod = xdith / xhecv / xtunk(ndime+1)

  rgacv = rgasc / xhecv
  adgam = xhecp / xhecv

  sound= sqrt(adgam * rgasc * xtemp)



  ! xvmsh is the mesh velocity, only different than zero when coupled to alefor

  do jdime=1,ndime        
     xconv(ndime+1,jdime  ,jdime)= 1.0_rp
     xconv(jdime  ,ndime+2,jdime)= rgacv
     xconv(jdime  ,ndime+1,jdime)= rgacv * 0.5_rp * velsq
     xconv(ndime+2,jdime  ,jdime)= ((1.0_rp + rgacv) * enepe - rgacv * 0.5_rp * velsq)
     xconv(ndime+2,ndime+1,jdime)= - xvelo(jdime) * ((1.0_rp + rgacv) * enepe - rgacv * velsq) + htrad(jdime)
     xconv(ndime+2,ndime+2,jdime)= ((1.0_rp + rgacv) * ( xvelo(jdime) - xvmsh(jdime) ))
     xdiff(ndime+2,ndime+1,jdime,jdime)= (dicod-visci) * velsq - dicod * enepe
     xdiff(ndime+2,ndime+2,jdime,jdime)= dicod
     do idime=1,ndime
        xconv(idime,idime  ,jdime)= xconv(idime,idime  ,jdime) + ( xvelo(jdime) - xvmsh(jdime) ) 
        xconv(jdime,idime  ,jdime)= xconv(jdime,idime  ,jdime) - rgacv * xvelo(idime)
        xconv(idime,jdime  ,jdime)= xconv(idime,jdime  ,jdime) + xvelo(idime)
        xconv(idime,ndime+1,jdime)= xconv(idime,ndime+1,jdime) - xvelo(idime) * xvelo(jdime)
        xconv(ndime+2,idime,jdime)= xconv(ndime+2,idime,jdime) - rgacv * xvelo(idime) * xvelo(jdime)
        xdiff(jdime,jdime,idime,idime)= visci
        xdiff(jdime,idime,idime,jdime)= xdiff(jdime,idime,idime,jdime) + visci
        xdiff(jdime,idime,jdime,idime)= xdiff(jdime,idime,jdime,idime) - 2.0_rp * visci / 3.0_rp
        xdiff(jdime,ndime+1,idime,idime)= - visci * xvelo(jdime)
        xdiff(jdime,ndime+1,idime,jdime)= xdiff(jdime,ndime+1,idime,jdime) - visci * xvelo(idime)
        xdiff(jdime,ndime+1,jdime,idime)= xdiff(jdime,ndime+1,jdime,idime) + 2.0_rp * visci * xvelo(idime) / 3.0_rp
        xdiff(ndime+2,idime,jdime,jdime)= (visci-dicod) * xvelo(idime)
        xdiff(ndime+2,jdime,jdime,idime)= xdiff(ndime+2,jdime,jdime,idime) + visci * xvelo(idime)
        xdiff(ndime+2,jdime,idime,jdime)= xdiff(ndime+2,jdime,idime,jdime) - 2.0_rp * visci * xvelo(idime) / 3.0_rp
        xdiff(ndime+2,ndime+1,jdime,idime)= xdiff(ndime+2,ndime+1,jdime,idime) + &
             0.5_rp * visci * xvelo(jdime) * xvelo(idime)
     end do
  end do

  ! skew symmetric
  xconv= xconv * factor_skew
  dflux_conv = dflux_conv * 0.5_rp

  !  do idime=1,ndime
  !     xconv_der(ndime+2,ndime+1,ndime+2,idime) = xvelo(idime)/xtunk(ndime+1)*(1.0_rp-rgacv)
  !     xconv_der(ndime+2,idime,ndime+2,idime) = (1.0_rp+rgacv)/xtunk(ndime+1)
  !     xconv_der(ndime+2,ndime+2,ndime+1,idime) = (1.0_rp+rgacv)/xtunk(ndime+1)*xvelo(idime)
  !     do kdime=1,ndime
  !        xconv_der(idime,kdime,idime,kdime) = (1/xtunk(ndime+1))
  !        xconv_der(idime,kdime,kdime,idime) =-rgacv/xtunk(ndime+1)
  !        xconv_der(idime,ndime+1,idime,kdime) = -xvelo(kdime)/xtunk(ndime+1)
  !        xconv_der(idime,ndime+1,kdime,idime) = &
  !             xconv_der(idime,ndime+1,idime,kdime) + rgacv/(xtunk(ndime+1)*xtunk(ndime+1))*xvelo(kdime)
  !        xconv_der(ndime+2,idime,idime,kdime) = -xvelo(kdime)/xtunk(ndime+1)*rgacv
  !        xconv_der(ndime+1,idime,kdime,idime) = &
  !             xconv_der(ndime+1,idime,idime,kdime) - xvelo(kdime)/xtunk(ndime+1)*rgacv
  !        xconv_der(ndime+2,ndime+1,kdime,idime) = 2.0_rp* rgacv/xtunk(ndime+1)*xvelo(idime)*xvelo(kdime)
  !        xconv_der(kdime,idime,ndime+1,idime) = -xvelo(kdime)/xtunk(ndime+1)
  !        xconv_der(ndime+2,kdime,ndime+1,idime) = 2.0_rp* rgacv/xtunk(ndime+1)*xvelo(idime)*xvelo(kdime)
  !        xconv_der(kdime,ndime+2,ndime+1,idime) = 2.0_rp/xtunk(ndime+1)*xvelo(idime)*xvelo(kdime)
  !        xconv_der(idime,kdime,ndime+1,idime) = rgacv/xtunk(ndime+1)*xvelo(kdime)
  !        xconv_der(kdime,kdime,ndime+1,idime) = -xvelo(idime)/xtunk(ndime+1)
  !     end do
  !     xconv_der(idime,idime,idime,idime) = xconv_der(idime,idime,idime,idime)+1/xtunk(ndime+1)
  !     xconv_der(ndime+2,ndime+1,idime,idime)= &
  !          xconv_der(ndime+2,ndime+1,idime,idime) &
  !          + 1/xtunk(ndime+1)*(xtunk(ndime+2)*(-1.0_rp-rgacv)/xtunk(ndime+1)+rgacv*velsq)
  !     xconv_der(ndime+2,ndime+2,idime,idime)=(1.0_rp+rgacv)/xtunk(ndime+1)
  !     xconv_der(ndime+2,idime,idime,idime) = xconv_der(ndime+2,idime,idime,idime) -rgacv*xvelo(idime)/xtunk(idime+1)
  !     xconv_der(idime,idime,ndime+1,idime) = xconv_der(idime,idime,idime,ndime+1)*(2.0_rp-rgacv)
  !     xconv_der(ndime+2,idime,ndime+1,idime) = xconv_der(ndime+2,idime,idime,ndime+1) &
  !          + velsq*rgacv/(xtunk(ndime+1))-xtunk(ndime+2)*(1.0_rp+rgacv)/(xtunk(ndime+1)*xtunk(ndime+1))
!!!!! ojoooo     xconv_der(idime,ndime+1,ndime+1,idime) = xconv_der(idime,ndime+1,idime,ndime+1) - rgacv*(velsq)/(*xtunk(ndime+1))
  !     xconv_der(ndime+2,ndime+1,ndime+1,idime) = &
  !          2.0_rp*xvelo(idime)*(xtunk(ndime+2)*(1.0_rp+rgacv)/(xtunk(ndime+1)*xtunk(ndime+1))-rgacv*(velmo*velmo)/xtunk(ndime+1))
  !     do kdime = 1,ndime
  !        xconv_der(idime,idime,kdime,kdime)= xconv_der(idime,idime,kdime,kdime)+1/xtunk(ndime+1)
  !        xconv_der(idime,ndime+1,kdime,kdime)= xconv_der(idime,ndime+1,kdime,kdime)-xvelo(idime)/(xtunk(ndime+1))
  !     end do
  !  end do


  xmach = 1.0e-10_rp          ! very small mach number when u=0 to avoid ill definition of xmach
  if (velmo > 0.0_rp) then
     xmach = velmo / sound
  end if

  !
  ! A LOCAL PRECONDITIONER IS APPLIED
  !
  if (kfl_lopre_nsa > 0) then

     do idofn=1,ndofn_nsa
        do jdofn=1,ndofn_nsa
           do idime=1,ndime
              xconv_ori(idofn,jdofn,idime) = xconv(idofn,jdofn,idime)
              xconv(idofn,jdofn,idime) = 0.0_rp
              do jdime=1,ndime
                 xdiff_ori(idofn,jdofn,idime,jdime) = xdiff(idofn,jdofn,idime,jdime)
                 xdiff(idofn,jdofn,idime,jdime) = 0.0_rp
              end do
           end do
        end do
     end do

     if (kfl_lopre_nsa == 1) then ! IDENTITY LOCAL PRECONDITIONER (JUST FOR TESTING PURPOUSES)

        do idofn=1,ndofn_nsa
           xlopr_conservative(idofn,idofn) = 1.0_rp           
        end do

        do idofn=1,ndofn_nsa
           do jdofn=1,ndofn_nsa
              if (kfl_pseud_nsa == 1) &
                   xtide(idofn) = xtide(idofn) &
                   + xlopr_conservative(idofn,jdofn) * (xunkn(jdofn,igaus,ITER_K)-xunkn(jdofn,igaus,TIME_N))
           end do
           if (kfl_pseud_nsa == 1) xtide(idofn) = xtide(idofn) * dtinv
        end do

        do idofn=1,ndofn_nsa
           do jdofn=1,ndofn_nsa
              do idime=1,ndime
                 xconv(idofn,jdofn,idime) = xconv_ori(idofn,jdofn,idime)
                 do jdime=1,ndime
                    xdiff(idofn,jdofn,idime,jdime) = xdiff_ori(idofn,jdofn,idime,jdime)
                 end do
              end do
           end do
        end do




     else if (kfl_lopre_nsa == 2) then  ! VAN LEER-LEE-ROE STEADY INVISCID LOCAL PRECONDITIONER IS APPLIED

        if (velsq == 0.0_rp) then

           ! when velsq is strictly zero and VLR is used, do not precondition

           xlopr_conservative= 1.0_rp           

           do idofn=1,ndofn_nsa
              do jdofn=1,ndofn_nsa
                 if (kfl_pseud_nsa == 1) &
                      xtide(idofn) = xtide(idofn) &
                      + xlopr_conservative(idofn,jdofn) * (xunkn(jdofn,igaus,ITER_K)-xunkn(jdofn,igaus,TIME_N))
              end do
              if (kfl_pseud_nsa == 1) xtide(idofn) = xtide(idofn) * dtinv
           end do

           xconv = xconv_ori
           xdiff = xdiff_ori

        else

           if (xmach < 1.0_rp-zensa_lopre) then
              xbeta = sqrt(1.0_rp-xmach*xmach)
              xtaup = xbeta
           else if (xmach < 1.0_rp) then
              xbeta = sqrt((2.0_rp-zensa_lopre)*zensa_lopre)
              xtaup = xbeta
           else if (xmach < 1.0_rp+zensa_lopre) then
              xbeta = sqrt((2.0_rp+zensa_lopre)*zensa_lopre)
              xtaup = xbeta / xmach
           else
              xbeta = sqrt(xmach*xmach-1.0_rp)
              xtaup = xbeta / xmach
           end if

!!! BEGIN: VLR preconditioner with conservative variables (any change of variables is applied)
           termA = xtaup / xbeta / xbeta
           termB = 1.0_rp - xmach * xmach
           termC = 0.5_rp * rgacv * xmach * xmach
           termD = 1.0_rp / sound / sound
           termE = termA * termD
           termF = rgacv * termD
           termG = termA * xmach * xmach
           termH = rgacv * termD * (1.0_rp - termG)
           termI = termC * (termG - 1.0_rp)

           term1 = (1.0_rp+termA-xtaup) / velsq + (rgacv*termA*termB-termA+rgacv) * termD
           term2 = - termA * termB * (1.0_rp+termC) - termC
           term3 = - rgacv * (termA*termB+1.0_rp) * termD
           term4 = - termE + termH
           term5 = 1.0_rp + termA -termA/rgacv - 3.0_rp*termG/2.0_rp + rgacv*termG + termC*(1.0_rp-termG)

           do idofn=1,ndime
              xlopr_conservative(idofn,ndime+1) = term2 * xvelo(idofn)
              xlopr_conservative(idofn,ndime+2) = term3 * xvelo(idofn)
              xlopr_conservative(ndime+1,idofn) = term4 * xvelo(idofn)
              xlopr_conservative(ndime+2,idofn) = term5 * xvelo(idofn)
              do jdofn=idofn,ndime
                 xlopr_conservative(idofn,jdofn) = term1 * xvelo(idofn) * xvelo(jdofn)
                 xlopr_conservative(jdofn,idofn) = xlopr_conservative(idofn,jdofn)
              end do
              xlopr_conservative(idofn,idofn) = xlopr_conservative(idofn,idofn) + xtaup 
           end do
           xlopr_conservative(ndime+1,ndime+1) = 1.0_rp + termG + termC * (termG-1.0_rp)
           xlopr_conservative(ndime+1,ndime+2) = - termH
           xlopr_conservative(ndime+2,ndime+1) = (1.0_rp/rgacv-termB) * termA * velsq - 0.5_rp * velsq &
                + 0.5_rp * termI * velsq - termC * velsq * termA
           xlopr_conservative(ndime+2,ndime+2) = (1.0_rp-rgacv) * termG + termI
!!! END: VLR preconditioner with conservative variables (any change of variables is applied)

           ! Compute P_conservative * A
           do idofn=1,ndofn_nsa
              do jdofn=1,ndofn_nsa
                 if (kfl_pseud_nsa == 1) &
                      xtide(idofn) = xtide(idofn) &
                      + xlopr_conservative(idofn,jdofn) * (xunkn(jdofn,igaus,ITER_K)-xunkn(jdofn,igaus,TIME_N))
                 do idime=1,ndime
                    do pdofn=1,ndofn_nsa
                       xconv(idofn,jdofn,idime) = &
                            xconv(idofn,jdofn,idime) + xlopr_conservative(idofn,pdofn)*xconv_ori(pdofn,jdofn,idime)
                    end do
                 end do
              end do
              if (kfl_pseud_nsa == 1) xtide(idofn) = xtide(idofn) * dtinv
           end do

        end if

     else if (kfl_lopre_nsa == 3) then  ! CHOI & MERKLE STEADY/UNSTEADY VISCOUS/INVISCID PRECONDITIONER IS APPLIED

        do idime=1,ndime
           gvemo(idime) = 0.0_rp
           greyn_cell_inv(idime) = 0.0_rp
           if (velmo > 0.0_rp) then
              gvemo(idime) = gvemo(idime) / velmo
              do jdime=1,ndime
                 gvemo(idime) = gvemo(idime) + xvelo(jdime) * gvelo(jdime,idime)
              end do
              !              greyn_cell_inv(idime) = - xvisc * (gtunk(ndime+1,idime)*velmo+xtunk(ndime+1)*gvemo(idime)) &
              !                   / hmaxi / xtunk(ndime+1) / xtunk(ndime+1) / velsq
              greyn_cell_inv(idime) = - xvisc * (gtunk(ndime+1,idime)*velmo+xtunk(ndime+1)*gvemo(idime)) &
                   / hmini / xtunk(ndime+1) / xtunk(ndime+1) / velsq
           end if

           gsoun(idime) = 0.5_rp * adgam * rgasc * gtemp(idime) / sound
           gmach(idime) = (gvemo(idime) - xmach * gsoun(idime)) / sound 
           gmach_ref(idime) = gmach(idime)

        end do

        xmach_ref = xmach  ! C&M STEADY CASE
        if (kfl_pseud_nsa == 1) then  ! C&M UNSTEADY CASE        
           cfl_sound = sound / dtinv / hmini
           xmach_ref = sqrt(xmach * xmach + 1.0_rp / cfl_sound / cfl_sound)
        end if

        if (xmach_ref < 0.00001_rp) then 
           xmach_ref = 0.00001_rp
           gmach_ref = 0.0_rp
        else if (xmach_ref > 1.0_rp) then
           xmach_ref = 1.0_rp
           gmach_ref = 0.0_rp
        end if

        if (kfl_visco_nsa == 0) then  ! C&M INVISCID CASE
           xbeta_visc = 1.0_rp
           gbeta_visc = 0.0_rp
        else  ! C&M VISCOUS CASE
           !           reyno_cell     = xtunk(ndime+1) * velmo * hmaxi / xvisc
           reyno_cell     = xtunk(ndime+1) * velmo * hmini / xvisc
           reyno_cell_inv = 1.0e10_rp   ! very large 1/Re when u=0
           if (reyno_cell > 0.0_rp) reyno_cell_inv = 1.0_rp / reyno_cell
           termA = xmach_ref * xmach_ref * (reyno_cell_inv - 1.0_rp + 1.0_rp / xmach / xmach)
           xbeta_visc = reyno_cell_inv * (reyno_cell_inv - 1.0_rp) / termA
           do idime=1,ndime
              termB = 2.0_rp * xmach_ref * gmach_ref(idime) * (reyno_cell_inv - 1.0_rp + 1.0_rp / xmach / xmach) &
                   + xmach_ref * xmach_ref * (greyn_cell_inv(idime) - 2.0_rp * gmach(idime) / xmach / xmach / xmach)
              gbeta_visc(idime) = (greyn_cell_inv(idime) * (2.0_rp * reyno_cell_inv - 1.0_rp) &
                   - reyno_cell_inv * (reyno_cell_inv - 1.0_rp) * termB / termA) / termA
           end do
           if (xbeta_visc < 1.0_rp) then
              xbeta_visc = 1.0_rp
              gbeta_visc = 0.0_rp
           end if
        end if

        xbeta = xbeta_visc * sound * sound
        xtaup = 1.0_rp
        do idime=1,ndime
           gbeta(idime) = gbeta_visc(idime) * sound * sound + 2.0_rp * sound * xbeta_visc * gsoun(idime)
        end do

!!! BEGIN: CHOI & MERKLE preconditioner with primitive variables (p, u T)
        termA = xbeta * xmach_ref * xmach_ref
        termB = 0.5_rp * velsq - xhecp * xtemp + termA * xtaup

        do idime=1,ndime
           xlopr(idime,idime) = 1.0_rp / xtunk(ndime+1)
           xlopr(idime,ndime+1) = - xvelo(idime) / xtunk(ndime+1)
           xlopr(ndime+2,idime) = - xvelo(idime) / xtunk(ndime+1) / xhecp
        end do
        xlopr(ndime+1,ndime+1) = termA
        xlopr(ndime+2,ndime+1) = termB / xtunk(ndime+1) / xhecp
        xlopr(ndime+2,ndime+2) = 1.0_rp / xtunk(ndime+1) / xhecp
!!! END: CHOI & MERKLE preconditioner with primitive variables (p, u T)

!!! BEGIN: CHOI & MERKLE preconditioner with conservative variables (we apply a change of variables)
        do idime=1,ndime
           xftra(idime,idime) = xtunk(ndime+1)
           xftra(idime,ndime+1) = xtunk(idime) / xpres
           xftra(idime,ndime+2) = - xtunk(idime) / xtemp
           xftra(ndime+2,idime) = xtunk(idime)
        end do
        xftra(ndime+1,ndime+1) = xtunk(ndime+1) / xpres
        xftra(ndime+1,ndime+2) = - xtunk(ndime+1) / xtemp
        xftra(ndime+2,ndime+1) = 0.5_rp * velsq / rgasc / xtemp + 1.0_rp / rgacv
        xftra(ndime+2,ndime+2) = - 0.5_rp * xtunk(ndime+1) * velsq / xtemp

        do idime=1,ndime
           xftra_inv(idime,idime) = 1.0_rp / xtunk(ndime+1)
           xftra_inv(idime,ndime+1) = - xvelo(idime) / xtunk(ndime+1)
           xftra_inv(ndime+1,idime) = - rgacv * xvelo(idime)
           xftra_inv(ndime+2,idime) = - xvelo(idime) / xhecv / xtunk(ndime+1)
        end do
        xftra_inv(ndime+1,ndime+1) = 0.5_rp * rgacv * velsq
        xftra_inv(ndime+1,ndime+2) = rgacv
        xftra_inv(ndime+2,ndime+1) = (velsq - xtunk(ndime+2) / xtunk(ndime+1)) / xhecv / xtunk(ndime+1)
        xftra_inv(ndime+2,ndime+2) = 1.0_rp / xhecv / xtunk(ndime+1)

!!$        do idofn=1,ndofn_nsa
!!$           do jdofn=1,ndofn_nsa
!!$              do ldofn=1,ndofn_nsa
!!$                 do mdofn=1,ndofn_nsa
!!$                    xlopr_conservative(idofn,jdofn) = xlopr_conservative(idofn,jdofn) + xftra_inv(idofn,ldofn) &
!!$                         *xlopr(ldofn,mdofn)*xftra(mdofn,jdofn)
!!$                 end do
!!$              end do
!!$           end do
!!$        end do


        do idofn=1,ndofn_nsa
           do jdofn=1,ndofn_nsa
              do ldofn=1,ndofn_nsa
                 xlopr_conservative(idofn,jdofn) = xlopr_conservative(idofn,jdofn) + xftra(idofn,ldofn) &
                      *xlopr(ldofn,jdofn)
              end do
           end do
        end do

!!! END: CHOI & MERKLE preconditioner with conservative variables (we apply a change of variables) 

!!! BEGIN: dP_conservative computation
        termB = 1.0_rp + 0.5_rp * velsq / xtemp / xhecp
        termC = xbeta * xmach_ref * xmach_ref / rgacv - 0.5_rp * velsq
        termD = 1.0_rp / xtemp / xhecp
        do kdime=1,ndime
           term1 = 0.0_rp
           do idime=1,ndime
              term1 = term1 + xvelo(idime) * gvelo(idime,kdime)
           end do
           termB_der = (term1 - 0.5_rp * velsq * gtemp(kdime) / xtemp) / xtemp / xhecp
           termC_der = (gbeta(kdime) * xmach_ref * xmach_ref + 2.0_rp * xbeta * xmach_ref * gmach_ref(kdime)) / rgacv - term1
           termD_der = - gtemp(kdime) / xtemp / xtemp / xhecp
           do idime=1,ndime
              termA = xvelo(idime) / xtemp / xhecp
              termA_der = (gvelo(idime,kdime)-xvelo(idime)*gtemp(kdime)/xtemp) / xtemp / xhecp
              do jdime=1,ndime
                 glopr(idime,jdime,kdime) = ( ( xvelo(idime)*gvelo(jdime,kdime)+xvelo(jdime)*gvelo(idime,kdime) ) &
                      - xvelo(idime)*xvelo(jdime)*gtemp(kdime)/xtemp ) / xtemp / xhecp
              end do
              glopr(ndime+1,idime,kdime) = termA_der
              glopr(ndime+2,idime,kdime) = gvelo(idime,kdime) * termB + xvelo(idime) * termB_der
              glopr(idime,ndime+1,kdime) = termA_der * termC + termA * termC_der
              glopr(idime,ndime+2,kdime) = - termA_der
           end do
           glopr(ndime+1,ndime+1,kdime) = termD_der * termC + termD * termC_der
           glopr(ndime+1,ndime+2,kdime) = - termD_der
           glopr(ndime+2,ndime+1,kdime) = termC_der * termB + termC * termB_der
           glopr(ndime+2,ndime+2,kdime) = - termB_der
        end do
!!! END: dP_conservative computation

        ! Compute: P_conservative * A + dP_conservative * K
        !          P_conservative * K

        do idofn=1,ndofn_nsa
           do jdofn=1,ndofn_nsa
              do jdime=1,ndime
                 do pdofn=1,ndofn_nsa
                    do idime=1,ndime
                       xconv_lopre_diff(idofn,jdofn,jdime) = xconv_lopre_diff(idofn,jdofn,jdime) + &
                            glopr(idofn,pdofn,idime) * xdiff_ori(pdofn,jdofn,idime,jdime)
                    end do
                 end do
              end do
           end do
        end do

        do idofn=1,ndofn_nsa
           do jdofn=1,ndofn_nsa
              if (kfl_pseud_nsa == 1) &
                   xtide(idofn) = xtide(idofn) &
                   + xlopr_conservative(idofn,jdofn) * (xunkn(jdofn,igaus,ITER_K)-xunkn(jdofn,igaus,TIME_N))
              do idime=1,ndime
                 do pdofn=1,ndofn_nsa
                    xconv(idofn,jdofn,idime) = &
                         xconv(idofn,jdofn,idime) + xlopr_conservative(idofn,pdofn)*xconv_ori(pdofn,jdofn,idime)
                    do jdime=1,ndime
                       xdiff(idofn,jdofn,idime,jdime) = xdiff(idofn,jdofn,idime,jdime) + &
                            xlopr_conservative(idofn,pdofn)*xdiff_ori(pdofn,jdofn,idime,jdime)
                    end do
                 end do
                 xconv(idofn,jdofn,idime) = xconv(idofn,jdofn,idime) + xconv_lopre_diff(idofn,jdofn,idime)
              end do
           end do
           if (kfl_pseud_nsa == 1) xtide(idofn) = xtide(idofn) * dtinv
        end do


!!$        ! Compute: P_conservative * A
!!$        !          P_conservative * K
!!$        do idofn=1,ndofn_nsa
!!$           do jdofn=1,ndofn_nsa
!!$              if (kfl_pseud_nsa == 1) &
!!$                   xtide(idofn) = xtide(idofn) &
!!$                   + xlopr_conservative(idofn,jdofn) * (xunkn(jdofn,igaus,ITER_K)-xunkn(jdofn,igaus,TIME_N))
!!$              do idime=1,ndime
!!$                 do pdofn=1,ndofn_nsa
!!$                    xconv(idofn,jdofn,idime) = &
!!$                         xconv(idofn,jdofn,idime) + xlopr_conservative(idofn,pdofn)*xconv_ori(pdofn,jdofn,idime)
!!$                    do jdime=1,ndime
!!$                       xdiff(idofn,jdofn,idime,jdime) = xdiff(idofn,jdofn,idime,jdime) + &
!!$                            xlopr_conservative(idofn,pdofn)*xdiff_ori(pdofn,jdofn,idime,jdime)
!!$                    end do
!!$                 end do
!!$              end do
!!$           end do
!!$           if (kfl_pseud_nsa == 1) xtide(idofn) = xtide(idofn) * dtinv
!!$        end do

     else if (kfl_lopre_nsa == 4) then  ! CHOI & MERKLE INVISCID UNSTEADY PRECONDITIONER IS APPLIED

        xmach_ref = xmach
        if (xmach_ref < 0.0001_rp) then 
           xmach_ref = 0.0001_rp
        else if (xmach_ref > 1.0_rp) then
           xmach_ref = 1.0_rp
        end if

        xbeta = sound * sound
        xtaup = 1.0_rp

!!! BEGIN: CHOI & MERKLE INVISCID UNSTEADY preconditioner with conservative variables (any change of variables is applied)
        termA = 1.0_rp / xhecp / xtemp
        termB = xbeta * xmach_ref * xmach_ref
        termC = 0.5_rp * velsq - xhecp * xtemp + termB * xtaup
        termD = 1.0_rp / rgasc / xtemp
        termE = 1.0_rp + 0.5_rp * termA * velsq

        do idime=1,ndofn_nsa
           do jdime=1,ndime
              xlopr_conservative(idime,jdime) = xvelo(idime) * xvelo(jdime) * termA
           end do
           xlopr_conservative(idime,idime) = xlopr_conservative(idime,idime) + 1.0_rp
           xlopr_conservative(idime,ndime+1) = xvelo(idime) * (termB*termD-1.0_rp-termA*termC)
           xlopr_conservative(idime,ndime+2) = - xvelo(idofn) * termA
           xlopr_conservative(ndime+1,idime) = xvelo(idime) * termA
           xlopr_conservative(ndime+2,idime) = xvelo(idime) * termE
        end do
        xlopr_conservative(ndime+1,ndime+1) = termB * termD - termA * termC
        xlopr_conservative(ndime+1,ndime+2) = - termA
        xlopr_conservative(ndime+2,ndime+1) = 1.0_rp / rgacv + 0.5_rp * adgam * xmach * xmach &
             - velsq * (1.0_rp + 0.5_rp * termA * termC)
        xlopr_conservative(ndime+2,ndime+2) = 0.5_rp * velsq * termA
!!! END: CHOI & MERKLE INVISCID UNSTEADY preconditioner with conservative variables (any change of variables is applied)

        ! Compute P_conservative * A
        do idofn=1,ndofn_nsa
           do jdofn=1,ndofn_nsa
              if (kfl_pseud_nsa == 1) xtide(idofn) = xtide(idofn) + xlopr_conservative(idofn,jdofn) * (xunkn(jdofn,igaus,ITER_K)-xunkn(jdofn,igaus,TIME_N))
              do idime=1,ndime
                 do pdofn=1,ndofn_nsa
                    xconv(idofn,jdofn,idime) = xconv(idofn,jdofn,idime) + xlopr_conservative(idofn,pdofn)*xconv_ori(pdofn,jdofn,idime)
                 end do
              end do

           end do
           if (kfl_pseud_nsa == 1) xtide(idofn) = xtide(idofn) * dtinv
        end do



     else ! OTHER LOCAL PRECONDITIONERS ARE APPLIED

        ! to program!! other local preconditioners!!! (MM)

        velmo_xy = sqrt(xvelo(1)*xvelo(1)+xvelo(2)*xvelo(2))
!!$     if (velmo < zensa) then           
!!$        velmo = zensa 
!!$     end if
!!$     if (velmo_xy < zensa) then           
!!$        velmo_xy = zensa 
!!$     end if

        if (xmach < 1.0_rp-zensa_lopre) then
           xbeta = sqrt(1.0_rp-xmach*xmach)
           xtaup = xbeta
        else if (xmach < 1.0_rp) then
           xbeta = sqrt((2.0_rp-zensa_lopre)*zensa_lopre)
           xtaup = xbeta
        else if (xmach < 1.0_rp+zensa_lopre) then
           xbeta = sqrt((2.0_rp+zensa_lopre)*zensa_lopre)
           xtaup = xbeta / xmach
        else
           xbeta = sqrt(xmach*xmach-1.0_rp)
           xtaup = xbeta / xmach
        end if

        ! program another local preconditioner HERE!!!!
!!! BEGIN: VLR preconditioner with symmetrizing variables and streamline coordinates
        xlopr(1,1) = 1.0_rp + xtaup / xbeta / xbeta
        do idime=2,ndime
           xlopr(idime,idime) = xtaup
        end do
        xlopr(1,ndime+1) = - xtaup * xmach / xbeta / xbeta
        xlopr(ndime+1,1) = xlopr(1,ndime+1)
        xlopr(ndime+1,ndime+1) = xtaup * xmach * xmach / xbeta / xbeta
        xlopr(ndime+2,ndime+2) = 1.0_rp
!!! END: VLR preconditioner with symmetrizing variables and streamline coordinates

!!! BEGIN: VLR preconditioner with conservative variables (we apply a change of variables)
!!$        do idofn=1,ndime
!!$           xftra(1,idofn) = xvelo(idofn) / xtunk(ndime+1) / velmo
!!$           xftra(ndime+1,idofn) = - rgacv * xvelo(idofn) / xtunk(ndime+1) / sound
!!$           xftra(ndime+2,idofn) = - rgasc * xvelo(idofn) / xpres
!!$        end do
!!$        xftra(1,ndime+1) = - velmo / xtunk(ndime+1)
!!$        xftra(2,1) = - xvelo(2) / xtunk(ndime+1) / velmo_xy
!!$        xftra(2,2) = xvelo(1) / xtunk(ndime+1) / velmo_xy
!!$        if (ndime==3) then
!!$           xftra(3,1) = - xvelo(1) * xvelo(3) /  xtunk(ndime+1) / velmo / velmo_xy
!!$           xftra(3,2) = - xvelo(2) * xvelo(3) /  xtunk(ndime+1) / velmo / velmo_xy
!!$           xftra(3,3) = velmo_xy /  xtunk(ndime+1) / velmo
!!$        end if
!!$        xftra(ndime+1,ndime+1) = 0.5_rp * rgacv * sound * xmach * xmach / xtunk(ndime+1)
!!$        xftra(ndime+2,ndime+1) = - xhecv * sound * sound / xpres + 0.5_rp * rgasc * sound * sound * xmach * xmach / xpres
!!$        xftra(ndime+1,ndime+2) = rgacv / xtunk(ndime+1) / sound 
!!$        xftra(ndime+2,ndime+2) = rgasc / xpres
!!$        
!!$        do idofn=1,ndime
!!$           xftra_inv(idofn,1) = xvelo(idofn) * xtunk(ndime+1) / velmo
!!$           xftra_inv(idofn,ndime+1) = xvelo(idofn) * xtunk(ndime+1) / sound
!!$           xftra_inv(idofn,ndime+2) = - xvelo(idofn) * xpres / xhecv / sound / sound
!!$        end do
!!$        xftra_inv(ndime+2,1) = velmo * xtunk(ndime+1)
!!$        xftra_inv(1,2) = - xvelo(2) * xtunk(ndime+1) / velmo_xy
!!$        xftra_inv(2,2) = xvelo(1) * xtunk(ndime+1) / velmo_xy
!!$        if (ndime==3) then
!!$           xftra_inv(1,3) = - xvelo(1) * xvelo(3) *  xtunk(ndime+1) / velmo / velmo_xy
!!$           xftra_inv(2,3) = - xvelo(2) * xvelo(3) *  xtunk(ndime+1) / velmo / velmo_xy
!!$           xftra_inv(3,3) = velmo_xy *  xtunk(ndime+1) / velmo
!!$        end if
!!$        xftra_inv(ndime+1,ndime+1) = xtunk(ndime+1) / sound
!!$        xftra_inv(ndime+2,ndime+1) = xtunk(ndime+1) * sound * (xhecv/rgasc+0.5_rp*xmach*xmach)
!!$        xftra_inv(ndime+1,ndime+2) = - xpres / xhecv / sound / sound 
!!$        xftra_inv(ndime+2,ndime+2) = - 0.5_rp * xpres * xmach * xmach / xhecv
!!$        
!!$        do idofn=1,ndofn_nsa
!!$           do jdofn=1,ndofn_nsa
!!$              do ldofn=1,ndofn_nsa
!!$                 do mdofn=1,ndofn_nsa
!!$                    xlopr_conservative(idofn,jdofn) = xlopr_conservative(idofn,jdofn) + xftra_inv(idofn,ldofn) &
!!$                         *xlopr(ldofn,mdofn)*xftra(mdofn,jdofn)
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
!!! END: VLR preconditioner with conservative variables (we apply a change of variables) 


!!! BEGIN: VLR preconditioner with conservative variables (we apply two change of variables)
        xttra(ndime+1,ndime+2) = rgacv / xtunk(ndime+1) / sound
        xttra(ndime+2,ndime+1) = - xhecv * sound * sound / xpres
        xttra(ndime+2,ndime+2) = rgasc / xpres

        xttra_inv(ndime+1,ndime+1) = xtunk(ndime+1) / sound
        xttra_inv(ndime+1,ndime+2) = - xpres / sound / sound / xhecv
        xttra_inv(ndime+2,ndime+1) = xtunk(ndime+1) * sound / rgacv

        xqtra(ndime+1,ndime+1) = 1.0_rp
        xqtra(ndime+2,ndime+2) = 1.0_rp
        xqtra(2,1) = - xvelo(2) / velmo_xy
        xqtra(2,2) = xvelo(1) / velmo_xy

        if (velmo_xy < zensa) then
           xqtra(2,1) = -1.0_rp !!(MM) Veure signe!!!
           xqtra(2,2) = 1.0_rp !!(MM) Veure signe!!!
        end if

        xqtra_inv(ndime+1,ndime+1) = 1.0_rp
        xqtra_inv(ndime+2,ndime+2) = 1.0_rp
        xqtra_inv(1,2) = xqtra(2,1)
        xqtra_inv(2,2) = xqtra(2,2)

        do idime=1,ndime
           xttra(idime,idime) = 1.0_rp / xtunk(ndime+1)
           xttra(idime,ndime+1) = - xvelo(idime) / xtunk(ndime+1)
           xttra(ndime+1,idime) = - rgacv * xvelo(idime) / xtunk(ndime+1) / sound
           xttra(ndime+1,ndime+1) = xttra(ndime+1,ndime+1) + 0.5_rp * rgacv * xvelo(idime) * xvelo(idime) / xtunk(ndime+1) / sound
           xttra(ndime+2,idime) = - rgasc * xvelo(idime) / xpres
           xttra(ndime+2,ndime+1) = xttra(ndime+2,ndime+1) + 0.5_rp * rgasc * xvelo(idime) * xvelo(idime) / xpres

           xttra_inv(idime,idime) = xtunk(ndime+1)
           xttra_inv(idime,ndime+1) = xtunk(idime) / sound
           xttra_inv(idime,ndime+2) = - xvelo(idime) * xpres / sound / sound / xhecv
           xttra_inv(ndime+2,idime) = xtunk(idime)
           xttra_inv(ndime+2,ndime+1) = xttra_inv(ndime+2,ndime+1) + 0.5_rp * xtunk(ndime+1) * xvelo(idime) * xvelo(idime) / sound
           xttra_inv(ndime+2,ndime+2) = xttra_inv(ndime+2,ndime+2) - 0.5_rp * xpres * xvelo(idime) * xvelo(idime) / xhecv / sound / sound

           xqtra(1,idime) = xvelo(idime) / velmo           
           if (velmo < zensa) then
              xqtra(1,idime) = 1.0_rp !!(MM) Veure signe!!!
           end if
           xqtra_inv(idime,1) = xqtra(1,idime)           
        end do

        if (ndime == 3) then
           xqtra(3,1) = - xvelo(1) * xvelo(3) / velmo / velmo_xy
           xqtra(3,2) = - xvelo(2) * xvelo(3) / velmo / velmo_xy
           xqtra(3,3) = velmo_xy / velmo

           if (velmo_xy < zensa) then
              xqtra(3,1) = -1.0_rp  !!(MM) Veure signe!!!
              xqtra(3,2) = -1.0_rp  !!(MM) Veure signe!!!
           end if
           if (velmo < zensa) then
              xqtra(3,3) = 1.0_rp
           end if

           xqtra_inv(1,3) = xqtra(3,1)
           xqtra_inv(2,3) = xqtra(3,2)
           xqtra_inv(3,3) = xqtra(3,3)
        end if

        do idofn=1,ndofn_nsa
           do jdofn=1,ndofn_nsa
              do kdofn=1,ndofn_nsa
                 do ldofn=1,ndofn_nsa
                    do mdofn=1,ndofn_nsa
                       do odofn=1,ndofn_nsa
                          xlopr_conservative(idofn,jdofn) &
                               = xlopr_conservative(idofn,jdofn) + xttra_inv(idofn,kdofn)*xqtra_inv(kdofn,ldofn) &
                               *xlopr(ldofn,mdofn)*xqtra(mdofn,odofn)*xttra(odofn,jdofn)
                       end do
                    end do
                 end do
              end do
           end do
        end do
!!! END: VLR preconditioner with conservative variables (we apply two change of variables) 

        ! Compute P_conservative * A
        do idofn=1,ndofn_nsa
           do jdofn=1,ndofn_nsa
              if (kfl_pseud_nsa == 1) &
                   xtide(idofn) = xtide(idofn) &
                   + xlopr_conservative(idofn,jdofn) * (xunkn(jdofn,igaus,ITER_K)-xunkn(jdofn,igaus,TIME_N))
              do idime=1,ndime
                 do pdofn=1,ndofn_nsa
                    xconv(idofn,jdofn,idime) = &
                         xconv(idofn,jdofn,idime) + xlopr_conservative(idofn,pdofn)*xconv_ori(pdofn,jdofn,idime)
                 end do
              end do
           end do
           if (kfl_pseud_nsa == 1) xtide(idofn) = xtide(idofn) * dtinv
        end do

!!$     do idofn=1,ndofn_nsa
!!$        do jdofn=1,ndofn_nsa
!!$           do idime=1,ndime
!!$              do kdofn=1,ndofn_nsa
!!$                 do ldofn=1,ndofn_nsa
!!$                    do mdofn=1,ndofn_nsa
!!$                       do odofn=1,ndofn_nsa
!!$                          do pdofn=1,ndofn_nsa
!!$                             xconv(idofn,jdofn,idime) = xconv(idofn,jdofn,idime) + xttra_inv(idofn,kdofn)*xqtra_inv(kdofn,ldofn) &
!!$                                  *xlopr(ldofn,mdofn)*xqtra(mdofn,odofn)*xttra(odofn,pdofn)*xconv_ori(pdofn,jdofn,idime)
!!$                          end do
!!$                       end do
!!$                    end do
!!$                 end do
!!$              end do
!!$           end do
!!$        end do
!!$     end do

     end if

  else if (kfl_pseud_nsa == 1) then

     do idofn=1,ndofn_nsa
        xtide(idofn) = xunkn(idofn,igaus,ITER_K)-xunkn(idofn,igaus,TIME_N)
        xtide(idofn) = xtide(idofn) * dtinv
     end do

  end if


  !
  ! volume forces: gravity, dissipation (if LES) and diffusion of enthalpy (if species)
  !
  xvofo(ndime +2,ndime + 1) = xvofo(ndime +2,ndime + 1) + sgsdi + dhtra

  xresi(ndime+2) = xresi(ndime+2) + heats  ! Source term of chemical reactions

  do idofn=1,ndofn_nsa
     if (kfl_pseud_nsa == 1) xresi(idofn) = xresi(idofn) - xtide(idofn)
     if (kfl_skews_nsa == 1) xresi(idofn) = xresi(idofn) - dflux_conv(idofn)

     do jdofn=1,ndofn_nsa
        xvofo(idofn,jdofn) = xvofo(idofn,jdofn) + gravm_nsa(idofn,jdofn)
        xresi(idofn)= xresi(idofn) + xvofo(idofn,jdofn)*xtunk(jdofn)
        do idime=1,ndime
           xresi(idofn)= xresi(idofn) - xconv(idofn,jdofn,idime) * gunkn(jdofn,idime)
           do jdime=1,ndime
              if (kfl_lopre_nsa < 2) then
                 xresi(idofn) = xresi(idofn) + ddiff(idofn,jdofn,1,idime) * gunkn(jdofn,jdime) + &
                      xdiff(idofn,jdofn,idime,jdime) * hunkn(jdofn,idime,jdime)
              else
                 ! (MM) This is an approximation. ddiff(idofn,jdofn,1,idime) * gunkn(jdofn,jdime) tendría que estar
                 ! Pero ddiff no tiene el precondicionador dentro !!!! hay que programarlo !!!!!
                 xresi(idofn) = xresi(idofn) + xdiff(idofn,jdofn,idime,jdime) * hunkn(jdofn,idime,jdime)
              end if
              if  (kfl_taudi_nsa >= 5) then   ! non-diagonal tau
                 conme(idofn,jdofn,idime) = conme(idofn,jdofn,idime) + &
                      xconv(idofn,jdofn,jdime) * xjaci(idime,jdime) / real(pgaus,rp)
                 !!he d'incloure a difme el canvi de variables del pas a l'espai parametric
                 !              difme(idofn,jdofn,idime,jdime) = difme(idofn,jdofn,idime,jdime) + &
                 !                   xdiff(idofn,jdofn,idime,jdime) / real(pgaus,rp)
                 do kdime=1,ndime
                    difme(idofn,jdofn,idime,jdime) = difme(idofn,jdofn,idime,jdime) + &
                         xdiff(idofn,jdofn,jdime,kdime) * xjaci(idime,kdime) / real(pgaus,rp)
                 end do
              end if
           end do
        end do

     end do
  end do

  !    if (ielem == 1 .and. igaus == 1) then
  !       write(6,*) xtide(1:ndofn_nsa)
  !       write(6,*) xunkn(1:ndofn_nsa,igaus,ITER_K)
  !       write(6,*) xunkn(1:ndofn_nsa,igaus,TIME_N)
  !   end if


end subroutine nsa_gauvalxy
