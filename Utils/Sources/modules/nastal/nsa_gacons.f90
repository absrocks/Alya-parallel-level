subroutine nsa_gacons(&
     ndime,mnode,pnode,npoin,ngaus,mgaus,ncomp,ntens,&
     ivert,lnods, &
     plapl, &
     elvel,elpre,elpol,&
     elden,&
     elent,&
     eltem,&
     elvis,&
     elvit,&
     elene,&
     elumo,&
     elfmo,&
     elcod,&
     elken,&
     eleps,&
     elsub,&
     elqva,&
     elqva_hy,&
     elqcl,&
     elqra,&
     xvelo,&
     xpres,&
     xphyd,&
     xdens,&
     xenth,&
     xdhyd,&
     xtemp,&
     xthyd,&
     xconc,&
     xvisc,&
     xvitu,&
     xener,&
     xumom,&
     xfmom,&
     xsube,&
     dgvel,&
     xqvap,xqvap_hy,&
     xqclo,&
     xqrai,&
     gvelo,&
     gpres,gpold,gphyd,&
     gdens,gdhyd,&
     gtemp,gthyd,&
     gvisc,gvitu,&
     gqvap,&
     gqclo,&
     gqrai,&
     gener,&
     gumom,gfmom,&
     gdvel,&
     xkene,xepsi,&
     gkene,gepsi,&
     gkevo,&
     xlapr,&
     xlake,&
     xlate,&
     xlade,&
     xldve,xcoor,zeroc, &
     dvolu,cartd,hessi,heslo,gpsha,deriv,weigp,frano,gpcod,gpcor,wvect,taray, &
     elcon,eldif,xconv,xdiff,dconv,ddiff)
!-----------------------------------------------------------------------
!
! Gather operations and simultaneous gauss point evaluation
!
! OJOOOOOOOO:
! - TODAS LAS COSAS EN PUNTOS DE GAUSS VIENEN INTERPOLADAS, NADA SE CALCULA LOCALMENTE...
!   NO ES COMO EN ALAMAK, HABRIA QUE VER QUE PASA CON VISCOSIDAD (DE LEY), VISCOSIDAD TUMBULENTA, ETC...
!
!-----------------------------------------------------------------------------------------------------------
  use      def_kintyp
  use      def_master, only : veloc,&
                              press,&
                              entha,&
                              densi,&
                              tempe,&
                              visco,&
                              turmu,&
                              energ,&
                              umome,&
                              untur,&
                              conce

  use      def_domain, only : coord,kfl_naxis

  use      def_nastal, only : stapa_nsa,&
                              brunt_nsa,&
                              brure_nsa,&
                              kfl_brunt_nsa,&
                              rekee_nsa,&
                              kfl_foreg_nsa,&
                              kfl_turbu_nsa,&
                              kfl_relat_nsa,&
                              kfl_hysta_nsa,&
                              kfl_inkee_nsa,&
                              kfl_inifi_nsa,&
                              kfl_rayle_nsa,&
                              kfl_coupl_nsa,&
                              cvcoe_nsa,&
                              cpcoe_nsa,&
                              rgasc_nsa,&
                              grnor_nsa,&
                              pbaro_nsa,&
                              tempe_nsa,&
                              adgam_nsa,& 
                              ndofn_nsa,&
                              rgacv_nsa,&
                              afact_nsa,&
                              kfl_infun_nsa,&
                              kfl_benme_nsa,&
                              kfl_physics_nsa,&
                              kdiff_nsa
  
  implicit none

  integer(ip),  intent(in)  :: ndime,mnode,pnode,npoin,ngaus,mgaus,ncomp,ntens,plapl,ivert
  integer(ip),  intent(in)  :: lnods(pnode)
  integer(ip)               :: inode,ipoin,ipaux,idome,igaus,itidi,idime,jdime,kdime,itime,nprev,idofn,jdofn

  real(rp), intent(in)  :: &
       gpsha(pnode,ngaus),deriv(ndime,pnode,ngaus),&
       weigp(ngaus),heslo(ntens,pnode,ngaus),frano(4,2),zeroc

  real(rp), intent(out) :: &
       elvel(ndime,mnode,ncomp),&
       elene(mnode,ncomp),&
       elpre(mnode,ncomp),&
       elpol(mnode,ncomp),&
       elden(mnode,ncomp),&
       elent(mnode,ncomp),&
       eltem(mnode,ncomp),&
       elqva(mnode,ncomp),&
       elqva_hy(mnode,ncomp),&
       elqcl(mnode,ncomp),&
       elqra(mnode,ncomp),&
       elvis(mnode,ncomp),&
       elcod(ndime,mnode),&
       elumo(ndime,mnode,ncomp),&
       elfmo(ndime,mnode,ncomp),&
       elken(mnode,ncomp),&
       eleps(mnode,ncomp),&
       elvit(mnode,ncomp),&
       xdens(mgaus,ncomp),&
       xenth(mgaus,ncomp),&
       xdhyd(mgaus),&
       xpres(mgaus,ncomp),&
       xphyd(mgaus),&
       xtemp(mgaus,ncomp),&
       xthyd(mgaus),&
       xconc(mgaus,ncomp),&
       xkene(mgaus,ncomp),&
       xqvap(mgaus,ncomp),&
       xqvap_hy(mgaus,ncomp),&
       xqclo(mgaus,ncomp),&
       xqrai(mgaus,ncomp),&
       xvisc(mgaus,ncomp),&
       xvitu(mgaus,ncomp),&
       xener(mgaus,ncomp),&
       xepsi(mgaus,ncomp),&
       xvelo(ndime,mgaus,ncomp),&
       xfmom(ndime,mgaus,ncomp),&
       xumom(ndime,mgaus,ncomp),&
       xlapr(mgaus,ncomp),&
       xlate(mgaus,ncomp),&
       xlake(mgaus,ncomp),&
       xlade(mgaus,ncomp),&
       xldve(ndime,mgaus,ncomp),&
       gdens(ndime,mgaus,ncomp),&
       gpres(ndime,mgaus,ncomp),&
       gpold(ndime,mgaus,ncomp),&
       gphyd(ndime,mgaus),&
       gdhyd(ndime,mgaus),&
       gthyd(ndime,mgaus),&
       gtemp(ndime,mgaus,ncomp),&
       gqvap(ndime,mgaus,ncomp),&
       gqclo(ndime,mgaus,ncomp),&
       gqrai(ndime,mgaus,ncomp),&
       gener(ndime,mgaus,ncomp),&
       gvisc(ndime,mgaus,ncomp),&
       gvitu(ndime,mgaus,ncomp),&
       gkene(ndime,mgaus,ncomp),&
       gepsi(ndime,mgaus,ncomp),&
       xcoor(ndime,mgaus),&
       taray(ndofn_nsa,ndofn_nsa,mgaus)
  
  real(rp), intent(out) :: &
       gvelo(ndime,ndime,mgaus,ncomp),gfmom(ndime,ndime,mgaus,ncomp),&
       gumom(ndime,ndime,mgaus,ncomp),gkevo(ndime,mgaus,ncomp),dvolu(mgaus),&
       cartd(ndime,mnode,mgaus),hessi(ntens,mnode,mgaus),gdvel(ndime,mgaus,ncomp),&
       dgvel(ndime,mgaus,ncomp),wvect(ndime,mgaus),gpcod(ndime,mgaus),gpcor(mgaus),&
       elcon(ndofn_nsa,ndofn_nsa,ndime,mnode),eldif(ndofn_nsa,ndofn_nsa,ndime,ndime,mnode), &
       xconv(ndofn_nsa,ndofn_nsa,ndime,mgaus),xdiff(ndofn_nsa,ndofn_nsa,ndime,ndime,mgaus), &
       dconv(ndofn_nsa,ndofn_nsa,mgaus),ddiff(ndofn_nsa,ndofn_nsa,2,ndime,mgaus)

  real(rp)              :: detjm,xshap,xcart,xhess,xdelo,xtelo,xprlo,xcolo
  real(rp)              :: xmvlo !SM: local mixing ratio of vapor
  real(rp)              :: d2sdx(ndime,ndime,ndime)
  real(rp)              :: auxi1(ndime,ndime,mnode)
  real(rp)              :: auxi2(ndime,ndime,mnode)
  real(rp)              :: xjaci(ndime,ndime),xjacm(ndime,ndime) 
  real(rp)              :: eldhy(mnode),elphy(mnode),elthy(mnode),elsub(ndime+2,mnode,2),xsube(ndime+2,mgaus,2)
  real(rp)              :: visci, dicod, velno(ndime), velsq
  real(rp)              :: xtunk(ndime+2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!miquel
!  integer*8 :: cpu_time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!miquel

  !
  ! Gathering of the current values 
  !

!     tstart_gath = cpu_time()
  nprev= 4

  !
  ! Hydrostatic density, temperature and pressure on nodes
  !
  elfmo= 0.0_rp
  do inode= 1,pnode
     eldhy(inode) = 0.0_rp
     elthy(inode) = 0.0_rp 
     elphy(inode) = 0.0_rp
     if (kfl_infun_nsa == 1) then
        ipoin= lnods(inode)
        eldhy(inode) = rekee_nsa(ndime+1,ipoin)
        elthy(inode) = rekee_nsa(ndime+2,ipoin)
        elphy(inode) = rekee_nsa(ndime+3,ipoin)
     end if
  end do
  inode_loop: do inode= 1,pnode
     ipoin= lnods(inode)
     ! Density on nodes
     elden(inode,1) = densi(ipoin,1)
     elden(inode,2) = densi(ipoin,2)
     elden(inode,3) = densi(ipoin,nprev)

     ! Temperature on nodes
     eltem(inode,1) = tempe(ipoin,1)
     eltem(inode,2) = tempe(ipoin,2)
     eltem(inode,3) = tempe(ipoin,nprev)
     ! Energy on nodes
     elene(inode,1) = energ(ipoin,1)        
     elene(inode,2) = energ(ipoin,2)        
     elene(inode,3) = energ(ipoin,nprev)        
     ! Viscosity on nodes
     elvis(inode,1) = visco(ipoin,1)
     !
     ! Kessler variables
     !
     if(kfl_benme_nsa >= 200 .and. kfl_physics_nsa > 0) then

        elqva(inode,1) = conce(ipoin,1,1)
        !elqva(inode,2) = conce(ipoin,1,2)
        !elqva(inode,3) = conce(ipoin,1,nprev)

        elqcl(inode,1) = conce(ipoin,2,1)
        !elqcl(inode,2) = conce(ipoin,2,2)
        !elqcl(inode,3) = conce(ipoin,2,nprev)

        elqra(inode,1) = conce(ipoin,3,1)
        !elqra(inode,2) = conce(ipoin,3,2)
        !elqra(inode,3) = conce(ipoin,3,nprev)

        elqva_hy(inode,1) = rekee_nsa(ndime+3 + 1, ipoin)
        !elqva_hy(inode,2) = rekee_nsa(ndime+3 + 2, ipoin)
        !elqva_hy(inode,3) = rekee_nsa(ndime+3 + 3, ipoin)
     end if

     ! Coordinates, Velocity and Momentum on nodes
     do idime=1,ndime
        elcod(idime,inode  ) = coord(idime,ipoin  )          
        elvel(idime,inode,1) = veloc(idime,ipoin,1)
        elvel(idime,inode,2) = veloc(idime,ipoin,2)
        elvel(idime,inode,3) = veloc(idime,ipoin,nprev)
        elumo(idime,inode,1) = umome(idime,ipoin,1)
        elumo(idime,inode,2) = umome(idime,ipoin,2)
        elumo(idime,inode,3) = umome(idime,ipoin,nprev)
        do idofn= 1,ndofn_nsa
           do jdofn= 1,ndofn_nsa
              elcon(idofn,jdofn,idime,inode) = 0.0_rp
              do jdime= 1,ndime
                 eldif(idofn,jdofn,idime,jdime,inode) = 0.0_rp
              end do
           end do
        end do        
     end do

!print*,elvel(idime,inode,1),elvel(idime,inode,2),elvel(idime,inode,3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! OJO: El eldif(ndime+2,.,.,.) esta puesto a cero pero no tiene que ser asi!!!!!!!!!
     visci = elvis(inode,1) / (elden(inode,1) + elsub(ndime+1,inode,1))
     do kdime=1,ndime        
        elcon(ndime+1,kdime  ,kdime,inode)= 1.0_rp
        elcon(kdime  ,ndime+2,kdime,inode)= afact_nsa*adgam_nsa*(elden(inode,1)*eltem(inode,1))**(adgam_nsa-1.0_rp)*elden(inode,1)
        elcon(kdime  ,ndime+1,kdime,inode)= afact_nsa*adgam_nsa*(elden(inode,1)*eltem(inode,1))**(adgam_nsa-1.0_rp)*eltem(inode,1)
        elcon(ndime+2,ndime+1,kdime,inode)= 0.0_rp
        elcon(ndime+2,ndime+2,kdime,inode)= elvel(kdime,inode,1)
        eldif(ndime+2,ndime+1,kdime,kdime,inode)= 0.0_rp
        eldif(ndime+2,ndime+2,kdime,kdime,inode)= 0.0_rp
        do idime=1,ndime
           elcon(idime,idime  ,kdime,inode)= elcon(idime,idime  ,kdime,inode) + velno(kdime) 
           elcon(idime,kdime  ,kdime,inode)= elcon(idime,kdime  ,kdime,inode) + velno(idime)
           elcon(idime,ndime+1,kdime,inode)= elcon(idime,ndime+1,kdime,inode) - velno(idime) * velno(kdime)
           elcon(ndime+2,idime,kdime,inode)= 0.0_rp
           eldif(kdime,kdime,idime,idime,inode)= visci
           eldif(kdime,idime,idime,kdime,inode)= eldif(kdime,idime,idime,kdime,inode) + visci
           eldif(kdime,idime,kdime,idime,inode)= eldif(kdime,idime,kdime,idime,inode) - 2.0_rp * visci / 3.0_rp
           eldif(kdime,ndime+1,idime,idime,inode)= - visci * velno(kdime)
           eldif(kdime,ndime+1,idime,kdime,inode)= eldif(kdime,ndime+1,idime,kdime,inode) - visci * velno(idime)
           eldif(kdime,ndime+1,kdime,idime,inode)= eldif(kdime,ndime+1,kdime,idime,inode) &
                + 2.0_rp * visci * velno(idime) / 3.0_rp
           eldif(ndime+2,idime,kdime,kdime,inode)= 0.0_rp
           eldif(ndime+2,kdime,kdime,idime,inode)= 0.0_rp
           eldif(ndime+2,kdime,idime,kdime,inode)= 0.0_rp
           eldif(ndime+2,ndime+1,kdime,idime,inode)= 0.0_rp
        end do
     end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


     elpol(inode,1) = press(ipoin,nprev)      ! old p, for fract. step  
     ! Pressure on nodes
     elpre(inode,1) = press(ipoin,1)          ! new p, for fract. step

     ! Turbulence terms
     if (kfl_turbu_nsa /= 0) then
        elvit(inode,1) = turmu(  ipoin  )
        elken(inode,1) = untur(1,ipoin,1)
        eleps(inode,1) = untur(2,ipoin,1)
     else
        elvit(inode,1) = 0.0_rp
        elken(inode,1) = 0.0_rp
        eleps(inode,1) = 0.0_rp
     end if
  end do inode_loop

!     tend_gath = cpu_time()-tstart_gath
!     ttotal_gath = ttotal_gath + tend_gath*69.84E-9



!!!!!!!!!!!!!!!!!!!---mirar---CREC QUE ES POT TREURE BASTANTA COSA
  xdelo= 0.0_rp
  xprlo= 0.0_rp
  xtelo= 0.0_rp
  xmvlo= 0.0_rp
  if (kfl_hysta_nsa == 1) then
     do inode=1,pnode
        ipoin=lnods(inode)
        xcolo= elcod(ivert,inode)
        if (kfl_brunt_nsa == 2) then
           stapa_nsa(1) = brunt_nsa(ipoin)
        else
           stapa_nsa(1) = brure_nsa
        end if

!        if ((kfl_inkee_nsa(2) + kfl_inkee_nsa(3) + kfl_inkee_nsa(4)) == 0) then
!           call nsa_stalaw(20,xdelo,xprlo,xtelo,xcolo,xmvlo)
!           call nsa_stalaw(30,xdelo,xprlo,xtelo,xcolo,xmvlo)
!           call nsa_stalaw( 1,xdelo,xprlo,xtelo,xcolo,xmvlo)       
!        else if (kfl_inkee_nsa(2) > 0 .and. kfl_inkee_nsa(4) > 0) then
!           if (kfl_inifi_nsa(2) == 2) then
!              xprlo= rekee_nsa(ndime+3    ,ipoin)
!              xdelo= rekee_nsa(ndime+1    ,ipoin)              
!           else
!              call runend('NSA_GACONS: PRESSURE INITIAL CONDITION MUST BE GIVEN!')
!           end if
!        end if
!        eldhy(inode)= xdelo
!        elphy(inode)= xprlo
     end do     
  end if
!!!!!!!!!!!!!!!!!!!!!---mirar---





  ! Rayleigh dumping in 
  if (kfl_rayle_nsa == 1) then
     call nsa_raydam(lnods,pnode,ngaus,ivert,gpsha,elcod,elumo(1,1,1),elden(1,1),taray)
  end if

  if (kfl_relat_nsa == 1) then
     do inode=1,pnode
        ipoin=lnods(inode)
        elent(inode,1) = entha(ipoin,1)
     end do
  end if

!     tstart_gacons_gaus_loop1  = cpu_time()


  gauss_points_loop: do igaus=1,ngaus
!     tstart_gacons_elmder2_loop = cpu_time()
     call elmder(pnode,ndime,deriv(1,1,igaus),elcod,cartd(1,1,igaus),detjm,xjacm,xjaci)
     dvolu(igaus)=weigp(igaus)*detjm        
!     tend_gacons_elmder2_loop = cpu_time()-tstart_gacons_elmder2_loop
!     ttotal_gacons_elmder2_loop = ttotal_gacons_elmder2_loop + tend_gacons_elmder2_loop*69.84E-9     
     hessi(1:ntens,1:mnode,igaus) = 0.0_rp
     if(plapl==1 .or. kdiff_nsa > 0.0_rp) then
!          tstart_gacons_elmhes2_loop = cpu_time()
          call elmhes(heslo(1,1,igaus),hessi(1,1,igaus),ndime,pnode,ntens,&
          xjaci,d2sdx,deriv(1,1,igaus),elcod)
!          tend_gacons_elmhes2_loop = cpu_time()-tstart_gacons_elmhes2_loop
!          ttotal_gacons_elmhes2_loop = ttotal_gacons_elmhes2_loop + tend_gacons_elmhes2_loop*69.84E-9
     end if
     
     wvect(1:ndime,igaus)=0.0_rp
     if(frano(1,2)>zeroc.or.frano(2,2)>zeroc.or.kfl_naxis==1) then 
        call mbvab0(gpcod(1,igaus),elcod,gpsha(1,igaus),ndime,pnode)
        wvect(1:ndime,igaus)=gpcod(1:ndime,igaus)
        gpcor(igaus)=gpcod(1,igaus)                    
     end if
!    tend_gacons_gaus_loop1 = cpu_time()-tstart_gacons_gaus_loop1
!    ttotal_gacons_gaus_loop1 = ttotal_gacons_gaus_loop1 + tend_gacons_gaus_loop1*69.84E-9
!    tstart_gacons_gaus_loop2  = cpu_time()

     ! Initialize gauss point values
     xdens(igaus,1)=0.0_rp 
     xdens(igaus,2)=0.0_rp 
     xdens(igaus,3)=0.0_rp 
     xdhyd(igaus)  =0.0_rp 

     xpres(igaus,1)=0.0_rp 
     xpres(igaus,2)=0.0_rp 
     xpres(igaus,3)=0.0_rp
     xphyd(igaus)  =0.0_rp
 
     xtemp(igaus,1)=0.0_rp 
     xtemp(igaus,2)=0.0_rp 
     xtemp(igaus,3)=0.0_rp 
     xthyd(igaus)  =0.0_rp

     xenth(igaus,1)=0.0_rp 
     xenth(igaus,2)=0.0_rp 
     xenth(igaus,3)=0.0_rp 

     
     !
     ! Kessler variables
     !
     if(kfl_benme_nsa >= 200 .and. kfl_physics_nsa > 0) then
        xqvap(igaus,1)=0.0_rp 
        !xqvap(igaus,2)=0.0_rp 
        !xqvap(igaus,3)=0.0_rp 

        xqclo(igaus,1)=0.0_rp 
        !xqclo(igaus,2)=0.0_rp 
        !xqclo(igaus,3)=0.0_rp 

        xqrai(igaus,1)=0.0_rp 
        !xqrai(igaus,2)=0.0_rp 
        !xqrai(igaus,3)=0.0_rp 

        xqvap_hy(igaus,1)=0.0_rp 
        !xqvap_hy(igaus,2)=0.0_rp 
        !xqvap_hy(igaus,3)=0.0_rp 
     end if

     xconc(igaus,1)=0.0_rp 
     xconc(igaus,2)=0.0_rp 
     xconc(igaus,3)=0.0_rp 
     xener(igaus,1)=0.0_rp 
     xener(igaus,2)=0.0_rp 
     xener(igaus,3)=0.0_rp 
     xvisc(igaus,1)=0.0_rp 
     xvitu(igaus,1)=0.0_rp 
     xkene(igaus,1)=0.0_rp 
     xepsi(igaus,1)=0.0_rp 
     xlapr(igaus,1)=0.0_rp 
     xlake(igaus,1)=0.0_rp 
     xlade(igaus,1)=0.0_rp 
     xlate(igaus,1)=0.0_rp 
     xcoor(1,igaus)=0.0_rp 
     xcoor(2,igaus)=0.0_rp 
     xcoor(ndime,igaus)=0.0_rp      
     do idime=1,ndime
        xvelo(idime,igaus,1)=0.0_rp 
        xvelo(idime,igaus,2)=0.0_rp 
        xvelo(idime,igaus,3)=0.0_rp 
        xldve(idime,igaus,1)=0.0_rp 
        xumom(idime,igaus,1)=0.0_rp 
        xumom(idime,igaus,2)=0.0_rp 
        xumom(idime,igaus,3)=0.0_rp 
        xfmom(idime,igaus,1)=0.0_rp 
        gdens(idime,igaus,1)=0.0_rp 
        gdens(idime,igaus,2)=0.0_rp
        gdhyd(idime,igaus)  =0.0_rp
        gpres(idime,igaus,1)=0.0_rp 
        gpres(idime,igaus,2)=0.0_rp 
        gpold(idime,igaus,1)=0.0_rp 
        gphyd(idime,igaus)  =0.0_rp 
        gtemp(idime,igaus,1)=0.0_rp 
        gtemp(idime,igaus,2)=0.0_rp
        gthyd(idime,igaus)  =0.0_rp
        gener(idime,igaus,1)=0.0_rp 
        gvisc(idime,igaus,1)=0.0_rp 
        gvitu(idime,igaus,1)=0.0_rp 
        gepsi(idime,igaus,1)=0.0_rp 
        gkene(idime,igaus,1)=0.0_rp 
        gkevo(idime,igaus,1)=0.0_rp 
        dgvel(idime,igaus,1)=0.0_rp 
        gdvel(idime,igaus,1)=0.0_rp 
        do jdime=1,ndime
           gvelo(idime,jdime,igaus,1)=0.0_rp 
           gfmom(idime,jdime,igaus,1)=0.0_rp 
           gumom(idime,jdime,igaus,1)=0.0_rp 
        end do
        do inode=1,pnode
           xshap = gpsha(inode,igaus)
           xcoor(idime,igaus) = xcoor(idime,igaus) + elcod(idime,inode) * xshap
        end do
     end do
     do idofn= 1,ndofn_nsa
        xsube(idofn,igaus,1) = 0.0_rp
        xsube(idofn,igaus,2) = 0.0_rp
        do jdofn=1,ndofn_nsa
           do idime=1,ndime
              xconv(idofn,jdofn,idime,igaus) = 0.0_rp           
!!$              ddiff(idofn,jdofn,1,idime,igaus) = 0.0_rp           
!!$              ddiff(idofn,jdofn,2,idime,igaus) = 0.0_rp           
              do jdime=1,ndime
!!$                 xdiff(idofn,jdofn,idime,jdime,igaus) = 0.0_rp           
             end do
           end do
           dconv(idofn,jdofn,igaus) = 0.0_rp
       end do
     end do

     ! Compute Hydrostatic Temperature at gauss point (xthyd)
     call nsa_compute_hyd(3,xdhyd(igaus),xphyd(igaus),xthyd(igaus),xcoor(ndime,igaus))

     ! Compute Hydrostatic Density at gauss point (xdhyd)
     call nsa_compute_hyd(1,xdhyd(igaus),xphyd(igaus),xthyd(igaus),xcoor(ndime,igaus))
     ! Compute Hydrostatic Pressure at gauss point (xphyd)
!     call nsa_compute_hyd(2,xdhyd(igaus),xphyd(igaus),xthyd(igaus),xcoor(ndime,igaus))
          
     ! Compute gauss values for press and veloc, which can be done now
     do inode=1,pnode
        xshap = gpsha(inode,igaus)
        do idime=1,ndime
           xcart=cartd(idime,inode,igaus) 
           xfmom(idime,igaus,1) = xfmom(idime,igaus,1) + &
                elfmo(idime,inode,1) * xshap
           xumom(idime,igaus,1) = xumom(idime,igaus,1) + &
                elumo(idime,inode,1) * xshap
           xumom(idime,igaus,2) = xumom(idime,igaus,2) + &
                elumo(idime,inode,2) * xshap
           xumom(idime,igaus,3) = xumom(idime,igaus,3) + &
                elumo(idime,inode,3) * xshap
           do jdime=1,ndime
              gfmom(idime,jdime,igaus,1) = gfmom(idime,jdime,igaus,1)  + &
                   elfmo(jdime,inode,1) * xcart
              gumom(idime,jdime,igaus,1) = gumom(idime,jdime,igaus,1)  + &
                   elumo(jdime,inode,1) * xcart
           end do
        end do
        if (plapl==1 .or. kdiff_nsa > 0.0_rp) then

           if(ndime.eq.2) then
              gdvel(1,igaus,1)=gdvel(1,igaus,1)+hessi(1,inode,igaus)*elvel(1,inode,1) &
                   +hessi(3,inode,igaus)*elvel(2,inode,1)
              gdvel(2,igaus,1)=gdvel(2,igaus,1)+hessi(3,inode,igaus)*elvel(1,inode,1) &
                   +hessi(2,inode,igaus)*elvel(2,inode,1)
           else if(ndime.eq.3) then
              gdvel(1,igaus,1)=gdvel(1,igaus,1)+hessi(1,inode,igaus)*elvel(1,inode,1) & 
                   +hessi(4,inode,igaus)*elvel(2,inode,1) &
                   +hessi(5,inode,igaus)*elvel(3,inode,1)
              gdvel(2,igaus,1)=gdvel(2,igaus,1)+hessi(4,inode,igaus)*elvel(1,inode,1) & 
                   +hessi(2,inode,igaus)*elvel(2,inode,1) &
                   +hessi(6,inode,igaus)*elvel(3,inode,1)
              gdvel(3,igaus,1)=gdvel(3,igaus,1)+hessi(5,inode,igaus)*elvel(1,inode,1) & 
                   +hessi(6,inode,igaus)*elvel(2,inode,1) &
                   +hessi(3,inode,igaus)*elvel(3,inode,1)
           end if
           do idime=1,ndime
              xhess = hessi(idime,inode,igaus)
              xlapr(igaus,1) = xlapr(igaus,1) + elpre(inode,1)*xhess
              xlade(igaus,1) = xlade(igaus,1) + elden(inode,1)*xhess
              xlate(igaus,1) = xlate(igaus,1) + eltem(inode,1)*xhess
              xlake(igaus,1) = xlake(igaus,1) + elken(inode,1)*xhess
              xldve(idime,igaus,1) = xldve(idime,igaus,1) + elden(inode,1)*xhess
              do jdime=1,ndime
                 dgvel(idime,igaus,1)=dgvel(idime,igaus,1) + &
                      hessi(idime,inode,igaus)*elvel(jdime,inode,1)
              end do
           end do
        end if
     end do
  end do gauss_points_loop

!     tend_gacons_gaus_loop1 = cpu_time()-tstart_gacons_gaus_loop1
!     ttotal_gacons_gaus_loop1 = ttotal_gacons_gaus_loop1 + tend_gacons_gaus_loop1*69.84E-9

  
  if (kfl_foreg_nsa == 0) then     

!     tstart_gacons_gaus_loop3  = cpu_time()

!------- BUCLE ORIGINAL -------------miquel
!!$     do inode=1,pnode
!!$        ! Gauss point values
!!$        do igaus=1,ngaus
!------- BUCLE OPTIMITZAT -------------miquel
     
     
     !
     ! Gauss point values
     !    
     gauss_points_loop_2:     do igaus=1,ngaus
        xdens(igaus,1) = xdhyd(igaus) 
        xtemp(igaus,1) = xthyd(igaus)
        ! Derivative of hydrostatic density (gdhyd)
        call nsa_compute_hyd_der(1,gdhyd(ndime,igaus),gphyd(ndime,igaus), &
             gthyd(ndime,igaus),xcoor(ndime,igaus))
        
        ! Derivative of hydrostatic temperature (gthyd)
        call nsa_compute_hyd_der(3,gdhyd(ndime,igaus),gphyd(ndime,igaus), &
             gthyd(ndime,igaus),xcoor(ndime,igaus))

        ! Derivative of hydrostatic pressure (gphyd)
        ! IT MUST BE ADDED TO nsa_compute_hyd_der, but with the addition
        ! of the missing variables to be passed (i.e. xthyd, xdhyd)
        gphyd(ndime,igaus) = pbaro_nsa*adgam_nsa*(rgasc_nsa/pbaro_nsa)**adgam_nsa * &
             (xthyd(igaus)*xdhyd(igaus))**(adgam_nsa-1.0_rp) * &
             (xthyd(igaus)*gdhyd(ndime,igaus) + xdhyd(igaus)*gthyd(ndime,igaus))

        gdens(ndime,igaus,1) = gdhyd(ndime,igaus)
        gtemp(ndime,igaus,1) = gthyd(ndime,igaus)
        do inode=1,pnode
           !Shape function:
           xshap = gpsha(inode,igaus)
           !DENS:
           xdens(igaus,1) = xdens(igaus,1) + xshap*(elden(inode,1)-eldhy(inode))
           xdens(igaus,2) = xdens(igaus,2) + xshap*elden(inode,2) 
           xdens(igaus,3) = xdens(igaus,3) + xshap*elden(inode,3) 
           !TEMPE:
           xtemp(igaus,1) = xtemp(igaus,1) + xshap*(eltem(inode,1)-elthy(inode)) 
           xtemp(igaus,2) = xtemp(igaus,2) + xshap*eltem(inode,2) 
           xtemp(igaus,3) = xtemp(igaus,3) + xshap*eltem(inode,3) 
           
           !
           ! Kessler variables
           !
           if(kfl_benme_nsa >= 200 .and. kfl_physics_nsa > 0) then
              xqvap(igaus,1) = xqvap(igaus,1) + xshap*elqva(inode,1) 
              !xqvap(igaus,2) = xqvap(igaus,2) + xshap*elqva(inode,2) 
              !xqvap(igaus,3) = xqvap(igaus,3) + xshap*elqva(inode,3)

              xqclo(igaus,1) = xqclo(igaus,1) + xshap*elqcl(inode,1) 
              !xqclo(igaus,2) = xqclo(igaus,2) + xshap*elqcl(inode,2) 
              !xqclo(igaus,3) = xqclo(igaus,3) + xshap*elqcl(inode,3)

              xqrai(igaus,1) = xqrai(igaus,1) + xshap*elqra(inode,1) 
              !xqrai(igaus,2) = xqrai(igaus,2) + xshap*elqra(inode,2) 
              !xqrai(igaus,3) = xqrai(igaus,3) + xshap*elqra(inode,3)

              xqvap_hy(igaus,1) = xqvap_hy(igaus,1) + xshap*elqva_hy(inode,1) 
              !xqvap_hy(igaus,2) = xqvap_hy(igaus,2) + xshap*elqva_hy(inode,2) 
              !xqvap_hy(igaus,3) = xqvap_hy(igaus,3) + xshap*elqva_hy(inode,3)
           end if
        
           ! xdens(igaus,1) = xdens(igaus,1) + xshap * elden(inode,1) 
          ! xdens(igaus,2) = xdens(igaus,2) + xshap * elden(inode,2) 
          ! xdens(igaus,3) = xdens(igaus,3) + xshap * elden(inode,3) 
          ! xdhyd(igaus  ) = xdhyd(igaus  ) + xshap * eldhy(inode  ) 
          ! xphyd(igaus  ) = xphyd(igaus  ) + xshap * elphy(inode  ) 
          ! xtemp(igaus,1) = xtemp(igaus,1) + xshap * eltem(inode,1) 
          ! xtemp(igaus,2) = xtemp(igaus,2) + xshap * eltem(inode,2) 
          ! xtemp(igaus,3) = xtemp(igaus,3) + xshap * eltem(inode,3) 
           !ENE:
           xener(igaus,1) = xener(igaus,1) + xshap * elene(inode,1) 
           xener(igaus,2) = xener(igaus,2) + xshap * elene(inode,2) 
           xener(igaus,3) = xener(igaus,3) + xshap * elene(inode,3) 
           !VISC Terms:
           xvisc(igaus,1) = xvisc(igaus,1) + xshap * elvis(inode,1) 
           xvitu(igaus,1) = xvitu(igaus,1) + xshap * elvit(inode,1) 
           
           do idime=1,ndime
              xcart=cartd(idime,inode,igaus)


              do idofn= 1,ndofn_nsa
                 do jdofn=1,ndofn_nsa
                    dconv(idofn,jdofn,igaus) = dconv(idofn,jdofn,igaus) + &
                         xcart * elcon(idofn,jdofn,idime,inode)
!!$                    do jdime=1,ndime
!!$                       ddiff(idofn,jdofn,1,idime) = ddiff(idofn,jdofn,1,idime) &
!!$                            + cartd(jdime,inode) * eldif(idofn,jdofn,jdime,idime,inode)
!!$                       ddiff(idofn,jdofn,2,idime) = ddiff(idofn,jdofn,2,idime) &
!!$                            + cartd(jdime,inode) * eldif(idofn,jdofn,idime,jdime,inode)
!!$                    end do
                 end do
              end do


              gdens(idime,igaus,1) = gdens(idime,igaus,1) + &
                   (elden(inode,1)-eldhy(inode)) * xcart
              gtemp(idime,igaus,1) = gtemp(idime,igaus,1) + &
                   (eltem(inode,1)-elthy(inode)) * xcart             
             ! gdens(idime,igaus,1) = gdens(idime,igaus,1) + &
             !      elden(inode,1) * xcart
             ! gtemp( idime,igaus,1) = gtemp(idime,igaus,1) + &
             !      eltem(inode,1) * xcart
              gpold(idime,igaus,1) = gpold(idime,igaus,1) + &
                   elpol(inode,1) * xcart
              gvisc(idime,igaus,1) = gvisc(idime,igaus,1) + &
                   elvis(inode,1) * xcart
              gvitu(idime,igaus,1) = gvitu(idime,igaus,1) + &
                   elvit(inode,1) * xcart
              gener(idime,igaus,1) = gener(idime,igaus,1) + &
                   elene(inode,1) * xcart
           end do
           
           !
           ! Kessler variables
           !
           if(kfl_benme_nsa >= 200) then
              do idime=1,ndime
                 xcart=cartd(idime,inode,igaus) 

                 gqvap(idime,igaus,1) = gqvap(idime,igaus,1) +&
                      elqva(inode,1) * xcart

                 gqclo(idime,igaus,1) = gqclo(idime,igaus,1) +&
                      elqcl(inode,1) * xcart

                 gqrai(idime,igaus,1) = gqrai(idime,igaus,1) +&
                      elqra(inode,1) * xcart
              end do
           end if


           if (kfl_relat_nsa == 1) then
              xenth(igaus,1) = xenth(igaus,1) + xshap*elent(inode,1) 
              xenth(igaus,2) = xenth(igaus,2) + xshap*elent(inode,2) 
              xenth(igaus,3) = xenth(igaus,3) + xshap*elent(inode,3) 
           end if

           if (kfl_turbu_nsa /= 0) then
              xkene(igaus,1) = xkene(igaus,1) + xshap*elken(inode,1) 
              xepsi(igaus,1) = xkene(igaus,1) + xshap*eleps(inode,1) 
              do idime=1,ndime
                 xcart=cartd(idime,inode,igaus) 
                 gkene(idime,igaus,1) = gkene(idime,igaus,1) + &
                      elken(inode,1) * xcart
                 gkevo(idime,igaus,1) = gkevo(idime,igaus,1) + &
                      elken(inode,1)/elden(inode,1) * xcart
                 gepsi(idime,igaus,1) = gepsi(idime,igaus,1) + &
                      eleps(inode,1) * xcart
              end do
           end if
        end do
        
     end do gauss_points_loop_2
     
!     tend_gacons_gaus_loop3 = cpu_time()-tstart_gacons_gaus_loop3
!     ttotal_gacons_gaus_loop3 = ttotal_gacons_gaus_loop3 + tend_gacons_gaus_loop3*69.84E-9


     gauss_points_loop3: do igaus=1,ngaus        
!!        do inode=1,pnode
          ! xshap = gpsha(inode,igaus)
           ! Pression on gauss points 
          ! xpres(igaus,1) = xpres(igaus,1) + elpre(inode,1) * xshap
           call nsa_stalaw(2,0,xdens(igaus,1),xpres(igaus,1),xtemp(igaus,1),xcoor(ndime,igaus),xmvlo,0.0_rp)
           do idime=1,ndime
             ! xcart=cartd(idime,inode,igaus) 
             ! xvelo(idime,igaus,1) = xvelo(idime,igaus,1) + elvel(idime,inode,1) * xshap
             ! xvelo(idime,igaus,2) = xvelo(idime,igaus,2) + elvel(idime,inode,2) * xshap
             ! xvelo(idime,igaus,3) = xvelo(idime,igaus,3) + elvel(idime,inode,3) * xshap
             ! gpres(idime,igaus,1) = gpres(idime,igaus,1) + elpre(inode,1) * xcart
             ! gphyd(idime,igaus)   = gphyd(idime,igaus) + elphy(inode) * xcart

              xvelo(idime,igaus,1) = xumom(idime,igaus,1) / xdens(igaus,1)
              xvelo(idime,igaus,2) = xumom(idime,igaus,2) / xdens(igaus,2)
              xvelo(idime,igaus,3) = xumom(idime,igaus,3) / xdens(igaus,3)           
              gpres(idime,igaus,1) = pbaro_nsa * adgam_nsa * (rgasc_nsa/pbaro_nsa) ** &
                   adgam_nsa * (xtemp(igaus,1) * xdens(igaus,1)) ** (adgam_nsa-1.0_rp) * &
                   (xtemp(igaus,1) * gdens(idime,igaus,1) + xdens(igaus,1) * gtemp(idime,igaus,1))
              do jdime=1,ndime
                ! gvelo(idime,jdime,igaus,1) = gvelo(idime,jdime,igaus,1)  + &
                !      elvel(jdime,inode,1) * xcart
                 gvelo(idime,jdime,igaus,1) = gumom(idime,jdime,igaus,1) / &
                      xdens(igaus,1) - xumom(jdime,igaus,1) / xdens(igaus,1) / &
                      xdens(igaus,1) * gdens(idime,igaus,1)
              end do
           end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! OJO: El xdiff(ndime+2,.,.,.) esta puesto a cero pero no tiene que ser asi!!!!!!!!!
  
     velsq = 0.0_rp
     do idime=1,ndime
        ! Total momentum unknown
        xtunk(idime) = xumom(idime,igaus,1) + xsube(idime,igaus,2)
        !Square velocity
        velsq = velsq + xvelo(idime,igaus,1) * xvelo(idime,igaus,1)
        
     end do
     ! Total density unknown
     xtunk(ndime+1) = xdens(igaus,1) + xsube(ndime+1,igaus,2)
     ! Total momentu temperature unknown
     xtunk(ndime+2) = xtemp(igaus,1) + xsube(ndime+2,igaus,2)
     
     visci = xvisc(igaus,1) / xtunk(ndime+1)
     do jdime=1,ndime        
        xconv(ndime+1,jdime  ,jdime,igaus)= 1.0_rp
        xconv(jdime  ,ndime+2,jdime,igaus)= afact_nsa*adgam_nsa*(xdens(igaus,1)*xtemp(igaus,1))**(adgam_nsa-1.0_rp)*xdens(igaus,1)
        xconv(jdime  ,ndime+1,jdime,igaus)= afact_nsa*adgam_nsa*(xdens(igaus,1)*xtemp(igaus,1))**(adgam_nsa-1.0_rp)*xtemp(igaus,1)
        xconv(ndime+2,ndime+1,jdime,igaus)= 0.0_rp
        xconv(ndime+2,ndime+2,jdime,igaus)= xvelo(jdime,igaus,1)
!!!!        xdiff(ndime+2,ndime+1,jdime,jdime,igaus)= 0.0_rp
!!!!        xdiff(ndime+2,ndime+2,jdime,jdime,igaus)= 0.0_rp
        do idime=1,ndime
           xconv(idime  ,idime  ,jdime,igaus)= xconv(idime,idime  ,jdime,igaus) + xvelo(jdime,igaus,1) 
           xconv(idime  ,jdime  ,jdime,igaus)= xconv(idime,jdime  ,jdime,igaus) + xvelo(idime,igaus,1)
           xconv(idime  ,ndime+1,jdime,igaus)= xconv(idime,ndime+1,jdime,igaus) - xvelo(idime,igaus,1) * xvelo(jdime,igaus,1)
           xconv(ndime+2,idime  ,jdime,igaus)= 0.0_rp
!!!!           xdiff(jdime,jdime,idime,idime,igaus)= visci
!!!!           xdiff(jdime,idime,idime,jdime,igaus)= xdiff(jdime,idime,idime,jdime,igaus) + visci
!!!!           xdiff(jdime,idime,jdime,idime,igaus)= xdiff(jdime,idime,jdime,idime,igaus) - 2.0_rp * visci / 3.0_rp
!!!!           xdiff(jdime,ndime+1,idime,idime,igaus)= - visci * xvelo(jdime,igaus,1)
!!!!           xdiff(jdime,ndime+1,idime,jdime,igaus)= xdiff(jdime,ndime+1,idime,jdime,igaus) - visci * xvelo(idime,igaus,1)
!!!!           xdiff(jdime,ndime+1,jdime,idime,igaus)= xdiff(jdime,ndime+1,jdime,idime,igaus) + &
!!!!                2.0_rp * visci * xvelo(idime,igaus,1) / 3.0_rp
!!!!           xdiff(ndime+2,idime,jdime,jdime,igaus)=  0.0_rp
!!!!           xdiff(ndime+2,jdime,jdime,idime,igaus)=  0.0_rp
!!!!           xdiff(ndime+2,jdime,idime,jdime,igaus)=  0.0_rp
!!!!           xdiff(ndime+2,ndime+1,jdime,idime,igaus)=  0.0_rp
        end do
     end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        end do
     end do gauss_points_loop3
     
  end if

end subroutine nsa_gacons


!-----------------------------------------------------------
! Like nsa_gacons above but it does not use the
! computation at the Gauss points because 
! this case does not have a reference background state
! given by an analytical function
!
! Simone Marras SM, March 2012
!-----------------------------------------------------------
subroutine nsa_gacons_kessler(&
     ndime,mnode,pnode,npoin,ngaus,mgaus,ncomp,ntens,&
     ivert,lnods,&
     kfl_foreg,kfl_turbu,kfl_relat,&
     kfl_hysta,kfl_inkee,kfl_inifi,&
     naxis,plapl,&
     elvel,&
     elpre,elphy,elpol,&
     elden,eldhy,&
     eltem,elthy,&
     elent,&
     elqva,elqva_hy,&
     elqcl,&
     elqra,&
     elvis,elvit,&
     elene,elumo,elfmo,elcod,&
     elken,eleps,elsub,&
     xvelo,&
     xpres,xphyd,&
     xdens,xdhyd,&
     xtemp,xthyd,&
     xenth,&
     xqvap,xqvap_hy,&
     xqclo,&
     xqrai,&     
     xvisc,xvitu,xener,xumom,xfmom,xsube,&
     dgvel,&
     gvelo,&
     gpres,gphyd,gpold,&
     gdens,gdhyd,&
     gtemp,gthyd,&
     gqvap,&
     gqclo,&
     gqrai,&
     gvisc,gvitu,gener,gumom,gfmom,gdvel,&
     xkene,xepsi,gkene,gepsi,gkevo,xlapr,xlake,xlate,xlade,xldve,xcoor,zeroc,&
     dvolu,cartd,hessi,heslo,gpsha,deriv,weigp,frano,gpcod,gpcor,wvect)
  !-----------------------------------------------------------------------
  !
  ! Gather operations and simultaneous gauss point evaluation
  !
  !-----------------------------------------------------------------------
  use      def_kintyp
  use      def_master, only : veloc,&
                              press,&
                              entha,&
                              densi,&
                              tempe,&
                              conce,&
                              visco,&
                              turmu,&
                              energ,&
                              umome,&
                              untur

  use      def_domain, only : coord
  use      def_nastal, only : kfl_infun_nsa, &
                              stapa_nsa,     &
                              brunt_nsa,     &
                              brure_nsa,     &
                              kfl_brunt_nsa, &
                              rekee_nsa,     &
                              kfl_benme_nsa
  
  implicit none
  integer(ip),  intent(in)  :: ndime,mnode,pnode,npoin,ngaus,mgaus,ncomp,ntens,plapl,naxis,ivert
  integer(ip),  intent(in)  :: lnods(pnode)
  integer(ip),  intent(in)  :: kfl_foreg,kfl_turbu,kfl_hysta,kfl_inkee,kfl_inifi,kfl_relat
  integer(ip)               :: inode,ipoin,idome,igaus,itidi,idime,jdime,itime,nprev

  real(rp), intent(in)  :: &
       gpsha(pnode,ngaus),deriv(ndime,pnode,ngaus),&
       weigp(ngaus),heslo(ntens,pnode,ngaus),frano(4,2),zeroc

  real(rp), intent(out) :: &
       elvel(ndime,mnode,ncomp),&
       elpre(mnode,ncomp),elphy(mnode),elpol(mnode,ncomp),&
       elden(mnode,ncomp),eldhy(mnode),&
       eltem(mnode,ncomp),elthy(mnode),&
       elent(mnode,ncomp),&
       elqva(mnode,ncomp),elqva_hy(mnode,ncomp),&
       elqcl(mnode,ncomp),&
       elqra(mnode,ncomp),&
       elvis(mnode,ncomp),elvit(mnode,ncomp),&
       elene(mnode,ncomp),elumo(ndime,mnode,ncomp),&
       elfmo(ndime,mnode,ncomp),elcod(ndime,mnode),&
       elken(mnode,ncomp),&
       eleps(mnode,ncomp),&
       xdens(mgaus,ncomp),xdhyd(mgaus),&
       xenth(mgaus,ncomp),&
       xpres(mgaus,ncomp),xphyd(mgaus),&
       xtemp(mgaus,ncomp),xthyd(mgaus),&
       xqvap(mgaus,ncomp),xqvap_hy(mgaus,ncomp),&
       xqclo(mgaus,ncomp),&
       xqrai(mgaus,ncomp),&
       xkene(mgaus,ncomp),&
       xvisc(mgaus,ncomp),xvitu(mgaus,ncomp),xener(mgaus,ncomp),xepsi(mgaus,ncomp),&
       xvelo(ndime,mgaus,ncomp),xfmom(ndime,mgaus,ncomp),xumom(ndime,mgaus,ncomp),&
       xlapr(mgaus,ncomp),xlate(mgaus,ncomp),xlake(mgaus,ncomp),xlade(mgaus,ncomp),&
       xldve(ndime,mgaus,ncomp),&
       gdens(ndime,mgaus,ncomp),gdhyd(ndime,mgaus),&
       gpres(ndime,mgaus,ncomp),gphyd(ndime,mgaus),gpold(ndime,mgaus,ncomp),&
       gtemp(ndime,mgaus,ncomp),gthyd(ndime,mgaus),&
       gqvap(ndime,mgaus,ncomp),gqclo(ndime,mgaus,ncomp),gqrai(ndime,mgaus,ncomp),&
       gener(ndime,mgaus,ncomp),gvisc(ndime,mgaus,ncomp),gvitu(ndime,mgaus,ncomp),&
       gkene(ndime,mgaus,ncomp),gepsi(ndime,mgaus,ncomp),xcoor(ndime,mgaus),&
       gvelo(ndime,ndime,mgaus,ncomp),gfmom(ndime,ndime,mgaus,ncomp),&
       gumom(ndime,ndime,mgaus,ncomp),gkevo(ndime,mgaus,ncomp),dvolu(mgaus),&
       cartd(ndime,mnode,mgaus),hessi(ntens,mnode,mgaus),gdvel(ndime,mgaus,ncomp),&
       dgvel(ndime,mgaus,ncomp),wvect(ndime,mgaus),gpcod(ndime,mgaus),gpcor(mgaus)

  real(rp)              :: detjm,xshap,xcart,xhess,xdelo,xtelo,xprlo,xcolo
  real(rp)              :: d2sdx(ndime,ndime,ndime)
  real(rp)              :: auxi1(ndime,ndime,mnode)
  real(rp)              :: auxi2(ndime,ndime,mnode)
  real(rp)              :: xjaci(ndime,ndime),xjacm(ndime,ndime) 
  real(rp)              :: elsub(ndime+2,mnode,2),xsube(ndime+2,mgaus,2),dummr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!miquel
!  integer*8 :: cpu_time
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!miquel

  !
  ! Gathering of the current values 
  !

!     tstart_gath = cpu_time()

  nprev= 4

  do inode= 1,pnode
     ipoin= lnods(inode)

     !
     !Element hydrostatic values 
     !(as they come from the initialization field or previous time-step):
     !
     elphy(inode) = rekee_nsa(ndime+3,ipoin)
     elthy(inode) = rekee_nsa(ndime+2,ipoin)
     eldhy(inode) = rekee_nsa(ndime+1,ipoin)

     !Density:
     elden(inode,1) = densi(ipoin,1)
     elden(inode,2) = densi(ipoin,2)
     elden(inode,3) = densi(ipoin,nprev)

     !Temperature:
     eltem(inode,1) = tempe(ipoin,1)
     eltem(inode,2) = tempe(ipoin,2)
     eltem(inode,3) = tempe(ipoin,nprev)
     
     !
     ! Kessler variables
     !
     if(kfl_benme_nsa >= 200) then
        elqva(inode,1) = conce(ipoin,1,1)
        elqva(inode,2) = conce(ipoin,1,2)
        elqva(inode,3) = conce(ipoin,1,3)

        elqcl(inode,1) = conce(ipoin,2,1)
        elqcl(inode,2) = conce(ipoin,2,2)
        elqcl(inode,3) = conce(ipoin,2,3)

        elqra(inode,1) = conce(ipoin,3,1)
        elqra(inode,2) = conce(ipoin,3,2)
        elqra(inode,3) = conce(ipoin,3,3)

        elqva_hy(inode,1) = rekee_nsa(ndime+3 + 1, ipoin)
     end if
     
     !Energy
     elene(inode,1) = energ(ipoin,1)        
     elene(inode,2) = energ(ipoin,2)        
     elene(inode,3) = energ(ipoin,nprev)        
     elvis(inode,1) = visco(ipoin,1)
     do idime=1,ndime
        elcod(idime,inode  ) = coord(idime,ipoin  )          
        elvel(idime,inode,1) = veloc(idime,ipoin,1)
        elvel(idime,inode,2) = veloc(idime,ipoin,2)
        elvel(idime,inode,3) = veloc(idime,ipoin,nprev)
        elumo(idime,inode,1) = umome(idime,ipoin,1)
        elumo(idime,inode,2) = umome(idime,ipoin,2)
        elumo(idime,inode,3) = umome(idime,ipoin,nprev)
        elfmo(idime,inode,1) = 0.0_rp
     end do
     elpol(        inode,1) = press(        ipoin,nprev)      ! old p, for fract. step  
     elpre(        inode,1) = press(        ipoin,1)          ! new p, for fract. step

     if (kfl_turbu /= 0) then
        elvit(inode,1) = turmu(  ipoin  )
        elken(inode,1) = untur(1,ipoin,1)
        eleps(inode,1) = untur(2,ipoin,1)
     else
        elvit(inode,1) = 0.0_rp
        elken(inode,1) = 0.0_rp
        eleps(inode,1) = 0.0_rp
     end if
  end do

!     tend_gath = cpu_time()-tstart_gath
!     ttotal_gath = ttotal_gath + tend_gath*69.84E-9

  xdelo= 0.0_rp
  xprlo= 0.0_rp
  xtelo= 0.0_rp
  
!     tstart_gacons_gaus_loop1  = cpu_time()

  do igaus=1,ngaus

!     tstart_gacons_elmder2_loop = cpu_time()

     call elmder(pnode,ndime,deriv(1,1,igaus),elcod,cartd(1,1,igaus),detjm,xjacm,xjaci)
     dvolu(igaus)=weigp(igaus)*detjm        

!     tend_gacons_elmder2_loop = cpu_time()-tstart_gacons_elmder2_loop
!     ttotal_gacons_elmder2_loop = ttotal_gacons_elmder2_loop + tend_gacons_elmder2_loop*69.84E-9
     
     hessi(1:ntens,1:mnode,igaus) = 0.0_rp
     if(plapl==1) then

!          tstart_gacons_elmhes2_loop = cpu_time()

          call elmhes(heslo(1,1,igaus),hessi(1,1,igaus),ndime,pnode,ntens,&
          xjaci,d2sdx,deriv(1,1,igaus),elcod)     

!          tend_gacons_elmhes2_loop = cpu_time()-tstart_gacons_elmhes2_loop
!          ttotal_gacons_elmhes2_loop = ttotal_gacons_elmhes2_loop + tend_gacons_elmhes2_loop*69.84E-9
       endif
     
     wvect(1:ndime,igaus)=0.0_rp
     if(frano(1,2)>zeroc.or.frano(2,2)>zeroc.or.naxis==1) then 
        call mbvab0(gpcod(1,igaus),elcod,gpsha(1,igaus),ndime,pnode)
        wvect(1:ndime,igaus)=gpcod(1:ndime,igaus)
        gpcor(igaus)=gpcod(1,igaus)                    
     end if
! end do

!    tend_gacons_gaus_loop1 = cpu_time()-tstart_gacons_gaus_loop1
!    ttotal_gacons_gaus_loop1 = ttotal_gacons_gaus_loop1 + tend_gacons_gaus_loop1*69.84E-9
  
!    tstart_gacons_gaus_loop2  = cpu_time()

! do igaus=1,ngaus
     ! Initialize gauss point values
     xdens(igaus,1)=0.0_rp 
     xdens(igaus,2)=0.0_rp 
     xdens(igaus,3)=0.0_rp
     xdhyd(igaus)  =0.0_rp 

     xpres(igaus,1)=0.0_rp 
     xpres(igaus,2)=0.0_rp 
     xpres(igaus,3)=0.0_rp 
     xphyd(igaus)  =0.0_rp 

     xtemp(igaus,1)=0.0_rp 
     xtemp(igaus,2)=0.0_rp 
     xtemp(igaus,3)=0.0_rp 
     xthyd(igaus)  =0.0_rp
     
     !
     ! Kessler variables
     !
     if(kfl_benme_nsa >= 200) then
        xqvap(igaus,1)=0.0_rp 

        xqclo(igaus,1)=0.0_rp 

        xqrai(igaus,1)=0.0_rp 
        
        xqvap_hy(igaus,1)=0.0_rp 
     end if
     
     xenth(igaus,1)=0.0_rp 
     xenth(igaus,2)=0.0_rp 
     xenth(igaus,3)=0.0_rp 
     xener(igaus,1)=0.0_rp 
     xener(igaus,2)=0.0_rp 
     xener(igaus,3)=0.0_rp 
     xvisc(igaus,1)=0.0_rp 
     xvitu(igaus,1)=0.0_rp 
     xkene(igaus,1)=0.0_rp 
     xepsi(igaus,1)=0.0_rp 
     xlapr(igaus,1)=0.0_rp 
     xlake(igaus,1)=0.0_rp 
     xlade(igaus,1)=0.0_rp 
     xlate(igaus,1)=0.0_rp 
     xcoor(1,igaus)=0.0_rp 
     xcoor(2,igaus)=0.0_rp 
     xcoor(ndime,igaus)=0.0_rp      
     do idime=1,ndime
        xvelo(idime,igaus,1)=0.0_rp 
        xldve(idime,igaus,1)=0.0_rp 
        xvelo(idime,igaus,2)=0.0_rp 
        xvelo(idime,igaus,3)=0.0_rp 
        xumom(idime,igaus,1)=0.0_rp 
        xumom(idime,igaus,2)=0.0_rp 
        xumom(idime,igaus,3)=0.0_rp 
        xfmom(idime,igaus,1)=0.0_rp 
        gdens(idime,igaus,1)=0.0_rp 
        gdhyd(idime,igaus)  =0.0_rp
        gpres(idime,igaus,1)=0.0_rp 
        gpres(idime,igaus,2)=0.0_rp 
        gpold(idime,igaus,1)=0.0_rp
        gphyd(idime,igaus)  =0.0_rp
        gtemp(idime,igaus,1)=0.0_rp
        gthyd(idime,igaus)  =0.0_rp 
        !
        ! Kessler variables
        !
        if(kfl_benme_nsa >= 200) then
           gqvap(idime,igaus,1)=0.0_rp 
           gqclo(idime,igaus,1)=0.0_rp 
           gqrai(idime,igaus,1)=0.0_rp 
        end if
        
        gener(idime,igaus,1)=0.0_rp 
        gvisc(idime,igaus,1)=0.0_rp 
        gvitu(idime,igaus,1)=0.0_rp 
        gepsi(idime,igaus,1)=0.0_rp 
        gkene(idime,igaus,1)=0.0_rp 
        gkevo(idime,igaus,1)=0.0_rp 
        dgvel(idime,igaus,1)=0.0_rp 
        gdvel(idime,igaus,1)=0.0_rp 
        do jdime=1,ndime
           gvelo(idime,jdime,igaus,1)=0.0_rp 
           gfmom(idime,jdime,igaus,1)=0.0_rp 
           gumom(idime,jdime,igaus,1)=0.0_rp 
        end do
     end do
     
     ! Compute gauss values for press and veloc, which can be done now
     do inode=1,pnode
        xshap = gpsha(inode,igaus)
        
        xvelo(1:ndime,igaus,1) = xvelo(1:ndime,igaus,1) + elvel(1:ndime,inode,1)*xshap
        xvelo(1:ndime,igaus,2) = xvelo(1:ndime,igaus,2) + elvel(1:ndime,inode,2)*xshap
        xvelo(1:ndime,igaus,3) = xvelo(1:ndime,igaus,3) + elvel(1:ndime,inode,3)*xshap
        xcoor(1:ndime,igaus  ) = xcoor(1:ndime,igaus)   + elcod(1:ndime,inode)*xshap
        do idime=1,ndime
           xcart=cartd(idime,inode,igaus) 
           gpres(idime,igaus,1) = gpres(idime,igaus,1) +                 & 
                elpre(        inode,1) * xcart
           gphyd(idime,igaus)   = gphyd(idime,igaus) +                   & ! hydrostatic gpres
                elphy(        inode) * xcart
           gpold(idime,igaus,1) = gpold(idime,igaus,1) +                 & 
                elpol(        inode,1) * xcart
           xfmom(        idime,igaus,1) = xfmom(        idime,igaus,1) + &
                elfmo(  idime,inode,1) * xshap
           xumom(        idime,igaus,1) = xumom(        idime,igaus,1) + &
                elumo(  idime,inode,1) * xshap
           xumom(        idime,igaus,2) = xumom(        idime,igaus,2) + &
                elumo(  idime,inode,2) * xshap
           xumom(        idime,igaus,3) = xumom(        idime,igaus,3) + &
                elumo(  idime,inode,3) * xshap
           do jdime=1,ndime
              gfmom(idime,jdime,igaus,1) = gfmom(idime,jdime,igaus,1)  + &
                   elfmo(jdime,inode,1) * xcart
              gumom(idime,jdime,igaus,1) = gumom(idime,jdime,igaus,1)  + &
                   elumo(jdime,inode,1) * xcart
              gvelo(idime,jdime,igaus,1) = gvelo(idime,jdime,igaus,1)  + &
                   elvel(jdime,inode,1) * xcart
           end do
        end do
        if (plapl==1) then
           if(ndime.eq.2) then
              gdvel(1,igaus,1)=gdvel(1,igaus,1)+hessi(1,inode,igaus)*elvel(1,inode,1) &
                   +hessi(3,inode,igaus)*elvel(2,inode,1)
              gdvel(2,igaus,1)=gdvel(2,igaus,1)+hessi(3,inode,igaus)*elvel(1,inode,1) &
                   +hessi(2,inode,igaus)*elvel(2,inode,1)
           else if(ndime.eq.3) then
              gdvel(1,igaus,1)=gdvel(1,igaus,1)+hessi(1,inode,igaus)*elvel(1,inode,1) & 
                   +hessi(4,inode,igaus)*elvel(2,inode,1) &
                   +hessi(5,inode,igaus)*elvel(3,inode,1)
              gdvel(2,igaus,1)=gdvel(2,igaus,1)+hessi(4,inode,igaus)*elvel(1,inode,1) & 
                   +hessi(2,inode,igaus)*elvel(2,inode,1) &
                   +hessi(6,inode,igaus)*elvel(3,inode,1)
              gdvel(3,igaus,1)=gdvel(3,igaus,1)+hessi(5,inode,igaus)*elvel(1,inode,1) & 
                   +hessi(6,inode,igaus)*elvel(2,inode,1) &
                   +hessi(3,inode,igaus)*elvel(3,inode,1)
           end if
           do idime=1,ndime
              xhess = hessi(idime,inode,igaus)
              xlapr(igaus,1) = xlapr(igaus,1) + elpre(inode,1)*xhess
              xlade(igaus,1) = xlade(igaus,1) + elden(inode,1)*xhess
              xlate(igaus,1) = xlate(igaus,1) + eltem(inode,1)*xhess
              xlake(igaus,1) = xlake(igaus,1) + elken(inode,1)*xhess
              xldve(idime,igaus,1) = xldve(idime,igaus,1) + elden(inode,1)*xhess
              do jdime=1,ndime
                 dgvel(idime,igaus,1)=dgvel(idime,igaus,1) + &
                      hessi(idime,inode,igaus)*elvel(jdime,inode,1)
              end do
           end do
        end if
     end do
  end do

!     tend_gacons_gaus_loop1 = cpu_time()-tstart_gacons_gaus_loop1
!     ttotal_gacons_gaus_loop1 = ttotal_gacons_gaus_loop1 + tend_gacons_gaus_loop1*69.84E-9

  
  if (kfl_foreg == 0) then     

!     tstart_gacons_gaus_loop3  = cpu_time()

!------- BUCLE ORIGINAL -------------miquel
!!$     do inode=1,pnode
!!$        ! Gauss point values
!!$        do igaus=1,ngaus
!------- BUCLE OPTIMITZAT -------------miquel
        ! Gauss point values
     do igaus=1,ngaus
        do inode=1,pnode


           xshap = gpsha(inode,igaus)
           
           !
           !   xvmas(inode,igaus) = xshap * elmas(inode)
           !
           !Densi:
           xdens(igaus,1) = xdens(igaus,1) + xshap*elden(inode,1) 
           xdens(igaus,2) = xdens(igaus,2) + xshap*elden(inode,2) 
           xdens(igaus,3) = xdens(igaus,3) + xshap*elden(inode,3) 
           !Densi hydro:
           xdhyd(igaus  ) = xdhyd(igaus  ) + xshap*eldhy(inode  ) 
 
           !Tempe:
           xtemp(igaus,1) = xtemp(igaus,1) + xshap*eltem(inode,1) 
           xtemp(igaus,2) = xtemp(igaus,2) + xshap*eltem(inode,2) 
           xtemp(igaus,3) = xtemp(igaus,3) + xshap*eltem(inode,3)
           xthyd(igaus  ) = xthyd(igaus  ) + xshap*elthy(inode  ) 
           
           !Press
           xpres(igaus,1) = xpres(igaus,1) + elpre(inode,1)*xshap
           !Press hydro:
           xphyd(igaus  ) = xphyd(igaus  ) + xshap*elphy(inode  ) 

           !
           ! Kessler variables
           !
           if(kfl_benme_nsa >= 200) then
              xqvap(igaus,1) = xqvap(igaus,1) + xshap*elqva(inode,1)
              !xqvap(igaus,2) = xqvap(igaus,2) + xshap*elqva(inode,2)
              !xqvap(igaus,3) = xqvap(igaus,3) + xshap*elqva(inode,3) 

              xqclo(igaus,1) = xqclo(igaus,1) + xshap*elqcl(inode,1) 
              !xqclo(igaus,2) = xqclo(igaus,2) + xshap*elqcl(inode,2)
              !xqclo(igaus,3) = xqclo(igaus,3) + xshap*elqcl(inode,3)

              xqrai(igaus,1) = xqrai(igaus,1) + xshap*elqra(inode,1) 
              !xqrai(igaus,2) = xqrai(igaus,2) + xshap*elqra(inode,2) 
              !xqrai(igaus,3) = xqrai(igaus,3) + xshap*elqra(inode,3) 

              xqvap_hy(igaus,1) = xqvap_hy(igaus,1) + xshap*elqva_hy(inode,1) 
           end if

           xener(igaus,1) = xener(igaus,1) + xshap*elene(inode,1) 
           xener(igaus,2) = xener(igaus,2) + xshap*elene(inode,2) 
           xener(igaus,3) = xener(igaus,3) + xshap*elene(inode,3) 
           xvisc(igaus,1) = xvisc(igaus,1) + xshap*elvis(inode,1) 
           xvitu(igaus,1) = xvitu(igaus,1) + xshap*elvit(inode,1) 
           
           do idime=1,ndime
              xcart=cartd(idime,inode,igaus) 
              gdens(        idime,igaus,1) = gdens(        idime,igaus,1) + &
                   elden(        inode,1) * xcart
              gtemp(        idime,igaus,1) = gtemp(        idime,igaus,1) + &
                   eltem(        inode,1) * xcart
              gvisc(        idime,igaus,1) = gvisc(        idime,igaus,1) + &
                   elvis(        inode,1) * xcart
              gvitu(        idime,igaus,1) = gvitu(        idime,igaus,1) + &
                   elvit(        inode,1) * xcart
              gener(        idime,igaus,1) = gener(        idime,igaus,1) + &
                   elene(        inode,1) * xcart
           end do
           !
           ! Kessler variables
           !
           if(kfl_benme_nsa >= 200) then
              do idime=1,ndime
                 xcart=cartd(idime,inode,igaus) 

                 gqvap(idime,igaus,1) = gqvap(idime,igaus,1) +&
                      elqva(inode,1) * xcart

                 gqclo(idime,igaus,1) = gqclo(idime,igaus,1) +&
                      elqcl(inode,1) * xcart

                 gqrai(idime,igaus,1) = gqrai(idime,igaus,1) +&
                      elqra(inode,1) * xcart
              end do
           end if

           if (kfl_turbu /= 0) then
              xkene(igaus,1) = xkene(igaus,1) + xshap*elken(inode,1) 
              xepsi(igaus,1) = xkene(igaus,1) + xshap*eleps(inode,1) 
              do idime=1,ndime
                 xcart=cartd(idime,inode,igaus) 
                 gkene(idime,igaus,1) = gkene(idime,igaus,1) + &
                      elken(inode,1) * xcart
                 gkevo(idime,igaus,1) = gkevo(idime,igaus,1) + &
                      elken(inode,1)/elden(inode,1) * xcart
                 gepsi(idime,igaus,1) = gepsi(idime,igaus,1) + &
                      eleps(inode,1) * xcart
              end do
           end if
        end do
        
     end do
     
!     tend_gacons_gaus_loop3 = cpu_time()-tstart_gacons_gaus_loop3
!     ttotal_gacons_gaus_loop3 = ttotal_gacons_gaus_loop3 + tend_gacons_gaus_loop3*69.84E-9


  end if
  
end subroutine nsa_gacons_kessler
