!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_elchea.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Mixed conservative and heat (rho,U,T) equations per-element operations:
!> @details Mixed conservative and heat (rho,U,T) equations per-element operations:
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_elchea
  use      def_master
  use      def_domain
  use      def_nastal

  implicit none

  real(rp)    :: elrhs(nevat_nsa)
  real(rp)    :: elsta(nevat_nsa)  !Element stabilization term (for postprocessing only)
  real(rp)    :: eldia(nevat_nsa)

  real(rp)    :: &
       elvel(ndime,mnode,ncomp_nsa),elcod(ndime,mnode),elene(mnode,ncomp_nsa),&
       elumo(ndime,mnode,ncomp_nsa),elfmo(ndime,mnode,ncomp_nsa),elpre(mnode,ncomp_nsa),&
       elpol(mnode,ncomp_nsa),elden(mnode,ncomp_nsa),eltem(mnode,ncomp_nsa),&
       eldhy(mnode),elphy(mnode),elthy(mnode),&  !Hydrostatic element density, pressure, and theta
       elqva(mnode,ncomp_nsa),elqva_hy(mnode,ncomp_nsa),elqcl(mnode,ncomp_nsa),elqra(mnode,ncomp_nsa),&
       elvis(mnode,ncomp_nsa),&
       elken(mnode,ncomp_nsa),eleps(mnode,ncomp_nsa),elvit(mnode,ncomp_nsa),elent(mnode,ncomp_nsa),&
       elsub(ndofn_nsa,mnode,2),elcon(ndofn_nsa,ndofn_nsa,ndime,mnode), &
       eldif(ndofn_nsa,ndofn_nsa,ndime,ndime,mnode)
  

  ! These are defined per-gauss point and dimensioned (... , mgaus)
  ! 
  ! -------------->>>>>>>>>>

  real(rp)    :: &
       xdens(mgaus,ncomp_nsa),xdhyd(mgaus),xphyd(mgaus),xthyd(mgaus),xpres(mgaus,ncomp_nsa),xtemp(mgaus,ncomp_nsa),&
       xqvap(mgaus,ncomp_nsa),xqvap_hy(mgaus,ncomp_nsa),xqclo(mgaus,ncomp_nsa),xqrai(mgaus,ncomp_nsa), &
       xconc(mgaus,ncomp_nsa), &
       xvisc(mgaus,ncomp_nsa),xvisc2(mgaus),&
       xener(mgaus,ncomp_nsa),xvelo(ndime,mgaus,ncomp_nsa),xenth(mgaus,ncomp_nsa),&
       xfmom(ndime,mgaus,ncomp_nsa),xumom(ndime,mgaus,ncomp_nsa),gdens(ndime,mgaus,ncomp_nsa),gdhyd(ndime,mgaus),&
       gpres(ndime,mgaus,ncomp_nsa),gpold(ndime,mgaus,ncomp_nsa),gphyd(ndime,mgaus),&
       gtemp(ndime,mgaus,ncomp_nsa),gthyd(ndime,mgaus),gkene(ndime,mgaus,ncomp_nsa),&
       gqvap(ndime,mgaus,ncomp_nsa),gqclo(ndime,mgaus,ncomp_nsa),gqrai(ndime,mgaus,ncomp_nsa),&
       gepsi(ndime,mgaus,ncomp_nsa),xvitu(mgaus,ncomp_nsa),gvitu(ndime,mgaus,ncomp_nsa),&
       gener(ndime,mgaus,ncomp_nsa),gvisc(ndime,mgaus,ncomp_nsa),gdvel(ndime,mgaus,ncomp_nsa),&
       dgvel(ndime,mgaus,ncomp_nsa),gvelo(ndime,ndime,mgaus,ncomp_nsa),xldve(ndime,mgaus,ncomp_nsa),&
       gfmom(ndime,ndime,mgaus,ncomp_nsa),gumom(ndime,ndime,mgaus,ncomp_nsa),difeq(ndofn_nsa)

  real(rp)    :: &
       xkene(mgaus,ncomp_nsa),xepsi(mgaus,ncomp_nsa),gkevo(ndime,mgaus,ncomp_nsa),&
       xlapr(mgaus,ncomp_nsa),xlate(mgaus,ncomp_nsa),xlade(mgaus,ncomp_nsa),xlake(mgaus,ncomp_nsa),xcoor(ndime,mgaus),&
       hessi(ntens,mnode,mgaus),cartd(ndime,mnode,mgaus),xjacm(ndime,ndime,mgaus),&
       xvcar(mnode,mgaus,ncomp_nsa),hconv,htran,hleti,velch,xmofr,vofor(ndime),dimst(3,3),taray(ndofn_nsa,ndofn_nsa,mgaus)

  real(rp)    :: &
       xconv(ndofn_nsa,ndofn_nsa,ndime,mgaus),gunkn(ndofn_nsa,ndime,mgaus),grasc(ndofn_nsa,ndime), &
       gsube(ndofn_nsa,ndime,mgaus), taudi(ndofn_nsa), &
       dconv(ndofn_nsa,ndofn_nsa,mgaus),xsube(ndofn_nsa,mgaus,3),xdtix(ndofn_nsa,mgaus,2), &
       xdiff(ndofn_nsa,ndofn_nsa,ndime,ndime,mgaus),ddiff(ndofn_nsa,ndofn_nsa,2,ndime,mgaus),xtime(ndofn_nsa,mgaus)

  real(rp)    :: &
       elunk(ndofn_nsa,mnode,ncomp_nsa),&
       conme(ndofn_nsa,ndofn_nsa,ndime),&
       difme(ndofn_nsa,ndofn_nsa,ndime,ndime),eldtt(ndofn_nsa,mnode,2),kapsh(ndofn_nsa),&
       xresi(ndofn_nsa,mgaus), xdith(mgaus),&
       xunkn(ndofn_nsa,mgaus,3)

  real(rp)    :: xshai,asfac,dtaux
  real(rp)    :: dvolu(mgaus)

  ! <<<<<<<<<<<<<-----------

  integer(ip) :: pnode,pgaus,plapl,itott
  integer(ip) :: pelty,ifsho,kfsho,ifcdr,pface,pblty,pnodb,pgaub

  integer(ip) :: iposi,jposi,imasl,ievat,jevat,kevat,ispec              ! Assembly
  integer(ip) :: idofn,itotv,ievab,ievav,ievac,ievae,jevav
  real(rp)    :: acont,adiag,totco 

  real(rp)    :: denac,visac,detjm,dvoti,gvcar,gtcar,cogra,cogri
  real(rp)    :: galte(ndofn_nsa),state(ndofn_nsa),galim(ndofn_nsa),staim(ndofn_nsa)
  real(rp)    :: shote(ndofn_nsa),resdc(ndofn_nsa),resta(ndofn_nsa),wvect(ndime,mgaus),&
       gpcod(ndime),waxis(mgaus)

  real(rp)    :: gpcor(3),vegau(ndime,3)

  real(rp)    :: shapw(mnode),deriw(ndime,mnode)         ! Center of gravity
  real(rp)    :: heslw(ntens,mnode)

  real(rp)    :: tragl(ndime,ndime),chave(ndime,2)       ! Stabilization
  real(rp)    :: chale(2),hleng(ndime)
  real(rp)    :: dshve,vidis,dundt(ndofn_nsa)
  real(rp)    :: ptesm(nevab_nsa,ndime),ptesd(nevab_nsa)
  real(rp)    :: rmom1(ndime,nevab_nsa)
  real(rp)    :: rmom2(ndime,nevab_nsa),bmome(ndime)
  real(rp)    :: rmom3(nevab_nsa),rpre1(ndime,mnode)

  integer(ip) :: i                           ! to define
  real(rp)    :: porac,frano(4,2)
  real(rp)    :: gshoc(3),condt,conco

  integer(ip) :: ielem,inode,inodb,ipoin,kpoin,igaus,igaub,idime         ! Indices and dimensions
  integer(ip) :: jdime,kdime,kdivi,kexvi,khior,ndofe,ndofc,iboun
  integer(ip) :: &
       iheat,imome,kmome,imrhs,imasr,imomr,kmomr,ihear,nfabo,ifabo,&
       ifibo,ibobo,ivars,kfoun
  real(rp)    :: ugrau(ndime),cterm(ndofn_nsa),vterm(ndofn_nsa),adiff(ndofn_nsa),spare(ndofn_nsa)
  real(rp)    :: tensi(ndime),tentu(ndime),tehea(ndime),vehea(ndime),shore(ndofn_nsa)
  real(rp)    :: baloc(ndime,ndime), bocod(ndime,mnodb),cartb(ndime,mnode)
  real(rp)    :: botem(mnodb,ncomp_nsa)
  real(rp)    :: banor(ndime),trapr(ndime)
  real(rp)    :: eucta,shapp(mnode),tauxi(3,5),dsurf

  real(rp)    :: sound,velmo,dvelo,cvaux,vitot,ditot,reate,tiint,tiine,tiinc,tiext,taupa,tauen,tauco
  real(rp)    :: vegte,gdgte,velix,xdila,xditu,gdito,sigmk
  real(rp)    :: dvene,cvene,prene,hoene,xshap,xcart,xauxi(2)
  real(rp)    :: fpert,qufac,dumom,dfmom,the12,thpre,thet2,thet1,theco,dotre,dvite,vweak,tweak,xdtau

  real(rp)    :: eplop,velop,vhemo,veige,vedim,seige,velvi,veaux,alpre,bepre,&
       gamso,xdenp,xdent,gpper

  !
  ! Kessler auxiliary variables:
  !
  real(rp)    :: qt_k,eps,es,rvs
  !
  ! Local auiliary variables
  !
  real(rp)    :: xdtot,xttot,xptot
  real(rp)    :: xdper,xtper,xpper
  real(rp)    :: gdtot,gttot,gptot
  
  integer(ip)::  icol, iz, nz, irestart, iloop, ifnp, j, istep
  character               :: fnp1*3, fnp*72

#ifdef _OPENMP
  INTEGER     :: OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
  INTEGER     :: thrid, numthr
  real(rp)    :: sufro_nsaOMP(16,4), fbody_nsaOMP(4,4,4)
#endif

  irestart = vtkrestart_nsa

  ccoef_nsa = 0.0_rp
  fbody_nsa = 0.0_rp
  surto_nsa = 0.0_rp
  sufro_nsa = 0.0_rp
  ifsho=0
  if (kfl_shock_nsa(1) > 0) ifsho=1
  if (kfl_shock_nsa(2) > 0) ifsho=1
  if (kfl_shock_nsa(3) > 0) ifsho=1

  thet1 = theta_nsa(1)
  thet2 = theta_nsa(2)
  the12 = thet1  * thet2
  thpre = thet1  - the12
  theco = 1.0_rp - thet2

  ndofc = ndime+1
  ndofe = ndime+2

  ! are there viscous terms around?
  kexvi = kfl_turbu_nsa + kfl_visco_nsa
  kdivi = kfl_turbu_nsa + lawvi_nsa

  dotre = 2.0_rp / 3.0_rp

  ! volume forces
  vofor = 0.0_rp

#ifdef _OPENMP
  sufro_nsaOMP = 0.0_rp
  fbody_nsaOMP = 0.0_rp
#endif

  !$*OMP  PARALLEL DO SCHEDULE (GUIDED)   &
  !$*OMP  DEFAULT (NONE)                 &
  !$*OMP  PRIVATE (alpre, baloc, banor, bepre, bocod, cartb, cartd, chale,  &
  !$*OMP           chave, cogra, cogri, cterm, cvene, dfmom, dgvel, ditot,  &
  !$*OMP           dshve, dsurf, dumom, dundt, dvelo, cvaux, dvene, dvite, dvolu,  &
  !$*OMP           elcod, elden, elent, eldia, elene, eleps, elfmo, elken, elpol,  &
  !$*OMP           elpre, elrhs, elsta, eltem, elqva, elqva_hy, elqcl, elqra,      $
  !$*OMP           elumo, elvel, elvis, elvit, eplop,  &
  !$*OMP           eucta, fpert, frano, galim, galte, gamso, gdens, gdhyd, gdgte,  &
  !$*OMP           gdito, gdvel, gener, gepsi, gfmom, gkene, gkevo, gpcod,  &
  !$*OMP           gpcor, gphyd, gpold, gpres, gpper, gshoc, condt, conco, &
  !$*OMP           gtemp, gthyd, gumom, gvcar, gtcar, gvelo,  &
  !$*OMP           gvisc, gvitu, hconv, hessi, hleng, hleti, hoene, htran,  &
  !$*OMP           iboun, idime, idofn, ielem, ievac, ievae, ievat, ievav,  &
  !$*OMP           ifabo, ifcdr, ifibo, igaub, igaus, inodb, inode, ipoin, ispec,  &
  !$*OMP           jdime, jevat, kfsho, khior, nfabo, pblty, pelty, pface,  &
  !$*OMP           pgaub, pgaus, plapl, pnodb, pnode, prene, qufac, reate,  &
  !$*OMP           resdc, resta, seige, shapp, shote, spare, staim, state,  &
  !$*OMP           tauco, tauen, taupa, tauxi, tehea, tensi, tentu, tiext,  &
  !$*OMP           tiinc, tiine, tiint, tragl, trapr, ugrau, veaux, vegte, &
  !$*OMP           vehea, veige, vedim, velch, velix, velmo, velop, velvi, vhemo,  &
  !$*OMP           vidis, vitot, vofor, dimst, vterm, adiff, vweak, tweak, wvect, xcart, xcoor, xdenp,  &
  !$*OMP           xdens, xenth,xdent, xdhyd, xphyd, xdila, xditu, xdtau, xener, xepsi, xfmom,  &
  !$*OMP           xkene, xlake, xlapr, xlate, xlade, xmofr, xpres, xshap, xtemp,xconc,elcon,eldif, &
  !$*OMP           xumom, xvcar, xvelo, xvisc, xvitu, thrid, xldve, taray, adiag, totco, xconv, xdiff, dconv, ddiff) & 
  !$*OMP   SHARED (adgam_nsa, afact_nsa, coord, cpcoe_nsa, cppra_nsa, cpprt_nsa,       & 
  !$*OMP           cvcoe_nsa, densi, dotre, dtieq_nsa, elmar, energ,        &
  !$*OMP           fbody_nsaOMP, gravi_nsa, grnor_nsa, hnatu, ifsho, ivert_nsa, &
  !$*OMP           kdivi, kexvi, kfl_advec_nsa, kfl_diagi_nsa, kfl_fixno_nsa, &
  !$*OMP           kfl_foreg_nsa, kfl_genal_nsa, kfl_hconv_nsa, kfl_hysta_nsa, kfl_lopre_nsa,&
  !$*OMP           kfl_naxis, kfl_shock_nsa, kfl_stabi_nsa, kfl_stafl_nsa, kfl_taufa_nsa,  &
  !$*OMP           kfl_turbu_nsa, kfl_visco_nsa, lawvi_nsa,  &
  !$*OMP           tempe_nsa, densi_nsa, rmach_nsa, speed_nsa, lawst_nsa, lboel, leobl_nsa, llapl, lnodb, lnods, ltypb, &
  !$*OMP           ltype, mgaus, mnode, nboun, ncomp_nsa, ndimb, ndime,     &
  !$*OMP           ndofn_nsa, ndofc, ndofe, nelem, nface, ngaus, nnode, npoin, &
  !$*OMP           ntens, press, rgasc_nsa, rhsid, safet_nsa, shock_nsa, sofac_nsa,    &
  !$*OMP           sound, sufro_nsaOMP, tempe, theco, thet1, theta_nsa,     &
  !$*OMP           thpre, turmu, umome, untur, vdiag_nsa, veloc, visco)     &
  !$*OMP   REDUCTION (+:surto_nsa)

  ccoef_nsa = 0.0_rp
  fbody_nsa = 0.0_rp
  surto_nsa     = 0.0_rp
  sufro_nsa     = 0.0_rp

  do ipoin = 1,npoin
     umoss_nsa(    1,ipoin,2) = 0.0_rp
     umoss_nsa(    2,ipoin,2) = 0.0_rp
     umoss_nsa(ndime,ipoin,2) = 0.0_rp
     denss_nsa(      ipoin,2) = 0.0_rp
     eness_nsa(      ipoin,2) = 0.0_rp

     if (kfl_diagi_nsa > 0) then
        do idime=1,ndime
           itott= (ipoin-1) * ndofn_nsa + idime
           vdiag_nsa(itott)= 0.0_rp
        end do
        vdiag_nsa(itott+1) = 0.0_rp
        vdiag_nsa(itott+2) = 0.0_rp     
     else if (kfl_lotim_nsa > 0) then
        do idime=1,ndime
           itott= (ipoin-1) * ndofn_nsa + idime
           dtaux= 1.0_rp/dtinv_nsa 
           if (dtieq_nsa(1,ipoin,1) < 1.2*dtaux) dtaux= dtieq_nsa(1,ipoin,1) 
           vdiag_nsa(itott) = dtaux
        end do
        dtaux= 1.0_rp/dtinv_nsa 
        if (dtieq_nsa(1,ipoin,1) < 1.2*dtaux) dtaux= dtieq_nsa(2,ipoin,1) 
        vdiag_nsa(itott+1) = dtaux
        dtaux= 1.0_rp/dtinv_nsa
        if (dtieq_nsa(1,ipoin,1) < 1.2*dtaux) dtaux= dtieq_nsa(3,ipoin,1) 
        vdiag_nsa(itott+2) = dtaux
     end if

  end do

  elements_loop: do ielem=1,nelem

#ifdef _OPENMP
     thrid = OMP_GET_THREAD_NUM() + 1
#endif

     ! Element properties and dimensions
     pelty = ltype(ielem)
     pnode = nnode(pelty)
     pgaus = ngaus(pelty)
     plapl = llapl(pelty)
     pface = nface(pelty)
     
     ! Initialize
     elrhs = 0.0_rp
     elsta = 0.0_rp
     eldia = 0.0_rp
     
     ! OJO: TO DEFINE
     frano = 0.0_rp
     khior = 0
     qufac = 1.0_rp
     if((ndime.eq.2).and.(pnode.ge.4)) then
        qufac = 0.5_rp 
        khior = 1
     end if
     if((ndime.eq.3).and.(pnode.ge.5))then
        qufac = 0.5_rp 
        khior = 1
     end if
     if (kfl_naxis==1) khior = 1

     !
     ! Gather operations and gauss point evaluation
     !
     if(kfl_benme_nsa < 200 .or. kfl_benme_nsa == 203 .or. kfl_benme_nsa == 204 .or.&
          kfl_benme_nsa == 210) then
        call nsa_gacons( &
             ndime,mnode,pnode,npoin,pgaus,mgaus,ncomp_nsa,ntens,ivert_nsa,lnods(1,ielem), &
             plapl,elvel,elpre,elpol,elden,elent,eltem,elvis,elvit, &
             elene,elumo,elfmo,elcod,elken,eleps,elsub, &  
             elqva,elqva_hy,&
             elqcl,&
             elqra,&
             xvelo,xpres,xphyd,xdens,xenth,xdhyd,xtemp,xthyd,xconc,xvisc,xvitu,xener,xumom,xfmom,xsube,dgvel, &
             xqvap,xqvap_hy,&
             xqclo,&
             xqrai,&
             gvelo,gpres,gpold,gphyd,gdens,gdhyd,gtemp,gthyd,gvisc,gvitu,&
             gqvap,&
             gqclo,&
             gqrai,&
             gener,gumom,gfmom,gdvel, &
             xkene,xepsi,gkene,gepsi,gkevo,xlapr,xlake,xlate,xlade,xldve,xcoor,zensa, &
             dvolu,cartd,hessi,elmar(pelty)%heslo,elmar(pelty)%shape,elmar(pelty)%deriv, &
             elmar(pelty)%weigp,frano,gpcod,gpcor,wvect,taray,elcon,eldif,xconv,dconv, &
             xdiff,ddiff)
     else
        !
        ! Gacons for problems where the background state
        ! is not defined by an analytic function and all
        ! quantities must be computed in the sense of FE (by interpolation) 
        ! Simone Marras SM, March 2012
        !
        call nsa_gacons_kessler(&
             ndime,mnode,pnode,npoin,pgaus,mgaus,ncomp_nsa,ntens,&
             ivert_nsa,lnods(1,ielem),&
             kfl_foreg_nsa,kfl_turbu_nsa,kfl_relat_nsa,&
             kfl_hysta_nsa,kfl_inkee_nsa,kfl_inifi_nsa(2),kfl_naxis,plapl,&
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
             xkene,xepsi,gkene,gepsi,gkevo,xlapr,xlake,xlate,xlade,xldve,xcoor,zensa,&
             dvolu,cartd,hessi,&
             elmar(pelty)%heslo,elmar(pelty)%shape,elmar(pelty)%deriv,&
             elmar(pelty)%weigp,frano,gpcod,gpcor,wvect)
     end if

     !
     ! Supplementary values
     !
     xvcar = 0.0_rp
     chave = 0.0_rp
     velch = 0.0_rp
     do idime=1,ndime
        do igaus=1,pgaus
           velix = xvelo(idime,igaus,1) 
           do inode=1,pnode
              xcart = cartd(idime,inode,igaus)
              xvcar(inode,igaus,1) = xvcar(inode,igaus,1) + velix * xcart
           end do
           chave(idime,1) = chave(idime,1) + velix / real(pgaus)           
        end do
        velch = velch + chave(idime,1)*chave(idime,1)
     end do

     ! hleng and tragl at center of gravity
     call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)

     htran = hleng(    1)
     hconv = hleng(ndime)

     ! Compute the convective characteristic length chale (MEAN, ELEMENT-WISE value!) 

     if (kfl_hconv_nsa == 1) then
        if (velch > zensa) then
           call nsa_elmchl(tragl,hleng,chave,chale,ndime,&
                pnode,hnatu(pelty),kfl_advec_nsa,zensa)
           hconv = chale(1)
        end if
     end if
     
     elemental_gauss_points: do igaus=1,pgaus

        velmo = 0.0_rp
        dvelo = 0.0_rp
        dumom = 0.0_rp
        dfmom = 0.0_rp
        vegte = 0.0_rp
        gdgte = 0.0_rp

        !Some auxiliary variables:
        !xdtot = xdens(igaus,1) + kfl_pertu_nsa*xdhyd(igaus)
        !xttot = xtemp(igaus,1) + kfl_pertu_nsa*xthyd(igaus)
        !xptot = xpres(igaus,1) + kfl_pertu_nsa*xphyd(igaus)
        xdtot = xdens(igaus,1) + xdhyd(igaus)
        xttot = xtemp(igaus,1) + xthyd(igaus)
        xptot = xpres(igaus,1) + xphyd(igaus)

        dvite = 0.0_rp

        if (lawvi_nsa > 0) call nsa_lawvis( -1,1,xvisc(igaus,1),xtemp(igaus,1),dvite)
        
        do idime=1,ndime
           velix = xvelo(idime,igaus,1)
           tehea(idime) = 0.0_rp
          ! dumom=dumom+gumom(idime,idime,igaus,1)
           dumom = dumom+velix*gdens(idime,igaus,1)+xdens(igaus,1)*gvelo(idime,idime,igaus,1)
           dfmom = dfmom +gfmom(idime,idime,igaus,1)
           dvelo = dvelo+gvelo(idime,idime,igaus,1)
           vegte = vegte+gtemp(idime,igaus,1)*velix
           gdgte = gdgte+gtemp(idime,igaus,1)*gdens(idime,igaus,1)
           velmo = velmo+velix*velix

           gvisc(idime,igaus,1) = gvisc(idime,igaus,1) * dvite
           vehea(idime) = velix
        end do
        velmo = sqrt(velmo)
        vhemo = velmo

        vitot = 0.0_rp
        ditot = 0.0_rp
        dvene = 0.0_rp
        prene = 0.0_rp
        cvene = 0.0_rp
        hoene = 0.0_rp
        vterm = 0.0_rp
        adiff = 0.0_rp
        cterm = 0.0_rp
        ugrau = 0.0_rp

        if(kfl_ncons_nsa == 0) then
           !
           ! Momentum in conservative formulation
           !
           do idime=1,ndime
              velix = xvelo(idime,igaus,1)
              do jdime=1,ndime
                 ugrau(idime) = ugrau(idime) + & 
                      xvelo(jdime,igaus,1) * gumom(jdime,idime,igaus,1)  
              end do
              cterm(idime) = dvelo * xumom(idime,igaus,1) + ugrau(idime)
           end do
        else
           !
           ! Momentum in non-conservative formulation
           !
           cterm(1) = &
                xvelo(1,igaus,1)*gvelo(1,1,igaus,1) + &
                xvelo(ndime,igaus,1)*gvelo(ndime,1,igaus,1)
           
           cterm(ndime) = &
                xvelo(1,igaus,1)*gvelo(1,ndime,igaus,1) + &
                xvelo(ndime,igaus,1)*gvelo(ndime,ndime,igaus,1)
           
        end if


        !
        ! Viscosity actions: NOTE: this is for PHYSICAL diffusion, not ARTIFICIAL DIFFUSION. 
        !
        if (kexvi > 0) then
           ! viscosity action           
           xdila = cppra_nsa * xvisc(igaus,1)
           xditu = cpprt_nsa * xvitu(igaus,1)
           vitot = xvisc(igaus,1) + xvitu(igaus,1)

           !if(kfl_adiff_nsa > 0) vitot = visco_nsa

           ditot = xdila          + xditu
           vhemo = 0.0_rp

           do idime=1,ndime
              velix = xvelo(idime,igaus,1)

              vehea(idime) = vehea(idime) &
                   - gdens(idime,igaus,1) * ditot / cvcoe_nsa / xdens(igaus,1) / xdens(igaus,1)
              vhemo = vhemo + vehea(idime) * vehea(idime)              
              if (kexvi > 0) then
                 !
                 ! 1. kexvi > 0  -->  viscosity contributions
                 !
                 ! momentum contributions for EXPLICIT scheme
                 vterm(idime) = - vitot*dgvel(idime,igaus,1) - vitot/3.0_rp*gdvel(idime,igaus,1)
                 !adiff(idime) = - vitot*dgvel(idime,igaus,1)

                 do jdime=1,ndime
                    cvene = cvene + gvelo(jdime,idime,igaus,1) * & 
                         (gvelo(jdime,idime,igaus,1) + gvelo(idime,jdime,igaus,1))  
                    vterm(idime) = vterm(idime) - (gvisc(jdime,igaus,1) + & 
                         gvitu(jdime,igaus,1)) * (gvelo(idime,jdime,igaus,1) + gvelo(jdime,idime,igaus,1))
                 end do

                 tehea(idime) = ditot * gtemp(idime,igaus,1)

                 ! both constant and non-constant mu constributions
                 vterm(idime) = vterm(idime) + (gvisc(idime,igaus,1) + gvitu(idime,igaus,1)) * dvelo * dotre

                 ! Heat equation terms:
                 cvene = vitot * ( cvene - dotre * dvelo)
                 cvene = cvene / cvcoe_nsa / xdens(igaus,1)

              end if

              if (kdivi > 0) then
                 !
                 ! 2. kdivi > 0  -->  non constant viscosity contributions
                 !
                 !....       DVENE =  gamma / Pr / rho * grad(visco)_i * grad(T)_i 
                 !
                 !              

                 gdito =  cppra_nsa * gvisc(idime,igaus,1) + cpprt_nsa * gvitu(idime,igaus,1) 
                 dvene = dvene + gdito * gtemp(idime,igaus,1)                  
              end if

           end do

           if (khior > 0) then                         ! high order elements
              hoene = hoene + ditot * xlate(igaus,1)
           end if
           dvene = (dvene + hoene) / cvcoe_nsa / xdens(igaus,1)
           vhemo = sqrt(vhemo)
        end if
        
        cvaux = cvcoe_nsa 
        prene = dvelo * xpres(igaus,1) / cvaux / xdens(igaus,1)
        if (lawst_nsa == 4) then  ! zero when solving potential temperature
           prene = 0.0_rp
        end if

        if ( kfl_isent_nsa == 0) then ! when flow is non-isentropic
           cterm(ndofe) = -vegte      !+ gdgte * ditot / cvcoe_nsa / xdens(igaus,1) / xdens(igaus,1) 
           spare(ndofe) =  cterm(ndofe) !+ cvene - prene + dvene
        end if
        
        !
        ! Compute volume forces        
        !
        xdper = xdens(igaus,1) - xdhyd(igaus)
        vofor = 0.0_rp

        if(kfl_ncons_nsa == 0) then
           !
           ! Momentum in conservative formulation:
           !
           do idime=1,ndime
              vofor(idime) = grnor_nsa * gravi_nsa(idime) * xdper
           end do

        else
           !
           ! Momentum in NON-conservative formulation:
           !
           do idime= 1,ndime
              vofor(idime)= grnor_nsa*gravi_nsa(idime)*xdper/(xdper + xdhyd(igaus))
           end do
           
        end if

        if(kfl_benme_nsa >= 200) then
           !
           ! Kessler microphysics:
           !
           ! Total rain water mixing ratio:
           ! qt_k = 0.608*(qv -qv_hy) + qc + qr
           !
           ! bq(ndime) = bq(ndime) + qt_k*gravity(ndime)
           !
           qt_k = 0.608_rp*(xqvap(igaus,1) - xqvap_hy(igaus,1)) + xqclo(igaus,1) + xqrai(igaus,1)
           
           if(kfl_ncons_nsa == 1) then
              do idime= 1,ndime
                 vofor(idime) = vofor(idime) + grnor_nsa*gravi_nsa(idime)*qt_k
              end do
           else
              do idime= 1,ndime
                 vofor(idime) = vofor(idime) + grnor_nsa*gravi_nsa(idime)*qt_k*xdtot
              end do
           end if
           
        end if

        ! Stabilization parameters
        reate = dvelo
        if (dvelo <= 0.0_rp) then
          ! reate = - dvelo
           reate = 0.0_rp
        end if
        sound = 0.0_rp
        if (kfl_foreg_nsa == 0) sound = sqrt(adgam_nsa * xpres(igaus,1) / xdens(igaus,1))

        taupa = 0.0_rp
        tauen = 0.0_rp
        tauco = 0.0_rp
        hleti = hconv * htran

        seige = sound
        veige = velmo

        if (kfl_lopre_nsa == 1) then
           !
           ! Weiss-Smith local preconditioning
           !
           reate = 0.0_rp
           eplop = 0.00001_rp                    ! as it is given in the paper...
           eplop = eplop * sound
           velop = velmo
           kexvi = kfl_turbu_nsa + kfl_visco_nsa           

           veaux  = velop
          !! if (velop .lt. eplop) veaux = eplop
           if (velop .lt. eplop) veaux = sound
           if (velop .gt. sound) veaux = sound 
           velop = veaux
           if (kexvi > 0) then
              velvi = 2.0_rp * vitot / hleng(ndime) / xdens(igaus,1) 
              if (velvi > velop) then
                 velop = velvi
              end if
           end if
           gamso = adgam_nsa / sound / sound    
           xdenp =   gamso ! derivative of densi w.r.t. press
           xdent = - gamso * xpres(igaus,1) / xtemp(igaus,1) ! derivative of densi w.r.t. tempe
           bepre = ( xdenp + xdent / xdens(igaus,1) / cpcoe_nsa )
           alpre = ( 1.0_rp - bepre * velop * velop) * 0.5_rp

           veige = velmo * (1.0_rp - alpre)
           seige = sqrt (alpre*alpre*velmo*velmo + velop*velop)

        end if        
        !
        ! Dimensionalization factors for VMS
        !   
        !  Momentum eq:     

        dimst(1,1) = 1.0_rp
       !! vedim = seige + veige

        vedim = speed_nsa / rmach_nsa + speed_nsa
        dimst(1,2) = vedim * vedim
        dimst(1,3) = dimst(1,2) * densi_nsa * densi_nsa / tempe_nsa / tempe_nsa
        !  Continuity eq:     
        dimst(2,1) = 1.0_rp / dimst(1,2)
        dimst(2,2) = 1.0_rp
        dimst(2,3) = dimst(1,3) / dimst(1,2)
        !  Heat eq:     
        dimst(3,1) = 1.0_rp / dimst(1,3)
        dimst(3,2) = 1.0_rp / dimst(2,3)
        dimst(3,3) = 1.0_rp

        tauxi = 0.0_rp
        xdtau = xdens(igaus,1)
        if (kfl_stabi_nsa >= 1) then
           if (kfl_taufa_nsa(1,2) == 1) then
              tauxi(1,1) = veige/hconv
              tauxi(2,1) = veige/hconv
           end if
           if (kfl_taufa_nsa(2,2) == 1) then
              tauxi(1,2) = reate
           end if
           if (kfl_taufa_nsa(3,2) == 1) then
              tauxi(1,3) = seige / hleng(ndime)
!!!!!!!PROVAR amb aquest terme també!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ! tauxi(2,3) = seige / hleng(ndime)
           end if
           if (kfl_taufa_nsa(4,2) == 1) then
              tauxi(1,4) = 4.0_rp * vitot / xdtau / hleti
              tauxi(2,4) = 4.0_rp * ditot / xdtau / cvcoe_nsa / hleti
           end if

           !           tauxi(1,1) = tauxi(1,1) * tauxi(1,1)  
           !           tauxi(1,2) = tauxi(1,2) * tauxi(1,2)   
           !           tauxi(1,3) = tauxi(1,3) * tauxi(1,3)  
           !           tauxi(1,4) = tauxi(1,4) * tauxi(1,4)  
           !           tauxi(2,1) = tauxi(2,1) * tauxi(2,1)  
           !           tauxi(2,2) = tauxi(2,2) * tauxi(2,2)   
           !           tauxi(2,3) = tauxi(2,3) * tauxi(2,3)  
           !           tauxi(2,4) = tauxi(2,4) * tauxi(2,4)  


           tiint = tauxi(1,1) + tauxi(1,2) + tauxi(1,3) + tauxi(1,4)
           tiine = tauxi(2,1) + tauxi(2,2) + tauxi(2,3) + tauxi(2,4)
           tiinc = tiint

           !           tiine= sqrt(tiine)
           !           tiint= sqrt(tiint)
           !           tiinc= sqrt(tiinc)

           if (tiint > zensa) taupa = 0.5_rp * qufac / tiint
           if (tiinc > zensa) tauco = 0.5_rp * qufac / tiinc
           if (tiine > zensa) tauen = 0.5_rp * qufac / tiine
          !! if (tiine > zensa) tauen = taupa

           if(kfl_tisch_nsa == 41) then
              tauen = parkb_nsa(itinn(modul),1) * tauen
              taupa = parkb_nsa(itinn(modul),1) * taupa
              tauco = parkb_nsa(itinn(modul),1) * tauco              
           end if
           

           if (kfl_stabi_nsa == 2) then ! local tau for pressure
              tauco = tauco * safet_nsa
           end if

        end if

        dundt = 0.0_rp

        elemental_nodes: do inode=1,pnode

           ipoin = lnods(inode,ielem)
           xshap = elmar(pelty)%shape(inode,igaus)

           tiext = dtieq_nsa(2,ipoin,1)

           dundt(    1  ) = (xumom(    1,igaus,3) - xumom(    1,igaus,2)) / tiext
           dundt(    2  ) = (xumom(    2,igaus,3) - xumom(    2,igaus,2)) / tiext
           dundt(ndime  ) = (xumom(ndime,igaus,3) - xumom(ndime,igaus,2)) / tiext
           dundt(ndime+1) = (xdens(      igaus,3) - xdens(      igaus,2)) / tiext
           dundt(ndime+2) = (xtemp(      igaus,3) - xtemp(      igaus,2)) / tiext

           if (kfl_stabi_nsa == 1) then
              tauco = tiext   ! use exterior time step, i.e., global pressure tau
           end if

           ! div(N u) = u grad(N)     +     N div(u)
           dshve = xvcar(inode,igaus,1) + xshap * dvelo

           ! Galerkin terms:
           galte = 0.0_rp
           ! Stabilization terms:
           state = 0.0_rp
           ! Shock capturing terms:
           shote = 0.0_rp
           resdc = 0.0_rp
           resta = 0.0_rp

           !
           ! Galerkin and residual of continuity equation:
           ! Atention!!!!!! resta(ndofc) está cambiado de signo!!!!!!!!!!! tener en cuenta!!!!!!!!


           ! For pseudo-artificial diffusion with constant tau:
           !
           if(kfl_adiff_nsa > 0 ) then
              tauco = kdiff_nsa
              taupa = kdiff_nsa
              tauen = kdiff_nsa
           end if

           !
           ! Galrkin and residual of continuity equation:
           !
           galte(ndofc)= - xshap * dumom

           resta(ndofc) = &
                + theta_nsa(1) * dfmom + (1.0_rp - theta_nsa(1)) * dumom  &
                - theta_nsa(1) * tauco * xlapr(igaus,1) 
           resdc(ndofc) =  dundt(ndofc) + resta(ndofc)
           if (kfl_shock_nsa(2) == 11) resdc(ndofc) =  resta(ndofc)
           if ( kfl_isent_nsa == 0) then ! when flow is non-isentropic
              resta(ndofe) = spare(ndofe)
              resdc(ndofe) = - dundt(ndofe) + resta(ndofe)
              if (kfl_shock_nsa(3) == 11) resdc(ndofe) = resta(ndofe)
              galte(ndofe) = xshap * spare(ndofe)
           end if
 
           !
           ! Residual momentum equations
           !
           do idime=1,ndime          
              resta(idime) = - cterm(idime) - (gpres(idime,igaus,1)-gphyd(idime,igaus)) + vofor(idime) !- vterm(idime) - adiff(idime)
           end do

           !
           ! Stabilization terms
           !
           if (kfl_stafl_nsa == 0) then
              if (kfl_stabi_nsa >= 5) then
                 if (lawst_nsa == 4) then                 
                    resta(ndime+1) = 0.0_rp
                    resta(ndime+2) = 0.0_rp
                    do jdime=1,ndime
                       resta(ndime+1) = resta(ndime+1) - gumom(jdime,jdime,igaus,1)
                       resta(ndime+2) = resta(ndime+2) - xvelo(jdime,igaus,1) * gtemp(jdime,igaus,1)
                    end do
                    do kdime=1,ndime
                       state(ndime+1) = state(ndime+1) + ( &
                            cartd(kdime,inode,igaus) * xconv(ndime+1,ndime+1,kdime,igaus) &
                            + dconv(ndime+1,ndime+1,igaus) * xshap &
                            ) * tauco * resta(ndime+1)
                       state(ndime+1) = state(ndime+1) + ( &
                            cartd(kdime,inode,igaus) * xconv(ndime+1,ndime+2,kdime,igaus) &
                            + dconv(ndime+1,ndime+2,igaus) * xshap &
                            ) * tauen * resta(ndime+2)
                       
                       state(ndime+2) = state(ndime+2) + ( &
                            cartd(kdime,inode,igaus) * xconv(ndime+2,ndime+1,kdime,igaus) &
                            + dconv(ndime+2,ndime+1,igaus) * xshap &
                            ) * tauco * resta(ndime+1)
                       state(ndime+2) = state(ndime+2) + ( &
                            cartd(kdime,inode,igaus) * xconv(ndime+2,ndime+2,kdime,igaus) &
                            + dconv(ndime+2,ndime+2,igaus) * xshap &
                            ) * tauen * resta(ndime+2)
                       do jdime=1,ndime                    
                          state(ndime+1) = state(ndime+1) + ( &
                               cartd(kdime,inode,igaus) * xconv(ndime+1,jdime,kdime,igaus) &
                               + dconv(ndime+1,jdime,igaus) * xshap &
                               ) * taupa * resta(jdime)
                          
                          state(ndime+2) = state(ndime+2) + ( &
                               cartd(kdime,inode,igaus) * xconv(ndime+2,jdime,kdime,igaus) &
                               + dconv(ndime+2,jdime,igaus) * xshap &
                               ) * taupa * resta(jdime)
                       end do
                    end do
                 end if
              end if
           else
              
              if (kfl_stabi_nsa < 5) then
                 state(ndofe) = tauen * dshve * spare(ndofe)
                ! state(ndofc) = xcart * (thet1 * xmofr - thpre * tauco * gpper)                 
                 state(ndofc) = 0.0_rp
              else
                 state(ndofe) = (xvcar(inode,igaus,1) + dvelo * xshap) * tauen * resta(ndofe)
                ! state(ndofe) =  xvcar(inode,igaus,1) * tauco * resta(ndofe)
              end if
           end if
          
           gvcar = 0.0_rp
           gtcar = 0.0_rp
           momentum_dimensions: do idime=1,ndime
              !gvcar = 0.0_rp
              xcart = cartd(idime,inode,igaus)
              xmofr = xfmom(idime,igaus,1) - xumom(idime,igaus,1)

              !
              ! dP'/dxi
              !
              gpper = (gpres(idime,igaus,1) - gphyd(idime,igaus))


!!! --->>> check this for runge-kuttas.. 
             ! cogra = thet2*gpres(idime,igaus,1) + theco*gpold(idime,igaus,1)

              cogra = - theco * gpper
              vweak = 0.0_rp
              tweak = 0.0_rp

              if ( kfl_isent_nsa == 0) then ! when flow is non-isentropic
                 galte(ndofe) = galte(ndofe) - xcart * tehea(idime)
              end if

              fpert = 0.0_rp

              ! momentum equation     
              if (kfl_stabi_nsa >= 1) then
                 fpert = taupa * dshve
              end if

              if(kfl_adiff_nsa > 0) then
                 gvcar = 0.0_rp
                 gtcar = 0.0_rp
                 do jdime=1,ndime
                    !Stress tensor
                    gvcar = gvcar + cartd(jdime,inode,igaus) * gvelo(jdime,idime,igaus,1)
                    gtcar = gtcar + cartd(jdime,inode,igaus) * gtemp(jdime,igaus,1)
                 end do
                 
                 vweak = vitot * gvcar
                 tweak = vitot * gtcar
              end if
              
              !gvcar = 0.0_rp
              !if (kexvi > 1) then
              !   do jdime=1,ndime
              !      !Stress tensor
              !      gvcar = gvcar + cartd(jdime,inode,igaus) * &
              !           (gvelo(jdime,idime,igaus,1) + gvelo(idime,jdime,igaus,1))
              !   end do
              !   vweak = vitot * gvcar &
              !        - (cartd(idime,inode,igaus) * (dotre * vitot * dvelo + xkene(igaus,1)) * 0.66666666666666666667_rp)
              !
              !
              !  
              !       
              !end if


              !
              ! Build momentum rhs:
              !
              if(kfl_ncons_nsa == 0) then
                 !
                 ! Momentum in conservative formulation:
                 !
                 cogra = - cterm(idime) - gpper + vofor(idime)
                 cogri = - gpper
                 
              else
                 !
                 ! Momentum in NON-conservative formulation:
                 !
                 cogra = - cterm(idime) - gpper/xdtot + vofor(idime)
                 cogri =                - gpper/xdtot
                 
              end if
              
              resdc(idime) = - dundt(idime) + resta(idime)
              if (kfl_shock_nsa(1) == 11) resdc(idime) = resta(idime)
              !
              ! Shock capturing terms:              
              !
              shote = 0.0_rp
              if (kfl_shock_nsa(1) > 0 ) then
                 ifcdr  = 1
                 vidis  = vitot
                 gshoc(    1)= gumom(    1,idime,igaus,1)
                 gshoc(    2)= gumom(    2,idime,igaus,1)
                 gshoc(ndime)= gumom(ndime,idime,igaus,1)
                 call nsa_shocap(ifcdr, &
                      kfl_shock_nsa(1),ndime,shock_nsa,gshoc,cartd(1,inode,igaus), &
                      resdc(idime),hconv,xvelo(1,igaus,1), &
                      taupa,vidis,1000._rp*zensa,qufac,shote(idime))
              end if !end shock capturing
                         
              !
              ! Galerkin term of momentum
              !
              galte(idime) = - vweak + xshap * cogra
              if(kfl_adiff_nsa > 0) &
                   galte(ndofe) = galte(ndofe) - tweak
              
              !
              ! Stabilization terms
              !
              if (kfl_stafl_nsa == 0) then                                 
                 if (kfl_stabi_nsa >= 5) then
                    if (lawst_nsa == 4) then 
                       state(idime) = state(idime) + grnor_nsa * real(gravi_nsa(idime)) * xshap * tauco * resta(ndime+1)
                       do kdime=1,ndime
                          do jdime=1,ndime
                            ! resta(jdime) =  - cterm(jdime) - (gpres(jdime,igaus,1)-gphyd(jdime,igaus)) + vofor(jdime)
                             state(idime) = state(idime) + ( &
                                  cartd(kdime,inode,igaus) * xconv(idime,jdime,kdime,igaus) &
                                 ! + dconv(idime,jdime,igaus) * xshap &
                                  ) * taupa * resta(jdime)
                          end do
                          state(idime) = state(idime) + ( &
                               cartd(kdime,inode,igaus) * xconv(idime,ndime+1,kdime,igaus) &
                              ! + dconv(idime,ndime+1,igaus) * xshap &
                               ) * tauco * resta(ndime+1)
                          state(idime) = state(idime) + ( &
                               cartd(kdime,inode,igaus) * xconv(idime,ndime+2,kdime,igaus) &
                              ! + dconv(idime,ndime+2,igaus) * xshap &
                               ) * taupa * resta(ndime+2)
                       end do
                    end if
                 end if
              else
                 if (kfl_stabi_nsa < 5) then
                    state(idime) = fpert * resta(idime)
                 else
                    state(idime) = xvcar(inode,igaus,1) * taupa * (resta(idime)) &
                         - xshap * grnor_nsa * gravi_nsa(idime) * taupa * resta(ndofc)
                    if (lawst_nsa == 1) then ! ideal gases
                       state(ndofe) =  state(ndofe) &
                            + xcart * xdhyd(igaus) * rgasc_nsa &
                            * taupa * (resta(idime))
                    else if (lawst_nsa == 4) then ! potential temperature equation
                       state(idime) = state(idime) &
                            + xcart * afact_nsa * xdens(igaus,1) * (xdens(igaus,1)*xtemp(igaus,1))**(adgam_nsa - 1.0_rp) &
                            * taupa * resta(ndofe) * 0.5_rp &  ! contribution for a possible linearization
                            - xcart * afact_nsa * xtemp(igaus,1) * (xdens(igaus,1)*xtemp(igaus,1))**(adgam_nsa - 1.0_rp) &
                            * tauco * resta(ndofc) * 0.5_rp ! contribution for a possible linearization
                          
                           ! * tauco * resta(ndofc) * adgam_nsa ! contribution for a possible linearization
                           ! * taupa * resta(ndofe) * adgam_nsa &  ! contribution for a possible linearization
                    end if
                    ! continuity equation:
                    state(ndofc) = state(ndofc) + cartd(idime,inode,igaus) * taupa * resta(idime)
                 end if
                 
              end if

              !Mapping the index from tensor to vector:
              ievav = (inode-1) * ndofn_nsa + idime

              !Build and integrate RHS of momentum to global system:
              elrhs(ievav) = elrhs(ievav)  + dvolu(igaus) * (state(idime) + shote(idime) + galte(idime))
              
              !
              ! Nodal stabilization term
              !
              elsta(ievav) = elsta(ievav) + dvolu(igaus) * (state(idime) )! + shote(idime))
              
              if (kfl_diagi_nsa > 0) then
                 
                ! dcart= cartd(1,inode,igaus)
                ! if (ndime == 2) then
                !    dcart=dcart+cartd(2,inode,igaus)
                ! else if (ndime == 3) then
                !    dcart=dcart+cartd(2,inode,igaus)+cartd(3,inode,igaus)
                ! end if
                 
                 galim(idime) =     xshap * (dshve + cogri)
                 staim(idime) = 0.0_rp
                ! staim(idime) =     taupa * 2.0_rp * xshap * xvcar(inode,igaus,1) * resta(idime)/elden(inode,1)
                 eldia(ievav) = eldia(ievav) + dvolu(igaus) * (galim(idime)+staim(idime))
              end if
           end do momentum_dimensions

           ! Shock capturing terms:              
           if (kfl_shock_nsa(2) > 0) then
              ! Continuity equation
              vidis        = 0.0_rp
              gshoc(    1) = gdens(    1,igaus,1)
              gshoc(    2) = gdens(    2,igaus,1)
              gshoc(ndime) = gdens(ndime,igaus,1)

              ifcdr        = 0
              call nsa_shocap(ifcdr,                                                             &
                   kfl_shock_nsa(2),ndime,shock_nsa,gshoc,cartd(1,inode,igaus),                  &
                   resdc(ndofc),hconv,xvelo(1,igaus,1),                                          &
                   tiext,vidis,1000._rp*zensa,qufac,shote(ndofc))
           end if

           if (kfl_shock_nsa(3) > 0) then

              if ( kfl_isent_nsa == 0) then                  ! when flow is non-isentropic

                 ! Heat equation

                 vidis        = ditot / cvcoe_nsa
                 gshoc(    1) = gtemp(    1,igaus,1)
                 gshoc(    2) = gtemp(    2,igaus,1)
                 gshoc(ndime) = gtemp(ndime,igaus,1)
                 ifcdr        = 1
                 call nsa_shocap(ifcdr,                                                             &
                      kfl_shock_nsa(3),ndime,shock_nsa,gshoc,cartd(1,inode,igaus),                  &
                      resdc(ndofe),hconv,vehea,                                                     &
                      tauen,vidis,1000._rp*zensa,qufac,shote(ndofe))
              end if

           end if

           !Build and integrate right hand side of continuity equation
           ievac = (inode-1) * ndofn_nsa + ndofc 
           elrhs(ievac)= elrhs(ievac)  + dvolu(igaus) * (galte(ndofc) + state(ndofc) + shote(ndofc))
           elsta(ievac)= elsta(ievac)  + dvolu(igaus) * (state(ndofc) )!+ shote(ndofc))

           ievae = (inode-1) * ndofn_nsa + ndofe 
           elrhs(ievae)= elrhs(ievae)  + dvolu(igaus) * (galte(ndofe) + state(ndofe) + shote(ndofe))
           elsta(ievae)= elsta(ievae)  + dvolu(igaus) * (state(ndofe) )!+ shote(ndofe))

           if (kfl_diagi_nsa > 0) then
              galim(ndofc) =   xshap * dshve
              galim(ndofe) =   xshap * dshve
              eldia(ievae)= eldia(ievae) + dvolu(igaus)*galim(ndofe)
              eldia(ievac)= eldia(ievac) + dvolu(igaus)*galim(ndofc)
           end if

           !
           ! Global SGS on ipoin:
           !
           ipoin= lnods(inode,ielem)
           xshai= elmar(pelty)%shape(inode,igaus)
           
           asfac= dvolu(igaus) * xshai / vmasc(ipoin)
           
           umoss_nsa(    1,ipoin,2) = umoss_nsa(    1,ipoin,2) + asfac * taupa*resta(1)
           umoss_nsa(    2,ipoin,2) = umoss_nsa(    2,ipoin,2) + asfac * taupa*resta(2)
           if (ndime == 3) &
                umoss_nsa(ndime,ipoin,2) = umoss_nsa(ndime,ipoin,2) + asfac * taupa*resta(ndime)
           
           denss_nsa(ipoin,2) = denss_nsa(ipoin,2) - asfac * tauco * resta(ndime+1)
           eness_nsa(ipoin,2) = eness_nsa(ipoin,2) + asfac * tauen * resta(ndime+2)

        end do elemental_nodes
           
     end do elemental_gauss_points

     elemental_assembly_nodes: do inode= 1,pnode
        ipoin= lnods(inode,ielem)
        do idofn= 1,ndofn_nsa
           ievat = (inode-1)*ndofn_nsa + idofn
           jevat = (ipoin-1)*ndofn_nsa + idofn

           !$*OMP      CRITICAL (lock_rhsid_vdiag)
           rhsid(jevat) = rhsid(jevat) + elrhs(ievat)
           
           !Assembly of Stabilization term
           !stabi_nsa(jevat) = stabi_nsa(jevat) + elsta(ievat)

           !           if (kfl_diagi_nsa > 0) then
           !              vdiag_nsa(idofn,ipoin)= vdiag_nsa(idofn,ipoin) + eldia(ievat)
           !           end if
           !$*OMP      END CRITICAL (lock_rhsid_vdiag)

        end do
     end do elemental_assembly_nodes

     !    For non-fractional schemes: 
     !    Compute boundary contributions when needed
     !    by looping over nfabo, as stored in leobl_nsa. if nfabo=0 (i.e. ielem has no
     !    faces which are boundary elements) the loop is skipped.
     !    Also compute ccoefs, forces and pitches. 
     !     

!!$     if (kfl_genal_nsa==0) then
!!$        nfabo = 0
!!$        if (nboun > 0) nfabo = leobl_nsa(ielem)%l(pface+1)
!!$        
!!$        ielem_faces: do ifabo= 1, nfabo
!!$           iboun = leobl_nsa(ielem)%l(ifabo)
!!$           pblty = ltypb(iboun) 
!!$           pnodb = nnode(pblty)
!!$           pgaub = ngaus(pblty)
!!$           
!!$           do inodb= 1,pnodb
!!$              inode = lboel(inodb,iboun)
!!$              bocod(1:ndime,inodb  ) = elcod(1:ndime,inode  )
!!$           end do
!!$           
!!$           iboun_gauss_points: do igaub=1,pgaub
!!$              
!!$              xvelo(1:ndime,1,1) =0.0_rp
!!$              xdens(        1,1) =0.0_rp
!!$              xpres(        1,1) =0.0_rp
!!$              xtemp(        1,1) =0.0_rp
!!$              xvisc(        1,1) =0.0_rp
!!$              xvitu(        1,1) =0.0_rp
!!$              do idime=1,ndime
!!$                 gvelo(1:ndime,idime,1,1)=0.0_rp
!!$                 gpres(        idime,1,1)=0.0_rp
!!$                 gdens(        idime,1,1)=0.0_rp
!!$                 gtemp(        idime,1,1)=0.0_rp
!!$                 gvisc(        idime,1,1)=0.0_rp
!!$                 gvitu(        idime,1,1)=0.0_rp
!!$              end do
!!$              
!!$              call bouder(pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),bocod,baloc,eucta)    
!!$              call chenor(pnode,baloc,bocod,elcod)    
!!$              dsurf =elmar(pblty)%weigp(igaub)*eucta 
!!$              
!!$              shapp= 0.0_rp
!!$              do inodb=1,pnodb
!!$                 shapp(lboel(inodb,iboun)) =  elmar(pblty)%shape(inodb,igaub)
!!$              end do
!!$
!!$              do inode=1,pnode
!!$                 cartb(1:ndime,inode) = 0.0_rp
!!$                 do inodb=1,pnodb
!!$                    do igaus=1,pgaus
!!$                       do idime=1,ndime
!!$                          cartb(idime,inode)=cartb(idime,inode)                   &
!!$                               +elmar(pelty)%shaga(igaus,lboel(inodb,iboun))          &
!!$                               *elmar(pblty)%shape(inodb,igaub)*cartd(idime,inode,igaus)
!!$                       end do
!!$                    end do                    
!!$                 end do
!!$              end do
!!$              
!!$              do inode=1,pnode
!!$                 xshap = shapp(inode)
!!$                 xdens(1,1) = xdens(1,1) + xshap * elden(inode,1)
!!$                 xpres(1,1) = xpres(1,1) + xshap * elpre(inode,1)
!!$                 xtemp(1,1) = xtemp(1,1) + xshap * eltem(inode,1)
!!$                 xvisc(1,1) = xvisc(1,1) + xshap * elvis(inode,1)
!!$                 xvitu(1,1) = xvitu(1,1) + xshap * elvit(inode,1)
!!$                 do idime=1,ndime
!!$                    xcart = cartb(idime,inode)
!!$                    xvelo(idime,1,1) = xvelo(idime,1,1) + elvel(idime,inode,1) * xshap
!!$                    gdens(idime,1,1) = gdens(idime,1,1) + elden(inode,1)       * xcart
!!$                    gtemp(idime,1,1) = gtemp(idime,1,1) + eltem(inode,1)       * xcart
!!$                    gpres(idime,1,1) = gpres(idime,1,1) + elpre(inode,1)       * xcart                 
!!$                    gvisc(idime,1,1) = gvisc(idime,1,1) + elvis(inode,1)       * xcart                 
!!$                    gvitu(idime,1,1) = gvitu(idime,1,1) + elvit(inode,1)       * xcart                 
!!$                    do jdime= 1,ndime
!!$                       gvelo(idime,jdime,1,1) = gvelo(idime,jdime,1,1) + elvel(jdime,inode,1) * xcart
!!$                    end do
!!$                 end do
!!$              end do
!!$              
!!$              dvelo = 0.0_rp
!!$              velmo = 0.0_rp
!!$              if (lawvi_nsa > 0) call nsa_lawvis( -1 , 1 ,xvisc(1,1),xtemp(1,1),dvite)
!!$              do idime=1,ndime
!!$                 dvelo = dvelo + gvelo(idime,idime,1,1)
!!$                 velmo = velmo + xvelo(idime,1,1) * xvelo(idime,1,1)                     
!!$                 gvisc(idime,1,1) = gvisc(idime,1,1) * dvite
!!$              end do
!!$              velmo=sqrt(velmo)
!!$              vitot= xvisc(1,1)+xvitu(1,1)
!!$              
!!$              ! Compute surface aerodynamic body forces              
!!$              call nsa_aebody(1,iboun,ielem,pnodb,ifibo)
!!$              
!!$              if (ifibo>0) then
!!$                 
!!$                 do idime= 1,ndime
!!$                    tensi(idime)=0.0_rp
!!$                    tentu(idime)=0.0_rp
!!$                    
!!$                    do jdime= 1,ndime
!!$                       tensi(idime)= tensi(idime) & 
!!$                            + baloc(jdime,ndime) * (gvelo(idime,jdime,1,1)+gvelo(jdime,idime,1,1))
!!$                    end do
!!$                    
!!$                    tensi(idime)=  xvisc(1,1) * tensi(idime) &
!!$                         - baloc(idime,ndime) * xvisc(1,1) * dvelo *  dotre             
!!$                 end do
!!$                 
!!$                 surto_nsa = surto_nsa + dsurf
!!$                 banor(1:ndime) = baloc(1:ndime,ndime)
!!$                 trapr(1:ndime) = xpres(1,1) * banor(1:ndime) 
!!$#ifdef _OPENMP
!!$                 do idime=1,ndime
!!$                    if (banor(idime) < 0.0_rp) & 
!!$                         sufro_nsaOMP(idime,thrid) = sufro_nsaOMP(idime,thrid) + dsurf * banor(idime)
!!$                    fbody_nsaOMP(1,idime,thrid)=fbody_nsaOMP(1,idime,thrid) + dsurf*trapr(idime)   ! pressure
!!$                    fbody_nsaOMP(2,idime,thrid)=fbody_nsaOMP(2,idime,thrid) + dsurf*tensi(idime)   ! viscous
!!$                 end do
!!$#else
!!$                 do idime=1,ndime
!!$                    if (banor(idime) < 0.0_rp) & 
!!$                         sufro_nsa(idime) = sufro_nsa(idime) + dsurf * banor(idime)
!!$                    fbody_nsa(1,idime)=fbody_nsa(1,idime) + dsurf*trapr(idime)   ! pressure
!!$                    fbody_nsa(2,idime)=fbody_nsa(2,idime) + dsurf*tensi(idime)   ! viscous
!!$                 end do
!!$#endif
!!$              end if
!!$              
!!$           end do iboun_gauss_points
!!$           
!!$        end do ielem_faces
!!$     
!!$     end if

  end do elements_loop

  !
  ! Write the subgrid scales to file:
  !
!!$   !
!!$  !Write vtk every n-steps
!!$  !
!!$  if(ittim == 1 .or. mod(ittim,irestart) == 0 .or. ittim == mitim) then
!!$     istep_nsa = istep_nsa + 1
!!$     write(fnp1,'(i3)')istep_nsa
!!$     iloop=2 - int(log10(real(istep_nsa)))
!!$     do j=1,iloop
!!$        fnp1(j:j)='0'
!!$     end do
!!$     !VTK:
!!$     fnp=trim('VTK_SGSmome') // '_' // trim(fnp1) // '.vtk'
!!$     print*, fnp
!!$     call nsa_outvtk_subscales(fnp)
!!$
!!$  end if

  !
  ! Distribute global subscale fields in parallel runs
  !
  call nsa_parall(7_ip) 
  do ipoin = 1,npoin
     umoss_nsa(    1,ipoin,1) = umoss_nsa(    1,ipoin,2)
     umoss_nsa(    2,ipoin,1) = umoss_nsa(    2,ipoin,2)
     if (ndime == 3) umoss_nsa(ndime,ipoin,1) = umoss_nsa(ndime,ipoin,2)
     denss_nsa(      ipoin,1) = denss_nsa(      ipoin,2)
     eness_nsa(      ipoin,1) = eness_nsa(      ipoin,2)
  end do


#ifdef _OPENMP
  numthr   = OMP_GET_MAX_THREADS()
  do thrid = 1, numthr
     do idime = 1, ndime
        sufro_nsa(idime) = sufro_nsa(idime) + sufro_nsaOMP(idime,thrid)
        fbody_nsa(1,idime) = fbody_nsa(1,idime) + fbody_nsaOMP(1,idime,thrid)
        fbody_nsa(2,idime) = fbody_nsa(2,idime) + fbody_nsaOMP(2,idime,thrid)
     enddo
  enddo
#endif
  !
  ! Aerodynamic body forces when needed (and was not done before)
  !
  call nsa_aebody(2,1,1,1,1)

end subroutine nsa_elchea
