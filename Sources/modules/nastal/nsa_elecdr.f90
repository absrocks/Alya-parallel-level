!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_elecdr.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   CDR-in-nastal 
!> @details CDR-in-nastal 
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_elecdr
  use      def_master
  use      def_domain
  use      def_nastal
  implicit none

  real(rp)    :: elrhs(nevat_nsa),diago(nevat_nsa),elmat(nevat_nsa,nevat_nsa)
  real(rp)    :: &
       elunk(ndofn_nsa,mnode,ncomp_nsa),elsub(ndofn_nsa,mnode), elort(ndofn_nsa,mnode), &
       elcod(ndime,mnode),elvel(ndime,mnode), &
       elpre(mnode),eltem(mnode), &
       elcon(ndofn_nsa,ndofn_nsa,ndime,mnode), &
       eldif(ndofn_nsa,ndofn_nsa,ndime,ndime,mnode),conme(ndofn_nsa,ndofn_nsa,ndime), &
       difme(ndofn_nsa,ndofn_nsa,ndime,ndime),eldtt(ndofn_nsa,mnode,2),kapsh(ndofn_nsa),elvis(mnode),elthe(mnode)
  real(rp)    :: &
       xconv(ndofn_nsa,ndofn_nsa,ndime,mgaus),gunkn(ndofn_nsa,ndime,mgaus),grasc(ndofn_nsa,ndime), &
       gsube(ndofn_nsa,ndime,mgaus), taudi(ndofn_nsa), &
       dconv(ndofn_nsa,ndofn_nsa,mgaus),xsube(ndofn_nsa,mgaus,3),xortp(ndofn_nsa,mgaus),xdtix(ndofn_nsa,mgaus,2), &
       xdiff(ndofn_nsa,ndofn_nsa,ndime,ndime,mgaus),ddiff(ndofn_nsa,ndofn_nsa,2,ndime,mgaus),xtime(ndofn_nsa,mgaus)
  real(rp)    :: &
       detjm,qufac,xshai,asfac,dtaux, &
       dvolu(mgaus),hessi(ntens,mnode),xresi(ndofn_nsa,mgaus), &
       cartd(ndime,mnode,mgaus),hesma(ndime,ndime,mnode,mgaus),tragl(ndime,ndime),hleng(ndime), &
       xjaci(ndime,ndime),xjacm(ndime,ndime), &
       xunkn(ndofn_nsa,mgaus,3), &
       dvelo(mgaus), &
       sound(mgaus),xpres(mgaus),xtemp(mgaus),xvisc(mgaus),xdith(mgaus),xlade(mgaus),xldve(ndime,mgaus), &
       xvelo(ndime,mgaus),gpres(ndime,mgaus),gtemp(ndime,mgaus),gvisc(ndime,mgaus),difeq(ndofn_nsa), &
       gvelo(ndime,ndime,mgaus),velmo(mgaus),d2sdx(ndime,ndime,ndime),chale(2), &
       shmet(ndime,ndime,ndofn_nsa),shote(ndofn_nsa), &
       xtide(ndofn_nsa,mgaus) ! Preconditioner * Unknown(last_iteration)-Unknown(last_time_step) / physical_time_step
       
  integer(ip) :: ielem,inode,jdofn,idofn,itott,idime,jdime,igaus,pelty,pnode,pgaus,plapl,pface,ievat,jevat, &
       ipoin,kpoin,ndaux,mfreq


  unkno = 0.0_rp
  if (kfl_stabi_nsa >= 1) call nsa_ortpro_cdr_consistent
 ! if (kfl_stabi_nsa >= 1) call nsa_ortpro_cdr
  rhsid = 0.0_rp
  unkno = 0.0_rp
  amatr = 0.0_rp

  ccoef_nsa = 0.0_rp
  fbody_nsa = 0.0_rp
  surto_nsa     = 0.0_rp
  sufro_nsa     = 0.0_rp
  mfreq = 0

  do ipoin = 1,npoin
     umoss_nsa(    1,ipoin,2) = 0.0_rp
     umoss_nsa(    2,ipoin,2) = 0.0_rp
     umoss_nsa(ndime,ipoin,2) = 0.0_rp
     denss_nsa(      ipoin,2) = 0.0_rp
     eness_nsa(      ipoin,2) = 0.0_rp
     frequ_nsa(ipoin) = 0.0_rp

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

  ndaux= ndofn_nsa
  if (kfl_isent_nsa == 1) ndaux= ndofn_nsa - 1

  elements_loop: do ielem=1,nelem

     ! Element properties and dimensions
     pelty=ltype(ielem)
     pnode=nnode(pelty)
     pgaus=ngaus(pelty)
     plapl=llapl(pelty)
     pface=nface(pelty)
     
     do idime=1,ndime
        do idofn= 1,ndofn_nsa
           do jdofn= 1,ndofn_nsa
              conme(idofn,jdofn,idime) = 0.0_rp
              do jdime= 1,ndime
                 difme(idofn,jdofn,idime,jdime) = 0.0_rp
              end do
           end do
        end do
     end do
     

     qufac = 1.0_rp
     if((ndime.eq.2).and.(pnode.ge.4)) then
        qufac = 0.5_rp 
     end if
     if((ndime.eq.3).and.(pnode.ge.5))then
        qufac = 0.5_rp 
     end if

     ! 1. gather

     call nsa_gaconsxy(&
          ielem,pnode,elcod,elunk,elsub,elcon,eldif,elvel,elpre,eltem,elthe,elvis,eldtt,elort)
!!$do inode=1,pnode
!!$print*,'elden',ielem,inode,elunk(ndime+1,inode,1)
!!$end do

     ! hleng and tragl at center of gravity
     call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)
     ! default values for chale

!!!     hleng(1) is the LARGEST
!!!     hleng(ndime) is the SMALLEST

!!!     chale(1) = hleng(1)          ! largest
!!!     chale(2) = hleng(ndime)      ! smallest

!!!! this looks like to be the best option for high aspect ratio elements. 
!!!! hconv is ALWAYS chale(1), so change chale to some testing
     chale(1) = hleng(ndime)      ! smallest
     chale(2) = hleng(1)          ! largest

     ! compute chale: stream-wise and cross-wise lengths
!     call velchl(pnode,elcod,elvel,chale,hleng) 



     elemental_gauss_points_residuals: do igaus=1,pgaus

        call elmder(pnode,ndime,elmar(pelty)%deriv(1,1,igaus),elcod,cartd(1,1,igaus),detjm,xjacm,xjaci)
        dvolu(igaus)=elmar(pelty)%weigp(igaus)*detjm                
        hessi(1:ntens,1:mnode) = 0.0_rp
        if(plapl==1) then
           call elmhes(elmar(pelty)%heslo(1,1,igaus),hessi,ndime,pnode,ntens,&
                xjaci,d2sdx,elmar(pelty)%deriv(1,1,igaus),elcod)     
        end if
        
        ! 2. calculo de todo en los gauss, incluido adjunto, salvo monyos

  
        call nsa_gauvalxy(ielem,igaus,pnode,pgaus,elmar(pelty)%weigp(igaus), &
             elcon,eldif,elunk,elsub,elort,eldtt,xunkn,xdtix,elpre, &
             xconv(1,1,1,igaus),dconv(1,1,igaus),xdiff(1,1,1,1,igaus), &
             ddiff(1,1,1,1,igaus),gunkn(1,1,igaus),gsube(1,1,igaus), &
             elmar(pelty)%shape(1,igaus),cartd(1,1,igaus),hesma(1,1,1,igaus),hessi,xsube,xortp(1,igaus),xresi(1,igaus), &
             sound(igaus),xpres(igaus),xtemp(igaus),xvisc(igaus),xdith(igaus),xvelo(1,igaus), &
             gpres(1,igaus),gtemp(1,igaus), &
             gvisc(1,igaus),gvelo(1,1,igaus),dvelo(igaus),velmo(igaus),xlade(igaus),xldve(1,igaus), &
             xjaci,conme,difme,dvolu(igaus),xtide(1,igaus))

        if (kfl_isent_nsa==1) xresi(ndime+2,igaus) = 0.0_rp




     end do elemental_gauss_points_residuals

     do ievat=1,nevat_nsa
        elrhs(ievat)= 0.0_rp
        diago(ievat)= 0.0_rp
        do jevat=1,nevat_nsa
           elmat(ievat,jevat)= 0.0_rp
        end do
     end do



     elemental_gauss_points_monyos_scatter: do igaus=1,pgaus

        do idime=1,ndime
           do idofn=1,ndofn_nsa
              grasc(idofn,idime) = gunkn(idofn,idime,igaus)
           end do
        end do
        
        if (kfl_stabi_nsa >= 1) then
           ! 2. compute subscales
           !
           if (kfl_taudi_nsa < 5) then   ! diagonal tau        
              
              call nsa_monyos_diag_cdr(&
                   ielem,igaus,xresi,xsube,xortp,chale,hleng,xunkn,gunkn,gpres,gvisc,dvelo,velmo,&
                   xvelo,xvisc,xdith,sound,xlade,xldve,xdtix,taudi,qufac,&
                   conme,difme)
              
           else if (kfl_taudi_nsa >= 5) then   ! non-diagonal tau

              !
              ! shock capturing using frequencies
              !
              mfreq = mfreq_nsa
              call nsa_shofreq(kapsh,grasc,xresi(1,igaus),chale,xvelo(1,igaus),velmo(igaus), &
                   10000000000000.0_rp*zensa,mfreq)

              !!if (ittim < 2) mfreq = 2

              call nsa_monyos_cdr(&
                   pelty,pgaus,dvolu,ielem,igaus,xresi,xsube,xortp, &
                   hleng,xunkn,gunkn,gpres,gvisc,dvelo,velmo, &
                   xvelo,xvisc,xdith,sound,xlade,xldve,xdtix,taudi,qufac, &
                   conme,difme,xtime,mfreq)
           end if

        end if

        ! xsube(..,..,1) = new subscale, coming from nsa_monyos
        ! xsube(..,..,2) = subscale of the last time step, coming from nsa_gauvalxy

        
        do idime=1,ndime
           umosg_nsa(idime,ielem,igaus,2) = umosg_nsa(idime,ielem,igaus,1) 
           umosg_nsa(idime,ielem,igaus,1) = xsube(idime,igaus,1)
        end do

        ! momentum equation viscosity
        difeq(1:ndime) = xvisc(igaus) / xunkn(ndime+1,igaus,1)   
        ! continuity (no viscosity terms)
        difeq(ndime+1) = 0.0_rp       
        ! energy (adding thermal diffusion)
        difeq(ndime+2) = xdith(igaus) / xunkn(ndime+1,igaus,1) / cpcoe_nsa 

        
        densg_nsa(ielem,igaus,2) = densg_nsa(ielem,igaus,1)
        enesg_nsa(ielem,igaus,2) = enesg_nsa(ielem,igaus,1)
        densg_nsa(ielem,igaus,1) = xsube(ndime+1,igaus,1)
        enesg_nsa(ielem,igaus,1) = xsube(ndime+2,igaus,1)

        ! 3. shock capturing: compute shock capturing metrics (cartd is a dummy argument)
        !

        call nsa_shocapxy(1_ip,&
             kapsh,grasc,&
             xresi(1,igaus),cartd(1,1,igaus),chale,xvelo(1,igaus),velmo(igaus),taudi,difeq,&
             1000.0_rp*zensa,qufac,shmet,shote)
        
        do inode=1,pnode
           ipoin= lnods(inode,ielem)
           xshai= elmar(pelty)%shape(inode,igaus)

           ! 3. shock capturing: compute numerical diffusion term  (now use cartd)
           !
           shote = 0.0_rp
           call nsa_shocapxy(2_ip, &
                kapsh,grasc, &
                xresi(1,igaus),cartd(1,inode,igaus),chale,xvelo(1,igaus),velmo(igaus),taudi,difeq, &
                1000.0_rp*zensa,qufac,shmet,shote)

           ! 4. compute the local rhs and assemble it to the global rhs (scatter)
           !

           if (kfl_diagi_nsa == 0 )  then
!!$print*
!!$print*,ielem,inode
              call nsa_elerhs_cdr(&
                   elrhs,            inode,igaus,ndaux,dvolu(igaus),             xvelo(1,igaus),xunkn, &
                   xconv(1,1,1,igaus),xdiff(1,1,1,1,igaus),dconv(1,1,igaus),ddiff(1,1,1,1,igaus), &
                   gunkn(1,1,igaus),xshai,cartd(1,1,igaus),hesma(1,1,1,igaus), &
                   xsube,shote)
           else if (kfl_diagi_nsa == 1)  then
              call nsa_elerhs_matrix(&
                   elrhs,elmat,diago,taudi,pnode,inode,igaus,ndaux,dvolu(igaus),elunk(1,1,1),xvelo(1,igaus),xunkn, &
                   xconv(1,1,1,igaus),xdiff(1,1,1,1,igaus),dconv(1,1,igaus),ddiff(1,1,1,1,igaus), &
                   gunkn(1,1,igaus),xshai,cartd(1,1,igaus),hesma(1,1,1,igaus), &
                   xsube,shote)
           end if

           asfac= dvolu(igaus) * xshai / vmass(ipoin)

           umoss_nsa(    1,ipoin,2) = umoss_nsa(    1,ipoin,2) + asfac * umosg_nsa(    1,ielem,igaus,1)
           umoss_nsa(    2,ipoin,2) = umoss_nsa(    2,ipoin,2) + asfac * umosg_nsa(    2,ielem,igaus,1)
           if (ndime == 3) &
                umoss_nsa(ndime,ipoin,2) = umoss_nsa(ndime,ipoin,2) + asfac * umosg_nsa(ndime,ielem,igaus,1)
           denss_nsa(      ipoin,2) = denss_nsa(      ipoin,2) + asfac * densg_nsa(      ielem,igaus,1)
           eness_nsa(      ipoin,2) = eness_nsa(      ipoin,2) + asfac * enesg_nsa(      ielem,igaus,1)
           frequ_nsa(ipoin) = frequ_nsa(ipoin) + asfac * mfreq
        end do

     end do elemental_gauss_points_monyos_scatter

     !
     ! Assembly
     !
     call assrhs(&
          ndofn_nsa,pnode,lnods(1,ielem),elrhs,rhsid)
     if (kfl_diagi_nsa == 1) call assrhs(&
          ndofn_nsa,pnode,lnods(1,ielem),diago,vdiag_nsa)

  end do elements_loop

!!$stop

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

end subroutine nsa_elecdr
