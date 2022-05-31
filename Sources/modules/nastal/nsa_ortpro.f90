subroutine nsa_ortpro
  !-----------------------------------------------------------------------
  !****f* nastal/nsa_elcons
  ! NAME 
  !    nsa_elcons
  ! DESCRIPTION
  !    Conservative set (rho,U,E) equations per-element operations:
  !    1. Compute elemental matrix and RHS 
  !    2. Compute boundary contributions
  !    3. Assemble 
  ! USES
  !    nsa_...
  ! USED BY
  !    nsa_gocons
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_nastal
  implicit none

  real(rp)    :: elrhs(nevat_nsa),diago(nevat_nsa)
  real(rp)    :: &
       elunk(ndofn_nsa,mnode,ncomp_nsa),elsub(ndofn_nsa,mnode), &
       elcod(ndime,mnode),elvel(ndime,mnode), &
       elpre(mnode),eltem(mnode), &
       elcon(ndofn_nsa,ndofn_nsa,ndime,mnode), &
       eldif(ndofn_nsa,ndofn_nsa,ndime,ndime,mnode),conme(ndofn_nsa,ndofn_nsa,ndime), &
       difme(ndofn_nsa,ndofn_nsa,ndime,ndime),eldtt(ndofn_nsa,mnode,2),kapsh(ndofn_nsa),elvis(mnode),elthe(mnode)
  real(rp)    :: &
       xconv(ndofn_nsa,ndofn_nsa,ndime,mgaus),gunkn(ndofn_nsa,ndime,mgaus),grasc(ndofn_nsa,ndime), &
       gsube(ndofn_nsa,ndime,mgaus), taudi(ndofn_nsa), &
       dconv(ndofn_nsa,ndofn_nsa,mgaus),xsube(ndofn_nsa,mgaus,3),xdtix(ndofn_nsa,mgaus,2), &
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
       shmet(ndime,ndime,ndofn_nsa),shote(ndofn_nsa)
       
  integer(ip) :: ielem,inode,jdofn,idofn,itott,idime,jdime,igaus,pelty,pnode,pgaus,plapl,pface,ievat,ipoin,ndaux,mfreq

  do ipoin=1,npoin
     do idofn=1,ndofn_nsa
        ortpr_nsa(idofn,ipoin) = 0.0_rp
     end do
  end do
  
  elements_loop: do ielem=1,nelem
     
     ! Element properties and dimensions
     pelty=ltype(ielem)
     pnode=nnode(pelty)
     pgaus=ngaus(pelty)
     plapl=llapl(pelty)
     pface=nface(pelty)

     
     ! 1. Gather
     call nsa_gacons_ortpro(&
          ielem,pnode,elcod,elunk,elsub,elcon,eldif,elvel,elpre,eltem,elthe,elvis,eldtt)
     
     ! hleng and tragl at center of gravity
     call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)

     ! compute chale: stream-wise and cross-wise lengths
     chale(1) = hleng(ndime)      ! smallest
     chale(2) = hleng(1)          ! largest
     
     elemental_gauss_points: do igaus=1,pgaus
        
        call elmder(pnode,ndime,elmar(pelty)%deriv(1,1,igaus),elcod,cartd(1,1,igaus),detjm,xjacm,xjaci)
        dvolu(igaus)=elmar(pelty)%weigp(igaus)*detjm                
        hessi(1:ntens,1:mnode) = 0.0_rp

        if(plapl==1) then
           call elmhes(elmar(pelty)%heslo(1,1,igaus),hessi,ndime,pnode,ntens,&
                xjaci,d2sdx,elmar(pelty)%deriv(1,1,igaus),elcod)     
        end if
        
        ! 2. Cálculo de todo en los gauss, incluido adjunto
        call nsa_gauval_ortpro(ielem,igaus,pnode,pgaus,elmar(pelty)%weigp(igaus), &
             elcon,eldif,elunk,elsub,eldtt,xunkn,xdtix,elpre, &
             xconv(1,1,1,igaus),dconv(1,1,igaus),xdiff(1,1,1,1,igaus), &
             ddiff(1,1,1,1,igaus),gunkn(1,1,igaus),gsube(1,1,igaus), &
             elmar(pelty)%shape(1,igaus),cartd(1,1,igaus),hesma(1,1,1,igaus),hessi,xsube,xresi(1,igaus), &
             sound(igaus),xpres(igaus),xtemp(igaus),xvisc(igaus),xdith(igaus),xvelo(1,igaus), &
             gpres(1,igaus),gtemp(1,igaus), &
             gvisc(1,igaus),gvelo(1,1,igaus),dvelo(igaus),velmo(igaus),xlade(igaus),xldve(1,igaus), &
             xjaci,conme,difme)

        ! Compute ortpro_nsa, the projection on the finite element space
        do inode=1,pnode
           ipoin= lnods(inode,ielem)
           xshai= elmar(pelty)%shape(inode,igaus)

           asfac = dvolu(igaus) * xshai / vmass(ipoin)
           do idofn=1,ndofn_nsa
              do jdofn=1,ndofn_nsa
                 do idime=1,ndime   
                    ortpr_nsa(idofn,ipoin) = ortpr_nsa(idofn,ipoin) + asfac * ( &
                         xconv(idofn,jdofn,idime,igaus) * gunkn(jdofn,idime,igaus) &
                         + xconv(idofn,jdofn,idime,igaus) * gsube(jdofn,idime,igaus) &
                         )
                 end do
              end do
           end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$           do idofn=1,ndofn_nsa
!!$              do jdofn=1,ndofn_nsa
!!$                 ortpr_nsa(idofn,ipoin) = ortpr_nsa(idofn,ipoin) - &
!!$                      xshai * dconv(idofn,jdofn,igaus) * xsube(jdofn,igaus,2) * dvolu(igaus) / vmass(ipoin)
!!$                 do idime=1,ndime   
!!$                    ortpr_nsa(idofn,ipoin) = ortpr_nsa(idofn,ipoin) + ( &
!!$                         xshai * xconv(idofn,jdofn,idime,igaus) * gunkn(jdofn,idime,igaus) &
!!$                         - cartd(idime,inode,igaus) * xconv(idofn,jdofn,idime,igaus) * xsube(jdofn,igaus,2) &
!!$                         ) * dvolu(igaus) / vmass(ipoin)
!!$                 end do
!!$              end do
!!$           end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        end do

     end do elemental_gauss_points
  end do elements_loop

end subroutine nsa_ortpro
