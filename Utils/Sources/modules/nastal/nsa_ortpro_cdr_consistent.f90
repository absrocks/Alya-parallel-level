subroutine nsa_ortpro_cdr_consistent
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
       shmet(ndime,ndime,ndofn_nsa),shote(ndofn_nsa),elmat(nevat_nsa,nevat_nsa),xshap(mnode)
       
  integer(ip) :: ielem,inode,jnode,jdofn,idofn,itott,idime,jdime,igaus,pelty,pnode,pgaus,plapl,pface, &
       pevat,ievat,jevat,ipoin,kpoin,jpoin,ndaux,mfreq,izdom

  do ipoin = 1,npoin
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
     pevat = ndofn_nsa*pnode

     
     ! 1. Gather
     call nsa_gacons_ortpro(&
          ielem,pnode,elcod,elunk,elsub,elcon,eldif,elvel,elpre,eltem,elthe,elvis,eldtt)
     
     ! hleng and tragl at center of gravity
     call elmlen(ndime,pnode,elmar(pelty)%dercg,tragl,elcod,hnatu(pelty),hleng)

     ! compute chale: stream-wise and cross-wise lengths
     chale(1) = hleng(ndime)      ! smallest
     chale(2) = hleng(1)          ! largest


     do ievat=1,nevat_nsa
        elrhs(ievat)= 0.0_rp
        do jevat=1,nevat_nsa
           elmat(ievat,jevat)= 0.0_rp
        end do
     end do
     
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do inode=1,pnode
           xshap(inode) = elmar(pelty)%shape(inode,igaus) 
        end do

        do inode=1,pnode
           ievat= (inode-1) * ndofn_nsa + ndime + 1
           do idime=1,ndime
              ! POSAR DIFUSIÓ !!!!!!!!!!!!!!!!!!!!!!!!!!!! 
              elrhs(ievat) = elrhs(ievat) + dvolu(igaus) * xshap(inode) * ( &
                   conve_nsa(idime) * gunkn(ndime+1,idime,igaus) &
                   ! + conve_nsa(idime) * gsube(ndime+1,idime,igaus) &
                   )
           end do
           do jnode=1,pnode
              do idofn=1,ndofn_nsa
                 ievat= (inode-1) * ndofn_nsa + idofn
                 jevat= (jnode-1) * ndofn_nsa + idofn
                 elmat(ievat,jevat) = elmat(ievat,jevat) + dvolu(igaus) * xshap(inode) * xshap(jnode)
              end do
           end do
        end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$        do inode=1,pnode
!!$           ipoin= lnods(inode,ielem)
!!$           xshai= elmar(pelty)%shape(inode,igaus)
!!$           asfac = dvolu(igaus) * xshai / vmasc(ipoin)
!!$           do idime=1,ndime
!!$              ! POSAR DIFISIÓ !!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!$              ortpr_nsa(ndime+1,ipoin) = ortpr_nsa(ndime+1,ipoin) + asfac * ( &
!!$                   conve_nsa(idime) * gunkn(ndime+1,idime,igaus) &
!!$                  ! + conve_nsa(idime) * gsube(ndime+1,idime,igaus) &
!!$                   )
!!$              end do
!!$        end do
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     end do elemental_gauss_points
!!$do inode=1,pnode
!!$ievat= (inode-1) * ndofn_nsa + ndime + 1
!!$print*
!!$print*,'elrhs',ielem,inode,elrhs(ievat)
!!$end do

!!$do inode=1,pnode
!!$do jnode=1,pnode
!!$ievat= (inode-1) * ndofn_nsa
!!$jevat= (jnode-1) * ndofn_nsa
!!$print*
!!$print*, 'elmat', ielem,inode,jnode,elmat(ievat+1,jevat+1),elmat(ievat+1,jevat+2),elmat(ievat+1,jevat+3),elmat(ievat+1,jevat+4)
!!$print*, 'elmat', ielem,inode,jnode,elmat(ievat+2,jevat+1),elmat(ievat+2,jevat+2),elmat(ievat+2,jevat+3),elmat(ievat+2,jevat+4)
!!$print*, 'elmat', ielem,inode,jnode,elmat(ievat+3,jevat+1),elmat(ievat+3,jevat+2),elmat(ievat+3,jevat+3),elmat(ievat+3,jevat+4)
!!$print*, 'elmat', ielem,inode,jnode,elmat(ievat+4,jevat+1),elmat(ievat+4,jevat+2),elmat(ievat+4,jevat+3),elmat(ievat+4,jevat+4)
!!$end do
!!$end do

     call assmat(&
          solve(1)%ndofn,pnode,pevat,solve(1)%nunkn,&
          ielem,solve(1)%kfl_algso,lnods(1,ielem),elmat,amatr)



     call assrhs(&
          ndofn_nsa,pnode,lnods(1,ielem),elrhs,rhsid)

  end do elements_loop

!!$do ipoin = 1,npoin
!!$ievat = (ipoin-1) * ndofn_nsa
!!$print*
!!$print*,'rhsid',ipoin,rhsid(ievat+1),rhsid(ievat+2),rhsid(ievat+3),rhsid(ievat+4)
!!$end do

!!$do ipoin = 1,npoin
!!$do jpoin=1,npoin
!!$ievat= (ipoin-1) * ndofn_nsa
!!$jevat= (jpoin-1) * ndofn_nsa
!!$print*
!!$print*, 'amatr',ipoin,jpoin,amatr(ievat+1,jevat+1),amatr(ievat+1,jevat+2),amatr(ievat+1,jevat+3),amatr(ievat+1,jevat+4)
!!$print*, 'amatr',ipoin,jpoin,amatr(ievat+2,jevat+1),amatr(ievat+2,jevat+2),amatr(ievat+2,jevat+3),amatr(ievat+2,jevat+4)
!!$print*, 'amatr',ipoin,jpoin,amatr(ievat+3,jevat+1),amatr(ievat+3,jevat+2),amatr(ievat+3,jevat+3),amatr(ievat+3,jevat+4)
!!$print*, 'amatr',ipoin,jpoin,amatr(ievat+4,jevat+1),amatr(ievat+4,jevat+2),amatr(ievat+4,jevat+3),amatr(ievat+4,jevat+4)
!!$end do
!!$end do

!!$do ipoin = 1,npoin
!!$do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
!!$jpoin = c_dom(izdom)
!!$print*
!!$print*,'amatr',ipoin,jpoin,amatr(izdom)
!!$end do
!!$end do


  call solver(rhsid,unkno,amatr,pmatr)

  do ipoin = 1,npoin
     do idofn=1,ndofn_nsa 
        ievat = (ipoin-1) * ndofn_nsa + idofn
        ortpr_nsa(idofn,ipoin) = unkno(ievat)
     end do
!!$print*
!!$print*,'ortpr',ipoin,ortpr_nsa(1,ipoin),ortpr_nsa(2,ipoin),ortpr_nsa(3,ipoin),ortpr_nsa(4,ipoin)
  end do



end subroutine nsa_ortpro_cdr_consistent
