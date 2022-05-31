subroutine nsa_gauval_ortpro(ielem,igaus,pnode,pgaus,weigh, &
     elcon,eldif,elunk,elsub,eldtt,xunkn,xdtix,elpre, &
     xconv,dconv,xdiff,ddiff,gunkn,gsube, &
     xshap,cartd,hesma,hessi,xsube,xresi, &
     sound,xpres,xtemp,xvisc,xdith,xvelo,gpres,gtemp, &
     gvisc,gvelo,dvelo,velmo,xlade,xldve,xjaci,conme,difme)
  use      def_master
  use      def_domain
  use      def_nastal
  implicit none
  integer(ip) :: ielem,igaus,idime,kdime,idofn,inode,pnode,pgaus,jdofn,jdime,itens
  real(rp)    :: &
       elunk(ndofn_nsa,mnode,ncomp_nsa),elsub(ndofn_nsa,mnode),eldtt(ndofn_nsa,mnode,2), &
       elpre(mnode),elcon(ndofn_nsa,ndofn_nsa,ndime,mnode),hessi(ntens,mnode), &
       eldif(ndofn_nsa,ndofn_nsa,ndime,ndime,mnode),ddiff(ndofn_nsa,ndofn_nsa,2,ndime), &
       xshap(mnode),cartd(ndime,mnode),hesma(ndime,ndime,mnode), &
       gunkn(ndofn_nsa,ndime),gsube(ndofn_nsa,ndime),dconv(ndofn_nsa,ndofn_nsa), &
       xconv(ndofn_nsa,ndofn_nsa,ndime),xdiff(ndofn_nsa,ndofn_nsa,ndime,ndime), &
       xsube(ndofn_nsa,mgaus,3),xresi(ndofn_nsa),xunkn(ndofn_nsa,mgaus,3),xshai, &
       xvelo(ndime),gpres(ndime),gpreo(ndime),gtemp(ndime),gvisc(ndime),gvelo(ndime,ndime),xldve(ndime), &
       xdtix(ndofn_nsa,mgaus,2),dvelo,sound,xpres,xpreo,xtemp,xvisc,xdith,velmo,xlade,dvite,velsq,enepe,dicod,visci, &
       xlimi,hunkn(ndofn_nsa,ndime,ndime),xtunk(ndofn_nsa), gtunk(ndofn_nsa,ndime), &
       xjaci(ndime,ndime),conme(ndofn_nsa,ndofn_nsa,ndime), &
       difme(ndofn_nsa,ndofn_nsa,ndime,ndime),weigh

  xlade     = 0.0_rp
  xpreo     = 0.0_rp       ! old pressure and pressure gradient: only 
  gpreo(1:ndime) = 0.0_rp  !     used as limiter when current pressure goes below zero
  do idofn= 1,ndofn_nsa
     xdtix(idofn,igaus,1) = 0.0_rp
     xdtix(idofn,igaus,2) = 0.0_rp
     do jdofn=1,ndofn_nsa
        do idime=1,ndime
           xconv(idofn,jdofn,idime) = 0.0_rp           
           ddiff(idofn,jdofn,1,idime) = 0.0_rp           
           ddiff(idofn,jdofn,2,idime) = 0.0_rp           
           do jdime=1,ndime
              xdiff(idofn,jdofn,idime,jdime) = 0.0_rp           
           end do
        end do
        dconv(idofn,jdofn) = 0.0_rp
     end do
     xunkn(idofn,igaus,1) = 0.0_rp
     xunkn(idofn,igaus,2) = 0.0_rp
     xunkn(idofn,igaus,3) = 0.0_rp
     xresi(idofn) = 0.0_rp
  end do
  do idime=1,ndime
     xldve(idime) = 0.0_rp
     gpres(idime) = 0.0_rp
     do jdime=1,ndime
        gvelo(idime,jdime) = 0.0_rp
     end do
     do idofn=1,ndofn_nsa
        gunkn(idofn,idime) = 0.0_rp
        gsube(idofn,idime) = 0.0_rp
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
        xunkn(idofn,igaus,1) = xunkn(idofn,igaus,1) &
             + xshai*elunk(idofn,inode,1)
        xunkn(idofn,igaus,3) = xunkn(idofn,igaus,3) + xshai*elunk(idofn,inode,3)
        xdtix(idofn,igaus,1) = xdtix(idofn,igaus,1) + xshai*eldtt(idofn,inode,1)     
        xdtix(idofn,igaus,2) = xdtix(idofn,igaus,2) + xshai*eldtt(idofn,inode,2)
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
     do idime=1,ndime
        gpreo(idime) = gpreo(idime) + elpre(inode) * cartd(idime,inode)
        do jdime=1,ndime
           hesma(idime,jdime,inode) = hessi(nindx_nsa(idime,jdime),inode)
         end do
        do idofn = 1,ndofn_nsa
           gunkn(idofn,idime) = gunkn(idofn,idime) &
                + cartd(idime,inode)*elunk(idofn,inode,1)
           gsube(idofn,idime) = gsube(idofn,idime) + cartd(idime,inode)*elsub(idofn,inode)
           do jdime=1,ndime
              hunkn(idofn,idime,jdime) = hunkn(idofn,idime,jdime) + hesma(idime,jdime,inode) * elunk(idofn,inode,1)
           end do
        end do
        xlade        = xlade        + hessi(idime,inode) * elunk(ndime+1,inode,1) 
        xldve(idime) = xldve(idime) + hessi(idime,inode) * elunk(ndime+1,inode,1)
     end do
  end do

  do idofn=1,ndofn_nsa
     xtunk(idofn) = xunkn(idofn,igaus,1) + xsube(idofn,igaus,2)
     do idime=1,ndime
        gtunk(idofn,idime) = gunkn(idofn,idime) + gsube(idofn,idime)
     end do
  end do
  if (kfl_track_nsa == 0) then
     do idofn=1,ndofn_nsa
        xtunk(idofn) = xunkn(idofn,igaus,1)
        do idime=1,ndime
           gtunk(idofn,idime) = gunkn(idofn,idime)
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
             - xtunk(idime) * gtunk(ndime+1, jdime) &
             / xtunk(ndime+1)) &
             / xtunk(ndime+1)
     end do
  end do

  do idime=1,ndime
     do jdime=1,ndime
        gpres(idime) = gpres(idime) + gvelo(jdime,idime) * xtunk(jdime) &
             + xvelo(jdime) * gtunk(jdime,idime)
     end do
     gpres(idime) = rgasc_nsa * (gtunk(ndime+2,idime) &
          - 0.5_rp * gpres(idime)) / cvcoe_nsa
  end do

  velmo= sqrt(velsq)
  xpres = rgasc_nsa * (xtunk(ndime+2) - 0.5_rp * xpres &
       / xtunk(ndime+1))  /   cvcoe_nsa

  if (xpres .lt. zensa) then
     xpres = xpreo
     do idime=1,ndime
        gpres(idime)= gpreo(idime)
     end do
  end if

  xtemp = xpres / xtunk(ndime+1) / rgasc_nsa

  dvite = 0.0_rp
  xvisc = 0.0_rp
  if (lawvi_nsa > 0) call nsa_lawvis(-1, 1,xvisc,xtemp,dvite)
  xdith = xvisc * cpcoe_nsa / prand_nsa
  dvelo= 0.0_rp
  do idime=1,ndime
     gtemp(idime) = (gpres(idime) - xpres * gtunk(ndime+1,idime) &
          / xtunk(ndime+1)) &
          / xtunk(ndime+1) / rgasc_nsa
     gvisc(idime) = dvite * gtemp(idime)
     dvelo = dvelo + gvelo(idime,idime)
  end do
  enepe = xtunk(ndime+2) / xtunk(ndime+1)
  visci = xvisc / xtunk(ndime+1)
  dicod = xdith / cvcoe_nsa / xtunk(ndime+1)

  do jdime=1,ndime        
     xconv(ndime+1,jdime  ,jdime)= 1.0_rp
     xconv(jdime  ,ndime+2,jdime)= rgacv_nsa
     xconv(jdime  ,ndime+1,jdime)= rgacv_nsa * 0.5_rp * velsq
     xconv(ndime+2,jdime  ,jdime)= ((1.0_rp + rgacv_nsa) * enepe - rgacv_nsa * 0.5_rp * velsq)
     xconv(ndime+2,ndime+1,jdime)= - xvelo(jdime) * ((1.0_rp + rgacv_nsa) * enepe - rgacv_nsa * velsq)
     xconv(ndime+2,ndime+2,jdime)= ((1.0_rp + rgacv_nsa) * xvelo(jdime) )
     xdiff(ndime+2,ndime+1,jdime,jdime)= (dicod-visci) * velsq - dicod * enepe
     xdiff(ndime+2,ndime+2,jdime,jdime)= dicod
     do idime=1,ndime
        xconv(idime,idime  ,jdime)= xconv(idime,idime  ,jdime) + xvelo(jdime) 
        xconv(jdime,idime  ,jdime)= xconv(jdime,idime  ,jdime) - rgacv_nsa * xvelo(idime)
        xconv(idime,jdime  ,jdime)= xconv(idime,jdime  ,jdime) + xvelo(idime)
        xconv(idime,ndime+1,jdime)= xconv(idime,ndime+1,jdime) - xvelo(idime) * xvelo(jdime)
        xconv(ndime+2,idime,jdime)= xconv(ndime+2,idime,jdime) - rgacv_nsa * xvelo(idime) * xvelo(jdime)
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

  do idime=1,ndime
     do idofn=1,ndofn_nsa
        do jdofn=1,ndofn_nsa
           xresi(idofn)= xresi(idofn) - xconv(idofn,jdofn,idime) * gunkn(jdofn,idime)
           do jdime=1,ndime
              xresi(idofn)= xresi(idofn) + ddiff(idofn,jdofn,1,idime) * gunkn(jdofn,jdime) + &
                   xdiff(idofn,jdofn,idime,jdime) * hunkn(jdofn,idime,jdime)
              conme(idofn,jdofn,idime) = conme(idofn,jdofn,idime) + &
                   xconv(idofn,jdofn,jdime) * xjaci(idime,jdime) / real(pgaus)
!!he d'incluir a difme el canvi de variables del pas a l'espai parametric


!              difme(idofn,jdofn,idime,jdime) = difme(idofn,jdofn,idime,jdime) + &
!                   xdiff(idofn,jdofn,idime,jdime) / real(pgaus)

              do kdime=1,ndime
                 difme(idofn,jdofn,idime,jdime) = difme(idofn,jdofn,idime,jdime) + &
                      xdiff(idofn,jdofn,jdime,kdime) * xjaci(idime,kdime) / real(pgaus)
              end do

           end do
        end do
     end do
  end do

  sound= sqrt(adgam_nsa * rgasc_nsa * xtemp)

end subroutine nsa_gauval_ortpro
