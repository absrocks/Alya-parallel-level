subroutine qualit(quali,qmaxi,qmini)
  !-----------------------------------------------------------------------
  !****f* domain/qualit
  ! NAME
  !    qualit
  ! DESCRIPTION
  !    This routine computes the quality of the mesh
  !    KFL_QUALI = 1 ... GAMMA
  !              = 2 ... Aspect ratio
  !              = 3 ... Condition number
  ! OUTPUT
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_elmtyp
  use def_domain
  use def_master
  use def_kermod
  use mod_messages, only : livinf
  implicit none
  real(rp),   intent(out) :: quali(*)
  integer(ip)             :: ielem,pelty,pnode,pface,ipoin,iface,inodb
  integer(ip)             :: idime,inode,pnodb,jpoin,iedge,pblty
  real(rp)                :: alpha,qmini,qmaxi,surfa,xjaci(9),xjacm(9)
  real(rp)                :: gpdet,dista,volum,baloc(9),hmaxi,xqual,hmini
  real(rp)                :: asrad
  real(rp)                :: gpcar(ndime,mnode) 
  real(rp)                :: bocod(ndime,mnodb)
  real(rp)                :: elcod(ndime,mnode)

  call livinf(0_ip,'COMPUTE ELEMENT QUALITY',0_ip)

  if( INOTMASTER ) then        

     if( ndime == 2 ) then
        alpha = sqrt(3.0_rp) / 6.0_rp 
     else
        alpha = sqrt(6.0_rp) / 12.0_rp
     end if

     qmini =  1e7_rp
     qmaxi = -1e7_rp

     do ielem = 1,nelem

        if( lelch(ielem) /= ELHOL ) then
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           pface = nface(pelty)
           surfa = 0.0_rp
           hmaxi = 0.0_rp
           hmini = 1.0e9_rp
           !
           ! SURFA: Surface
           !
           do iface = 1,pface
              pnodb = nnodf(pelty) % l(iface)
              do inodb = 1,pnodb
                 ipoin = lnods(lface(pelty) % l(inodb,iface),ielem)
                 do idime = 1,ndime
                    bocod(idime,inodb) = coord(idime,ipoin)
                 end do
              end do
              pblty = ltypf(pelty) % l(iface)
              call bouder(&
                   pnodb,ndime,ndimb,elmar(pblty) % dercg,&               ! Cartesian derivative
                   bocod,baloc,gpdet)                                   ! and Jacobian
              surfa = surfa + gpdet * elmar(pblty) % weicg
           end do
           !
           ! VOLUM: Volume
           !
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              do idime = 1,ndime
                 elcod(idime,inode) = coord(idime,ipoin)
              end do
           end do
           call jacobi(&
                ndime,pnode,elcod,elmar(pelty) % dercg,&
                xjacm,xjaci,gpcar,gpdet)
           volum = gpdet * elmar(pelty) % weicg
           !
           ! HMAXI: Largest edge
           !
           do iedge = 1,needg(pelty) 
              ipoin = lnods(leedg(1,iedge,pelty),ielem) 
              jpoin = lnods(leedg(2,iedge,pelty),ielem) 
              dista = 0.0_rp
              do idime = 1,ndime
                 dista = dista &
                      + ( coord(idime,ipoin)-coord(idime,jpoin) ) &
                      * ( coord(idime,ipoin)-coord(idime,jpoin) )
              end do
              hmaxi = max(hmaxi,dista)
              hmini = min(hmini,dista)
           end do
           hmaxi = sqrt(hmaxi)
           hmini = sqrt(hmini)

           if( kfl_quali == 1 ) then
              !
              ! Gamma team criteria
              !
              xqual = alpha * hmaxi * surfa / (3.0_rp*abs(volum))

           else if( kfl_quali == 2 ) then
              !
              ! Aspect ratio
              !
              xqual = hmaxi/hmini

           else if( kfl_quali == 3 ) then
              !
              ! Condition number
              !
              call cndnum(ielem,xqual,asrad)
              
           else

              call runend('QUALIT: NOT CODED')

           end if

           qmini = min(qmini,xqual)
           qmaxi = max(qmaxi,xqual)
           quali(ielem) = xqual
        else
           quali(ielem) = 0.0_rp
        end if
     end do

  end if
  !
  ! AllReduce in parallel
  !
  call pararr('MIN',0_ip,1_ip,qmini)
  call pararr('MAX',0_ip,1_ip,qmaxi)

end subroutine qualit

subroutine cndnum(ielem,kappaS,asrad)
  use def_kintyp
  use def_elmtyp
  use def_parame
  use def_master
  use def_domain
  implicit none  
  integer(ip),intent(in)      :: ielem
  real(rp) ,intent(out)       :: kappaS,asrad
  integer(ip)                 :: idime,inode,tetra,ii,jj,pnode,pelty
  real(rp)                    :: invW(ndime,3),tinvW(3,ndime),invWt(3,ndime)
  real(rp)                    :: A(ndime,3),S(ndime,3),invS(ndime,3)
  real(rp)                    :: normS,norminvS,StS(ndime,3)
  real(rp)                    :: invSt(ndime,3)
  real(rp)                    :: aux(ndime,3),aux2(ndime,3)
  real(rp)                    :: detinvWt,detS,detSt
  real(rp)                    :: elcod(ndime,mnode),volum
  real(rp)                    :: xjaci(ndime,ndime),xjacm(ndime,ndime) ,detjm
  real(rp)                    :: cartd(ndime,mnode) ,sumsi
  real(rp)                    :: side1(ndime),side2(ndime),side3(ndime)
  real(rp)                    :: side4(ndime),side5(ndime),side6(ndime)

  detinvWt  =  0.0_rp
  detS      =  0.0_rp
  pnode     =  lnnod(ielem)
  pelty     =  ltype(ielem)

  if (ndime==3)then
     invW(1,1) =  1.0_rp
     invW(1,2) = -1.0_rp/3.0_rp*sqrt(3.0_rp)            !!!-5.77350269189626e-01!!!
     invW(1,3) = -1.0_rp/6.0_rp*sqrt(6.0_rp)             !!-4.08248290463863e-01 
     invW(2,1) =  0.0_rp
     invW(2,2) =  2.0_rp/3.0_rp*sqrt(3.0_rp)           !!!1.15470053837925e+00 
     invW(2,3) = -1.0_rp/6.0_rp*sqrt(6.0_rp)         !!!-4.08248290463863e-01
     invW(3,1) =  0.0_rp
     invW(3,2) =  0.0_rp
     invW(3,3) =  1.0_rp/2.0_rp*sqrt(6.0_rp)           !!!!1.22474487139159e+00
  else
     invW(1,1) =  1.0_rp
     invW(1,2) =  1.0_rp/2.0_rp              
     invW(2,1) =  0.0_rp
     invW(2,2) =  sqrt(3.0_rp)/2.0_rp           !!!1.15470053837925e+00 
  end if

  normS     =  0.0_rp
  norminvS  =  0.0_rp

  do inode=1,pnode-1
     do idime=1,ndime
        A(idime,inode)    = 0.0_rp
        S(idime,inode)    = 0.0_rp
        StS (idime,inode) = 0.0_rp
        invS(idime,inode) = 0.0_rp
     end do
  end do
  do inode=1,pnode-1
     do idime = 1,ndime
        tinvW (idime,inode)= invW(inode,idime)
     end do
  end do
  call invmtx(tinvW,invWt,detinvWt,ndime)
  ii = 0
  do inode=2,pnode
     ii = ii + 1
     do idime=1,ndime
        A(idime,ii) = coord(idime,lnods(inode,ielem)) - coord(idime,lnods(1,ielem))
     end do
  end do
  do inode =1,pnode-1
     do idime=1,ndime
        do jj = 1,pnode-1
           S(idime,inode) = S(idime,inode) + A(idime,jj) * invW(jj,inode)
        end do
     end do
  end do

  do inode=1,pnode-1
     do idime=1,ndime
        normS = normS + S(idime,inode) * S(idime,inode)  
        do jj = 1,pnode-1
           StS (idime,inode) = StS(idime,inode) + S(jj,idime)* S(jj,inode)
        end do
     end do
  end do
  normS = sqrt(normS)

  do inode=1,pnode-1
     do idime=1,ndime
        aux(idime,inode) = S(idime,inode)
     end do
  end do
  call invmtx(aux,aux2,detS,ndime)
  do inode=1,pnode-1
     do idime=1,ndime
        invS(idime,inode) = aux2(idime,inode) 
     end do
  end do

  do inode=1,pnode-1
     do idime=1,ndime
        norminvS = norminvS + invS(idime,inode) * invS(idime,inode)  
     end do
  end do

  norminvS = sqrt(norminvS)
  kappaS   = normS* norminvS
  do idime=1,ndime
     elcod(idime,1) = coord(idime,lnods(1,ielem))
  end do
  do idime=1,ndime
     elcod(idime,2) = coord(idime,lnods(2,ielem))
  end do
  do idime=1,ndime
     elcod(idime,3) = coord(idime,lnods(3,ielem))
  end do
  if(ndime==3)then
     do idime=1,ndime
        elcod(idime,4) = coord(idime,lnods(4,ielem))
     end do
  end if
  volum = 0.0_rp
  sumsi = 0.0_rp
  !!pelty = TET04
  call jacobi(&
       ndime,pnode,elcod,elmar(pelty)%dercg,&
       xjacm,xjaci,cartd,detjm)

  volum = volum + elmar(pelty)%weicg*detjm
  if( ndime == 3 ) then
     do idime=1,ndime
        side1(idime) = coord(idime,lnods(4,ielem))-coord(idime,lnods(1,ielem))
        side2(idime) = coord(idime,lnods(3,ielem))-coord(idime,lnods(1,ielem))
        side3(idime) = coord(idime,lnods(2,ielem))-coord(idime,lnods(1,ielem))
        side4(idime) = coord(idime,lnods(4,ielem))-coord(idime,lnods(2,ielem))
        side5(idime) = coord(idime,lnods(3,ielem))-coord(idime,lnods(2,ielem))
        side6(idime) = coord(idime,lnods(4,ielem))-coord(idime,lnods(3,ielem)) 
     end do
  else
     do idime=1,ndime
        side1(idime) = coord(idime,lnods(3,ielem))-coord(idime,lnods(1,ielem))
        side2(idime) = coord(idime,lnods(2,ielem))-coord(idime,lnods(1,ielem))
        side3(idime) = coord(idime,lnods(3,ielem))-coord(idime,lnods(2,ielem))
        side4(idime) = 0.0_rp
        side5(idime) = 0.0_rp
        side6(idime) = 0.0_rp
     end do
  end if

  do idime=1,ndime
     sumsi = sumsi + side1(idime)*side1(idime) + side2(idime)*side2(idime)+ &
          side3(idime)*side3(idime)+ side4(idime)*side4(idime)+             &
          side5(idime)*side5(idime)+ side6(idime)*side6(idime)
  end do

  asrad = sumsi / (volum+epsilon(1.0_rp))  


end subroutine cndnum
