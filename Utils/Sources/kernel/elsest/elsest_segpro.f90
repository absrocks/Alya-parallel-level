subroutine elsest_segpro(&
     ithre,mnode,nelem,ndime,nnode,lnods,ltype,ltopo,lchec,&
     ipara,coord,point_x,shapt,derit,coloc,ifoun)
  !
  ! Maximum segment length strategy
  !
  !use def_elsest, only     : ip,rp,semax,lmini,lmaxi,elcod,i1p
  use def_elsest, only     : ip,rp,lmini,lmaxi,elcod,i1p
  implicit none
  integer(ip), intent(in)  :: ithre,mnode,nelem,ndime,nnode(*)
  integer(ip), intent(in)  :: lnods(mnode,nelem)
  integer(ip), intent(in)  :: ltype(*),ipara(*)
  integer(ip), intent(in)  :: ltopo(*),lchec(*)
  real(rp),    intent(in)  :: coord(ndime,*),point_x(*)
  integer(ip), intent(out) :: ifoun
  real(rp),    intent(out) :: coloc(*),shapt(*),derit(*)
  integer(ip)              :: ielem,inode,pnode,ptype,ipoin,idime,ilook
  real(rp)                 :: dista,dummr

  ielem = 0
  ifoun = 0
  do while( ifoun == 0 .and. ielem < nelem )
     ielem=ielem+1
     ptype=abs(ltype(ielem))
     if( ptype < 100 .and. ltype(ielem) > 0 ) then
        pnode=nnode(ptype)
        inode=0
        do while(inode<pnode)
           inode=inode+1
           ipoin=lnods(inode,ielem)
           dista=0.0_rp
           do idime=1,ndime
              dummr=coord(idime,ipoin)-point_x(idime)
              dista=dista+dummr*dummr
           end do
           !if(dista>lmaxi*semax(ielem)) inode=100
        end do

        ilook=1
        if(ipara(14)/=0) then
           if(lchec(ielem)/=ipara(14)) ilook=0
           if(ltype(ielem)<=0) ilook=0
        end if

        if(inode/=100.and.ilook==1) then
           do inode=1,pnode
              ipoin=lnods(inode,ielem)
              do idime=1,ndime
                 elcod(idime,inode,ithre)=coord(idime,ipoin)
              end do
           end do
           call elsest_chkelm(&
                ndime,ltopo(ptype),pnode,elcod(:,:,ithre),shapt,&
                derit,point_x,coloc,ifoun,lmini,lmaxi)
           if(ifoun>0) ifoun=ielem
        end if
     end if
  end do

end subroutine elsest_segpro
