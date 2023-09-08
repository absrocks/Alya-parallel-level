subroutine rad_heatfl(gradt)
  !------------------------------------------------------------------------
  !
  ! Case of a scalar
  !
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  implicit none
  real(rp),    intent(out), target :: gradt(ndime,npoin)
  integer(ip)                      :: ipoin,idime,inode,ielem
  integer(ip)                      :: pnode,pelty,pmate,dummi,jnode
  real(rp)                         :: detjm,gpvol,cartc(ndime,mnode) 
  real(rp)                         :: eltem(mnode,2),elcod(ndime,mnode)
  real(rp)                         :: eledd(mnode),elvel(ndime,mnode)
  real(rp)                         :: xjaci(9),xjacm(9),gpdif,gpcon,dummr
  real(rp)                         :: gptem(mnode),fact1,fact2

!!$  if(kfl_paral/=0) then
!!$     !
!!$     ! Initialization
!!$     !
!!$     do ipoin=1,npoin
!!$        do idime=1,ndime
!!$           gradt(idime,ipoin)=0.0_rp
!!$        end do
!!$     end do
!!$     !
!!$     ! Loop over elements
!!$     !
!!$     elements: do ielem=1,nelem
!!$        pelty=ltype(ielem) 
!!$        pnode=nnode(pelty)
!!$        pmate=1
!!$        if(nmate>1) pmate=lmate(ielem)
!!$        !
!!$        ! Gather vectors
!!$        !
!!$        call rad_elmgat(&
!!$             pnode,lnods(1,ielem),eltem,elvel,elcod,eledd)
!!$        !
!!$        ! Compute GPTEM
!!$        !
!!$        call gather(&
!!$             2_ip,pnode,pnode,1_ip,dummi,&
!!$             elmar(pelty)%shapc,eltem,gptem)
!!$        !
!!$        ! Loop over Gauss points (which are nodes)
!!$        !
!!$        gauss_points: do inode=1,pnode
!!$           !
!!$           !  Properties: GPDEN, GPDIF, GPGRD and GPREA
!!$           ! 
!!$           call rad_elmpro(&
!!$                ielem,pmate,pnode,1_ip,1_ip,1_ip,&
!!$                elmar(pelty)%shapc(1,inode),gpcar,gpdif,&   !!F Careful check what should be dummy
!!$                dummr,dummr,dummr)
!!$
!!$           ipoin=lnods(inode,ielem)
!!$           call elmder(&
!!$                pnode,ndime,elmar(pelty)%deric(1,1,inode),&
!!$                elcod,cartc,detjm,xjacm,xjaci)
!!$           gpvol=elmar(pelty)%weigc(inode)*detjm
!!$           if(kfl_naxis==1) gpvol=gpvol*twopi*elcod(1,inode)
!!$           !
!!$           ! Gradients
!!$           !
!!$           fact1=gpvol*gpcon
!!$           do jnode=1,pnode
!!$              fact2=fact1*eltem(jnode,1)
!!$              do idime=1,ndime
!!$                 gradt(idime,ipoin)=gradt(idime,ipoin)&
!!$                      +fact2*cartc(idime,jnode)
!!$              end do
!!$           end do
!!$        end do gauss_points
!!$     end do elements
!!$     !
!!$     ! Parall service
!!$     !
!!$     if(kfl_paral>=1) then
!!$        parr2 => gradt
!!$        call vocabu(NPOIN_REAL_2DIM,ndime,0_ip)
!!$        call Parall(400_ip)
!!$     end if
!!$     !
!!$     ! Solve diagonal system
!!$     !
!!$     do ipoin=1,npoin
!!$        do idime=1,ndime
!!$           gradt(idime,ipoin)=gradt(idime,ipoin)/vmasc(ipoin)
!!$        end do
!!$     end do
!!$  end if

end subroutine rad_heatfl
