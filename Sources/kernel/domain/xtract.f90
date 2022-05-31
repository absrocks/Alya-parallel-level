subroutine xtract
  use def_master   ! general global variables
  use def_domain   ! geometry information

  implicit none

  integer(ip)   :: ielem,inode,ipoin,jpoin,idime,pelty,pnode,kinsi,nelex,ipoix,&
       ielex,npoix,kfl_ielex,kfl_ielem,lnoxe(50)
  integer(ip) , pointer :: lnxtr(:),lnodx(:,:),ipxtr(:),kpxtr(:)
  real(rp)      :: vradi,vorig(3),vdist

  !
  ! THIS SUBROUTINE ONLY RUNS IN SEQUENTIAL!!!!!!!!
  !

  ! the central element is kfl_ielem, computed somewhere else (such as *_INIVAR) or set here
  
  vradi= 1.0_rp   ! 
  kfl_ielem= 55748  ! element origin

!!  vradi= 5.0_rp   ! rayon
!!  kfl_ielem= 1  ! element origine

  ! moh-daniel
  !vradi= 1.5_rp   ! rayon
  !kfl_ielem= 3875068  ! element origine 

  pelty= ltype(kfl_ielem)
  pnode= nnode(pelty)  
  vorig= 0.0_rp  
  do inode=1,pnode
     ipoin= lnods(inode,kfl_ielem)
     do idime=1,ndime
        vorig(idime)= vorig(idime) + coord(idime,ipoin)/real(pnode,rp)
     end do
  end do

  allocate(lnxtr(nelem))
  allocate(ipxtr(npoin))
  lnxtr= 0
  ipxtr= 0

  write(993,*) 'TYPE'
  nelex= 0
  npoix= 0
  do ielem=1,nelem
     pelty= ltype(ielem)
     pnode= nnode(pelty)
     kinsi= 0
     do inode=1,pnode
        ipoin= lnods(inode,ielem)
        vdist= 0.0_rp
        do idime=1,ndime
           vdist= vdist + (coord(idime,ipoin)-vorig(idime))*(coord(idime,ipoin)-vorig(idime))
        end do
        vdist= sqrt(vdist)
        if (vdist < vradi) kinsi=kinsi+1
     end do
     if (kinsi==pnode) then
        nelex= nelex+1
        if (ielem == kfl_ielem) kfl_ielex= nelex
        lnxtr(ielem) = pnode
        write(993,*) pelty
        do inode=1,pnode
           ipoin= lnods(inode,ielem)
           if (ipxtr(ipoin)== 0) then
              npoix= npoix+1
              ipxtr(ipoin)= npoix
           end if
        end do
     end if
  end do
  write(993,*) 'END_TYPE'

  write(6,*) 
  write(6,*) 
  write(6,*) '+|------------> DEBUGGING:'
  write(6,*) '+|------------> nelex & nelem  =  ', nelex, nelem
  write(6,*) '+|------------> npoix & npoin  =  ', npoix, npoin
  write(6,*) '+|------------> center  =  ', vorig(1:ndime)
  write(6,*) '+|------------> kfl_ielem=  ', kfl_ielem
  write(6,*) '+|------------> kfl_ielex=  ', kfl_ielex 
  write(6,*) 

  allocate(kpxtr(npoix))
  kpxtr= 0

  write(993,*) 'ELEMENTS'
  ielex= 0
  do ielem= 1,nelem
     if (lnxtr(ielem) > 0) then
        pnode= lnxtr(ielem)
        ielex= ielex+1
        do inode= 1,pnode
           ipoin= lnods(inode,ielem)
           lnoxe(inode)= ipxtr(ipoin)
           kpxtr(ipxtr(ipoin))= ipoin
        end do
        write(993,300) ielex,lnoxe(1:pnode)
        !write(995,300) ielem, '   2'
     end if
  end do
  write(993,*) 'END_ELEMENTS'

  write(993,*) 'COORDINATES'
  do ipoix= 1,npoix
     write(993,200) ipoix,coord(1:ndime,kpxtr(ipoix))
  end do
  write(993,*) 'END_COORDINATES'

  fiber => xfiel(1) % a(:,:,1)
  do ipoix= 1,npoix
     write(994,200) ipoix,fiber(1:ndime,kpxtr(ipoix))
  end do


  call runend('XTRACT: DEBUGGING SUBROUTINE, COMMENT ITS CALL IN READOM AND RECOMPILE!!!!')

200 format(i8,6(f14.8))
300 format(40(2x,i8))
end subroutine xtract
