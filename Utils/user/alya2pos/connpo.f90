subroutine connpo(ielty,ltype,nelem,npoin,mnode,lnods)
  !-----------------------------------------------------------------------
  !
  ! NAME
  !    connpo
  ! DESCRIPTION
  !    This routine compute the node/element connectivity arrays.
  ! OUTPUT
  !    NEPOI(NPOIN) ............ # of elements connected to nodes
  !    PELPO(NPOIN+1) .......... Pointer to list of element LELPO
  !    LELPO(PELPO(NPOIN+1)) ... List of elements connected to nodes
  !
  !-----------------------------------------------------------------------
  !use def_parame
  !use def_master
  !use def_domain
  !use mod_memchk
  
  use def_kintyp
  implicit none
  integer(ip), intent(in)             :: ielty,nelem,npoin,mnode
  integer(ip), intent(in)             :: ltype(*)
  integer(ip), intent(in)             :: lnods(mnode,*)
  
  integer(ip) :: inode,ipoin,ielem,nlelp,pnode
  integer(ip)  :: istat

  
  call elmtyp()
  
  pnode=nnode(ielty)
  
  
  !
  ! Allocate memory for NEPOI and compute it
  !
  allocate(nepoi(npoin),stat=istat)
  do ipoin=1,npoin
     nepoi(ipoin)=0
  end do
  !call memchk(zero,istat,memor_dom,'NEPOI','connpo',nepoi)
  do ielem=1,nelem
     do inode=1,nnode(ltype(ielem))
        ipoin=lnods(inode,ielem)
        !write(*,*)"ipoin",ipoin
        nepoi(ipoin)=nepoi(ipoin)+1
        !write(*,*)"ipoin,nepoin(ipoin)",ipoin,nepoi(ipoin)
     end do
  end do
  !
  ! Allocate memory for PELPO and compute it
  !
  allocate(pelpo(npoin+1),stat=istat)
  do ipoin=1,npoin+1
     pelpo(ipoin)=0
  end do
  !call memchk(zero,istat,memor_dom,'PELPO','connpo',pelpo)
  pelpo(1)=1
  do ipoin=1,npoin
     pelpo(ipoin+1)=pelpo(ipoin)+nepoi(ipoin)
     !write(*,*)"ipoin,pelpo",ipoin,pelpo(ipoin)
  end do
  !
  ! Allocate memory for LELPO and construct the list
  !
  nlelp=pelpo(npoin+1)
  !write(*,*)"nlelp",nlelp
  
  allocate(lelpo(nlelp),stat=istat)

  do ipoin=1,nlelp
     lelpo(ipoin)=0
  end do
  !call memchk(zero,istat,memor_dom,'LELPO','connpo',lelpo)
  mpopo=0
  do ielem=1,nelem
     mpopo = mpopo + nnode(ltype(ielem))*nnode(ltype(ielem))
     !write(*,*)"mpopo",mpopo
     do inode=1,nnode(ltype(ielem))
        ipoin=lnods(inode,ielem)
        lelpo(pelpo(ipoin))=ielem
        !write(*,*)"lelpo(pelpo(ipoin))",lelpo(pelpo(ipoin))
        pelpo(ipoin)=pelpo(ipoin)+1
     end do
  end do
  !
  ! Recompute PELPO and maximum number of element neighbors MEPOI
  !
  pelpo(1)=1
  mepoi=-1
  do ipoin=1,npoin
     pelpo(ipoin+1)=pelpo(ipoin)+nepoi(ipoin)
     !write(*,*)"ipoin,pelpo",ipoin,pelpo(ipoin)
     mepoi=max(mepoi,nepoi(ipoin))
     !write(*,*)"mepoi",mepoi
  end do
  !
  ! Deallocate memory for temporary node/element connectivity
  !
  !call memchk(two,istat,memor_dom,'NEPOI','connpo',nepoi)
  deallocate(nepoi,stat=istat)

  !if(istat/=0) call memerr(two,'NEPOI','connpo',0_ip)

end subroutine connpo

