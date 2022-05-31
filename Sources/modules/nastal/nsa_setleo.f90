subroutine nsa_setleo
!-----------------------------------------------------------------------
!****f* Nastal/nsa_setlbe
! NAME
!    setlbe
! DESCRIPTION
!    This routine constucts LEOBL_NSA, correspondence list of 
!    boundaries to volume elements
! USED BY
!    Domain
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal
  use      mod_memchk
  implicit none
  integer(ip) :: iboun,ielem,ielty,iblty,icoun
  integer(4)  :: istat
  
  if(kfl_paral/=0) then

     if(nboun==0) return
     
     allocate(leobl_nsa(nelem),      stat=istat)
     call memchk(zero,istat,memor_dom,'LEOBL_NSA','nsa_setlbe',leobl_nsa)
     !
     ! Compute leobl_nsa:
     !
     ! LEOBL_NSA(IFACE   , IELEM) = IBOUN . It is the correspondence list of boundaries to elemets
     !
     ! LEOBL_NSA(NFACE+1 , IELEM) = 1 or 0: With or without boundaries.
     !       
     ! ... remember that:      
     !       
     ! LBOEL(INODB   ,     IBOUN) = INODE . Node in element IELEM equal to node
     !                                      INODB in the boundary IBOUN
     ! LBOEL(NNODB+1 ,     IBOUN) = IELEM . Element to which IBOUN is connected
     !
     do ielem=1,nelem
        ielty=ltype(ielem)
        allocate(leobl_nsa(ielem)%l(nface(ielty)+1),stat=istat)
        call memchk(zero,istat,memor_dom,'LEOBL_NSA','nsa_setlbe',leobl_nsa(ielem)%l)          
        leobl_nsa(ielem)%l(1:nface(ielty)+1) = 0   ! no boundaries in ielem
     end do
     do iboun=1,nboun
        iblty                                = ltypb(iboun)
        ielem                                = lelbo(iboun)
        if(ielem>nelem) then
           write(*,*) 'error',ielem,iboun
           stop
        end if
        ielty                                = ltype(ielem)
        icoun                                = leobl_nsa(ielem)%l(nface(ielty)+1)+1
        leobl_nsa(ielem)%l(1:nface(ielty)+1) = icoun
        leobl_nsa(ielem)%l(icoun)            = iboun
     end do

  end if
  
end subroutine nsa_setleo
