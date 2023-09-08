subroutine renum0(lpont,lpntn,lword,nelem,nnode,npoin,lures,iwork)
!-----------------------------------------------------------------------
!                 
! This routine sets up array LPNTN if node renumbering for
! profile minimization is desired
!                
!-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  integer(ip) :: lword,nelem,nnode,npoin
  integer(ip) :: lures
  integer(ip) :: lpont(nnode,nelem), lpntn(npoin), iwork(lword)
  logical(lg) :: stamp
  integer(ip) :: n1,n2,n3,n4,ipoin,nposi,ns
!
! Computes the connections between degrees of freedom
!
  n1 = 1
  call renum1(lpont,iwork(n1:),nnode,npoin,nelem,nposi)
!
! Auxiliar memory for renumbering
!
  n2 = n1 + nposi          !  nodal connections
  n3 = n2 + npoin          !  nstart
  n4 = n3 + npoin - n1     !  level
  if(n4>lword) then
     write(lures,900) n4, lword
     call runend('RENUM0: INSUFFICIENT WORK SPACE')
  end if
!
! Renumbering routines
!
  call renum2(iwork(n1:),iwork(n2:),iwork(n3:),lpntn,npoin,ns)
  call renum3(iwork(n1:),iwork(n2:),iwork(n3:),lpntn,npoin,ns)
!
! Obtain the inverse of array lpntn
!
  do ipoin=1,npoin
     iwork(lpntn(ipoin))=ipoin
  end do
!
! Print renumbered nodes
!
  stamp=.false.
  if(stamp) then
     write(lures,901)
     do ipoin=1,npoin
        write(lures,902) ipoin,lpntn(ipoin),ipoin,iwork(ipoin)
     end do
  end if
!
! Formats
!
900 format(//,5x,67('#'),//,5x,&
         & '>>>> RENUMBERING MODULE REQUIRES MORE WORKING SPACE:',/,&
         & 11x,'*** NUMBER OF REQUIRED  INTEGERS :',i10,/,&
         & 11x,'*** NUMBER OF ALLOCATED INTEGERS :',i10   )
901 format(//,5x,67('#'),//,5x,&
         & '>>>> RENUMBERED NODES :',//,&
         & 11x,'OLD',10x,'NEW',10x,'NEW',10x,'OLD',/)
902 format( 11x,   i10,3x,   i10,3x,   i10,3x,i10)
  
end subroutine renum0
