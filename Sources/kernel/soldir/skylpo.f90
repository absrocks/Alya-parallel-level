subroutine skylpo(ndofn,npoin,neqns,lpont,lpdof)
!-----------------------------------------------------------------------
!
!     This routine constructs the array LPDOF corresponding to the
!     matrix with NDOFN x NDOFN block components from LPONT
!
!-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  integer(ip) :: ndofn,npoin,neqns
  integer(ip) :: lpdof(neqns), lpont(npoin)
  integer(ip) :: ipoin,idofn,ieqns,ndof2,ndofh,acdof,lprev,lpon0,heigh
!
! Initial operations
!
  ndof2 = ndofn*ndofn
  ndofh = ndofn*(ndofn-1)/2
!
! First node
!
  acdof = 0
  do idofn = 1,ndofn
     lpdof(idofn) = acdof
     acdof = acdof + idofn 
  end do
!
! Rest of nodes
!
  do ipoin = 2,npoin
     lprev = lpont(ipoin-1)
     lpon0 = ndof2*lprev + (ipoin-1)*ndofh 
     heigh = (lpont(ipoin)-lprev)*ndofn
     acdof = 0
     do idofn = 1,ndofn
        ieqns = (ipoin-1)*ndofn + idofn
        lpdof(ieqns) = lpon0 + heigh*idofn + acdof
        acdof = acdof + idofn 
     end do
  end do
  
end subroutine skylpo
