!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_elsest.f90
!> @author  Guillaume Houzeaux
!> @date    27/02/2013
!> @brief   Prepare elsest for subdomains
!> @details Prepare elsest for subdomains
!> @} 
!-----------------------------------------------------------------------

subroutine dod_elsest()

  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_dodeme 
  use mod_memory
  use def_kermod
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: isubd,dummi
  real(rp)    :: dummr

  if( INOTSLAVE ) then
     call livinf(0_ip,'CONSTRUCT BIN/OCT TREE FOR EACH SUBDOMAIN USING ELSEST',0_ip)
     ielse(5) = nsubd
     call runend('DODEME: NOT CODED')
     !call elsest(&
     !     0_ip,0_ip,ielse,mnode,ndime,npoin,nelem,nnode(1:),&
     !     lnods,ltype,ltopo,coord,dummr,relse,dummi,&
     !     shape_dod,deriv_dod,dummr,dummi)
     do isubd = 1,nsubd
        current_subdomain => subdomain(isubd)
        call runend('DODEME: NOT CODED')
        !call elsest(&
        !     1_ip,isubd,ielse,mnode,ndime,current_subdomain % npoin,current_subdomain % nelem,&
        !     nnode(1:),current_subdomain % lnods,current_subdomain % ltype,&
        !     ltopo,current_subdomain % coord,dummr,relse,dummi,shape_dod,deriv_dod,dummr,dummi)
     end do
  end if

end subroutine dod_elsest
