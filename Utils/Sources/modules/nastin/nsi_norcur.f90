subroutine nsi_norcur
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_norcur
  ! NAME 
  !    nsi_norcur
  ! DESCRIPTION
  !    This routine computes normal and curvature for surface tension
  ! USES
  ! USED BY
  !    nsi_norcur
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_memchk
  use mod_postpr
  use mod_gradie
  implicit none
  integer(ip)             :: ipoin,idime
  integer(4)              :: istat
  real(rp)                :: normn

  if ( kfl_modul(ID_LEVELS) == 0 ) call runend(' NSI_NORCUR: only makes sense with LEVEL SET coupling')

  allocate(norle_nsi(ndime,npoin),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'NORLE_NSI','nsi_norcur',norle_nsi)
  allocate(curle_nsi(npoin),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'CURLE_NSI','nsi_norcur',curle_nsi)

  call gradie(fleve(1:npoin,1),norle_nsi)    ! grad LS
  
  do ipoin = 1,npoin  ! Normal = grad LS / | grad LS |
     normn = sqrt (dot_product ( norle_nsi(1:ndime,ipoin) , norle_nsi(1:ndime,ipoin) ) )
     if (normn /= 0.0_rp) then
        do idime = 1,ndime
           norle_nsi(idime,ipoin) = norle_nsi(idime,ipoin) / normn
        end do
     end if
  end do

  call divvec(norle_nsi,curle_nsi)    ! curvature = div NORMAL

end subroutine nsi_norcur
