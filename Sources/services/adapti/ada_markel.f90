subroutine ada_markel(keywo,idele,imark)
!-----------------------------------------------------------------------
!****f* adapti/ada_markel
! NAME 
!    ada_markel
! DESCRIPTION
!    This routine marks the element for division
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain

  use      def_adapti
  implicit none
  integer(ip) , intent(in) :: idele
  character(4), intent(in) :: keywo
  integer(ip)              :: imark,ielem,ifacs

  integer(ip)              :: iasva,ierre,inodb,nasva,pnodb
  real(rp)                 :: esnor


  if (kfl_redgr_ada < 10) return

  if (keywo == 'FACES') then
     
     ifacs= idele
     pnodb= lfacs(mnodb+3,ifacs)
     ielem= lfacs(mnodb+3,ifacs)

     do inodb= 1,pnodb

        ! FALTA VER ESO PARA LOS EDGES...

     end do


  else if (keywo == 'CELLS') then

     ielem= idele
     nasva= ndime+2     

     esnor = 0.0_rp
     do iasva= 1,ndime
        ierre= (ielem-1)*nasva+iasva
        esnor = esnor +  erres(ierre) * erres(ierre)
     end do
     
     esnor= sqrt(esnor)
     
  end if

  imark= 0
  if (esnor > 0.1) imark= 1

end subroutine ada_markel
