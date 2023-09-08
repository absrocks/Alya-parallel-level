!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_intfor.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Save matrix to compute internal force F
!> @details Save matrix to compute internal force F:
!>          F = bu - Auu u - Aup p 
!>          F is the force of the fluid on the solid
!>          F has unit [RHO*U*^2L]
!>          sig.n = F / mass
!> @} 
subroutine nsi_intfor(Auu,Aup,bu)
  use def_kintyp
  use def_domain
  use def_elmtyp
  use def_master
  use def_nastin
  implicit none
  real(rp),    intent(in) :: Auu(ndime,ndime,*) !> Auu
  real(rp),    intent(in) :: Aup(ndime,*)       !> Aup
  real(rp),    intent(in) :: bu(ndime,*)        !> bu
  integer(ip)             :: ipoin,jzdom,idime,jdime,izdom
  integer(4)              :: istat

  if( kfl_intfo_nsi >= 1 ) then

     if( INOTMASTER ) then

        !----------------------------------------------------------------
        !
        ! Determine where force should be computed
        !
        !----------------------------------------------------------------

        do ipoin = 1,npoin
           if( intfo_nsi(ipoin) % kfl_exist /= 0 ) then
              if( kfl_intfo_nsi == 2 ) then
                 deallocate( intfo_nsi(ipoin) % Auu , stat = istat )
                 deallocate( intfo_nsi(ipoin) % Aup , stat = istat )
                 deallocate( intfo_nsi(ipoin) % bu  , stat = istat )
              end if
           end if
        end do
        !
        ! On boundary nodes
        !
        do ipoin = 1,npoin
           if( lpoty(ipoin) > 0 ) then
              intfo_nsi(ipoin) % kfl_exist = 1
           else
              intfo_nsi(ipoin) % kfl_exist = 0
           end if
        end do
        !
        ! On zone intersection with SOLIDZ
        !
        if( kfl_coupl(ID_SOLIDZ,ID_NASTIN) /= 0 ) then
           do ipoin = 1,npoin
              if( lnoch(ipoin) == NODE_CONTACT_FLUID ) intfo_nsi(ipoin) % kfl_exist = 2
           end do
        end if
        if( kfl_coupl(ID_IMMBOU,ID_NASTIN) /= 0 ) then
           do ipoin = 1,npoin
              if( lntib(ipoin) < 0 ) intfo_nsi(ipoin) % kfl_exist = 1
           end do
        end if

        !----------------------------------------------------------------
        !
        ! Save Auu, Aup and bu in INTFO_NSI(:) % 
        !        
        !----------------------------------------------------------------

        do ipoin = 1,npoin           

           if( intfo_nsi(ipoin) % kfl_exist /= 0 ) then

              jzdom = r_dom(ipoin+1)-r_dom(ipoin)
              allocate( intfo_nsi(ipoin) % Auu(ndime,ndime,jzdom) , stat = istat )
              allocate( intfo_nsi(ipoin) % Aup(ndime,jzdom)       , stat = istat )
              allocate( intfo_nsi(ipoin) % bu(ndime)              , stat = istat )
              jzdom = 0
              do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
                 jzdom = jzdom + 1
                 do idime = 1,ndime
                    do jdime = 1,ndime
                       intfo_nsi(ipoin) % Auu(jdime,idime,jzdom) = Auu(jdime,idime,izdom)
                    end do
                 end do
                 do idime = 1,ndime
                    intfo_nsi(ipoin) % Aup(idime,jzdom) = Aup(idime,izdom)
                 end do
              end do
              do idime = 1,ndime
                 intfo_nsi(ipoin) % bu(idime) = bu(idime,ipoin)
              end do
           end if
        end do
     end if

     kfl_intfo_nsi = 2

  end if

end subroutine nsi_intfor
