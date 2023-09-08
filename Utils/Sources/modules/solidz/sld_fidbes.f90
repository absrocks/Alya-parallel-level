subroutine sld_fidbes
  !------------------------------------------------------------------------
  !****f* Solidz/sld_fidbes
  ! NAME 
  !    sld_fidbes
  ! DESCRIPTION
  !    Compute fibers for output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solidz
  implicit none
  integer(ip)             :: ipoin,idime,jdime
  real(rp)                :: dummr


  if( INOTMASTER ) then

     do ipoin = 1,npoin
        do idime = 1,ndime
           fibde_sld(idime,ipoin) = 0.0_rp           
           do jdime= 1,ndime
              fibde_sld(idime,ipoin) = fibde_sld(idime,ipoin) &
                   + gdepo(idime,jdime,ipoin) * fiber(jdime,ipoin)
           end do
        end do        
     end do

  end if

!  if( INOTMASTER ) then
!     do ipoin = 1,npoin
!        do idime = 1,ndime
!           fibde_sld(idime,ipoin) = 0.0_rp           
!        end do
!     end do
!     call sld_elmope(5_ip)
!  end if
!  
!  call rhsmod(ndime,    fibde_sld)
!  
!  if( INOTMASTER ) then
!     do ipoin = 1,npoin
!        dummr = 1.0_rp/vmass(ipoin)
!        do idime = 1,ndime
!           fibde_sld(idime,ipoin) = dummr * fibde_sld(idime,ipoin)
!        end do
!     end do
!  end if
  



end subroutine sld_fidbes
