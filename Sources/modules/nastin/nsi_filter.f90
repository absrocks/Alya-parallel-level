subroutine nsi_filter()
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_filter
  ! NAME 
  !    nsi_filter
  ! DESCRIPTION
  !    Filter for postprocess
  ! USES
  ! USED BY
  !    nsi_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  implicit none
  integer(ip) :: ipoin,idime,ppoin
  real(rp)    :: u

  ppoin = 0

  select case ( kfl_filte )

  case (1_ip) 
     !
     ! Cut from above
     !
     do ipoin = 1,npoin
        u = 0.0_rp
        do idime = 1,ndime
           u = u + veloc(idime,ipoin,1) * veloc(idime,ipoin,1)
        end do
        u = sqrt(u)
        if( u < 0.4_rp ) then
           ppoin = ppoin + 1
           gefil(ipoin) = ppoin
        end if
     end do

  case (2_ip)
     !
     ! Coordinates
     !
     do ipoin = 1,npoin
        if(coord(1,ipoin)<0.0_rp) then
           ppoin = ppoin + 1
           gefil(ipoin) = ppoin
        end if
     end do

  case (3_ip)
     !
     ! Coordinates
     !
     do ipoin = 1,npoin
        if(coord(1,ipoin)>4.0_rp)  then
           ppoin = ppoin + 1
           gefil(ipoin) = ppoin
        end if
     end do

  case default

     call runend('NSI_FILTER: THIS FILTER WAS NOT CODED')

  end select

end subroutine nsi_filter
