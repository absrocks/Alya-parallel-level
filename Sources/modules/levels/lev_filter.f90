subroutine lev_filter()
  !------------------------------------------------------------------------
  !****f* Levels/lev_filter
  ! NAME 
  !    lev_filter
  ! DESCRIPTION
  !    Filter for postprocess
  ! USES
  ! USED BY
  !    lev_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_levels
  implicit none
  integer(ip) :: ipoin

  select case ( kfl_filte )

  case (1_ip) 
     !
     ! Test: take all
     !
     do ipoin = 1,npoin
        if( fleve(ipoin,1) > -0.1 .and. fleve(ipoin,1) < 0.1 ) then
           gefil(ipoin) = 1
        end if
     end do

  case default

     call runend('LEV_FILTER: THIS FILTER WAS NOT CODED')

  end select

end subroutine lev_filter
