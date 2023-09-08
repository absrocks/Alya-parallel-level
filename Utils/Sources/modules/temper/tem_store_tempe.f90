subroutine tem_store_tempe(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_store_tempe
  ! NAME 
  !    tem_store_tempe
  ! DESCRIPTION
  !    This routine performs stores the updated values of temperature
  !    after computed from the enthalpy.
  ! USED BY
  !    tem_endite (itask=1, inner loop) 
  !    tem_endste (itask=2)
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: itime,ipoin

  if( INOTMASTER ) then

     select case (itask)

     case(1_ip)
       ! 
       ! Store current temperature (:,2) <=  (:,1)
       !
       do ipoin = 1,npoin
          tempe(ipoin,2) = tempe(ipoin,1)
       end do 

     case(2_ip)
       ! 
       ! High-order temporal schemes 
       !
       if(kfl_tisch_tem==2) then
          !
          ! BDF scheme
          !
          do itime=2+kfl_tiaor_tem,4,-1
             do ipoin=1,npoin
                tempe(ipoin,itime) = tempe(ipoin,itime-1)
             end do
          end do
       end if
       ! 
       ! Store current temperature (:,3) <=  (:,1)
       !
       do ipoin=1,npoin
          tempe(ipoin,3) = tempe(ipoin,1)
       end do 

     end select

  end if

end subroutine tem_store_tempe
