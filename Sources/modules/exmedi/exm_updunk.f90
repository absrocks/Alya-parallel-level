subroutine exm_updunk(itask)
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_updunk
  ! NAME 
  !    exm_updunk
  ! DESCRIPTION
  !    This routine performs several types of updates 
  ! USED BY
  !    exm_begste (itask=1)
  !    exm_begite (itask=2)
  !    exm_endite (itask=3, inner loop) 
  !    exm_endite (itask=4, outer loop) 
  !    exm_endste (itask=5)
  !***
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  !
  ! This routine performs several types of updates for the unknowns
  !
  ! REVISAR EL 9 Y EL 13, CON LO DEL VOLTS... DEBERÍA VOLARE ESO
  !
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_exmedi
  use mod_memory
  implicit none
  integer(ip), intent(in) :: itask !> case to update
  integer(ip)             :: ipoin,itime,icomp

  ! parameters defined in def_master:
  ! ITER_K   = 1
  ! ITER_AUX = 2
  ! TIME_N   = 3

  if(INOTMASTER) then

     select case (itask)

     case ( ITASK_INIUNK )
        !
        ! (:,*) <= (:,1)
        !
        do icomp = 2,memory_size(elmag,2_ip)
           do ipoin = 1,npoin
              elmag(ipoin,icomp) = elmag(ipoin,1)
           end do
        end do
        
     case( ITASK_ENDSTE )
        !
        ! (:,5) <= (:,4)
        ! (:,4) <= (:,3)
        ! (:,3) <= (:,1)
        !
        if( kfl_tiacc_exm == 2 ) then
           do itime = 2+kfl_tiacc_exm,4,-1
              do ipoin = 1,npoin
                 elmag(ipoin,itime) = elmag(ipoin,itime-1)
              end do
           end do
        end if
        do ipoin = 1,npoin
           elmag(ipoin,TIME_N) = elmag(ipoin,ITER_K)
        end do
        
    case( ITASK_BEGSTE ) 
       !
       ! (:,2) <= (:,1)
       !
       do ipoin = 1,npoin
          elmag(ipoin,2) = elmag(ipoin,1)
       end do
       
    case( ITASK_BEGITE ) 
       !
       ! UNKNO <= (:,1)
       !
       do ipoin = 1,npoin
          unkno(ipoin) = elmag(ipoin,1)
       end do

     case( ITASK_ENDINN ) 
       !
       ! UNKNO <= (:,1)
       !
       do ipoin = 1,npoin
          elmag(ipoin,1) = unkno(ipoin) 
       end do

     case( ITASK_ENDITE ) 
       !
       ! (:,2) <= (:,1)
       !
       do ipoin = 1,npoin
          elmag(ipoin,2) = elmag(ipoin,1)
       end do
       
     end select

  end if

end subroutine exm_updunk

