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
  use      def_parame
  use      def_master
  use      def_domain
  use      def_exmedi

  implicit none
  integer(ip), intent(in) :: itask !> case to update
  integer(ip)             :: ipoin,ncomp,idofn,igate,itime

  ! parameters defined in def_master:
  ! ITER_K   = 1
  ! ITER_AUX = 2
  ! TIME_N   = 3

  if(INOTMASTER) then

     select case (itask)
        !
        ! Assign U(n,0,*) <-- U(n-1,*,*), initial guess for outer iterations
        !
     case(ITERAUX_EQ_TIMEN)
        do ipoin=1,npoin
           elmag(ipoin,ITER_AUX) = elmag(ipoin,TIME_N)
        end do
        !
        ! Assign U(n,i,0) <-- U(n,i-1,*), initial guess for inner iterations
        !
     case(ITERK_EQ_ITERAUX)
        do ipoin=1,npoin
           elmag(ipoin,ITER_K)     = elmag(ipoin,ITER_AUX)
        end do
     case(ITERK_EQ_UNKNO)
        !
        ! Assign u(n,i,j-1) <-- u(n,i,j)
        !        
        ! u(,ITER_K) <-- unkno  
        elmag(1:npoin,ITER_K) = unkno(1:npoin)

     case(ITERAUX_EQ_ITERK)
        !
        ! Assign U(n,i-1,*) <-- U(n,i,*)
        !        
        do ipoin=1,npoin
           elmag(ipoin,ITER_AUX) = elmag(ipoin,ITER_K)
        end do

     case(TIMEN_EQ_ITERK)

        if(kfl_tiacc_exm==2) then
           do itime=2+kfl_tiacc_exm,4,-1
              do ipoin=1,npoin
                 elmag(ipoin,itime) = elmag(ipoin,itime-1)
              end do
           end do
        end if
        do ipoin=1,npoin
           elmag(ipoin,TIME_N) = elmag(ipoin,ITER_K)
        end do
        
     end select

  end if

end subroutine exm_updunk

