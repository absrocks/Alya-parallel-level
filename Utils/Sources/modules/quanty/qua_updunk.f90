subroutine qua_updunk(itask)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_updunk
  ! NAME 
  !    qua_updunk
  ! DESCRIPTION
  !    This routine performs several types of updates of the solution
  ! USED BY
  !    qua_begste (itask=1)
  !    qua_begite (itask=2)
  !    qua_endite (itask=3, inner loop) 
  !    qua_endite (itask=4, outer loop) 
  !    qua_endste (itask=5)
  !    qua_restar (itask=6)
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_quanty
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,ieiva,kk
  real(rp)                :: denom 
  character(30) :: filefin 

  if( INOTMASTER ) then

     select case (itask)

     case(1_ip)
        !
        ! Assign (n,0,*) <-- (n-1,*,*), initial guess for outer iterations
        !

        if(kfl_timei_qua==1) then

           do ipoin=1,npoin
              rhoon(ipoin,2) = rhoon(ipoin,1)
           end do

        else

           do ipoin=1,npoin
              rhoon(ipoin,2) = 0.0_rp
           end do

        endif


     case(2_ip)
        !
        ! Assign T(n,i,0) <-- T(n,i-1,*), initial guess for inner iterations
        !
        do ipoin=1,npoin
           rhoon(ipoin,1) = rhoon(ipoin,2)
        end do

     case(3_ip)
        !
        ! Assign rho(1) <-- rho(2), update of the density
        !
        !filefin = 'data\\'//'rho'//'.3D'
        !OPEN(UNIT=111, FILE=filefin,STATUS='UNKNOWN')
        !WRITE(111,*) 'rho_1  rho_2  diff'
        !DO KK=1,npoin
        !   WRITE(111,'(3E15.8)') rhoon(KK,1),rhoon(KK,2), DABS(rhoon(kk,1)-rhoon(kk,2)) 
        !ENDDO
        !close(111)

        do ipoin = 1,npoin           
           rhoon(ipoin,1) = rhoon(ipoin,2)
           rhoon(ipoin,2) = 0.0_rp
        end do

     case(4_ip)
        !
        ! Update density
        !     
        do ipoin = 1,npoin           
           rhoon(ipoin,2) = 0.0_rp
        end do
        kk = 0

        denom=1.0
        if(kfl_spher==1)  denom = 16.0_rp*atan(1.0_rp)

        do ieiva = 1,eigen_qua(1)%neiva
           do ipoin = 1,npoin
              kk = kk + 1
              rhoon(ipoin,2) = rhoon(ipoin,2) + noccupa(ieiva)*eigen(kk)*eigen(kk)/denom
           end do
        end do
        do ipoin = 1,npoin
           rhoon(ipoin,2) = mezcla * rhoon(ipoin,2) + (1.0_rp-mezcla)*rhoon(ipoin,1)
        end do


     end select

  end if

end subroutine qua_updunk

