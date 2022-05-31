subroutine rupclo(ndime,ngaus,posgp,weigp,ierro)

  !-----------------------------------------------------------------------
  ! 
  !     This routine sets up the integration constants of closed rules for
  !     PRISMS
  ! 
  !             NDIME = 3    
  ! 
  !          NGAUS  EXACT POL.
  !          -----  ----------  
  !            3       p1       
  ! 
  !-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: ndime,ngaus
  integer(ip), intent(out) :: ierro
  real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
  ierro=0
  !
  ! Area integral
  !
  if(ndime==3) then
     if(ngaus==6) then
        posgp(1,1)= 0.0_rp
        posgp(2,1)= 0.0_rp
        posgp(3,1)= 0.0_rp
        posgp(1,2)= 1.0_rp
        posgp(2,2)= 0.0_rp
        posgp(3,2)= 0.0_rp
        posgp(1,3)= 0.0_rp
        posgp(2,3)= 1.0_rp
        posgp(3,3)= 0.0_rp
        posgp(1,4)= 0.0_rp
        posgp(2,4)= 0.0_rp
        posgp(3,4)= 1.0_rp
        posgp(1,5)= 1.0_rp
        posgp(2,5)= 0.0_rp
        posgp(3,5)= 1.0_rp
        posgp(1,6)= 0.0_rp
        posgp(2,6)= 1.0_rp
        posgp(3,6)= 1.0_rp
        weigp(  1)= 1.0_rp/12.0_rp
        weigp(  2)= 1.0_rp/12.0_rp
        weigp(  3)= 1.0_rp/12.0_rp
        weigp(  4)= 1.0_rp/12.0_rp
        weigp(  5)= 1.0_rp/12.0_rp
        weigp(  6)= 1.0_rp/12.0_rp

     else if(ngaus==18) then

        posgp(1,1)= 0.0_rp
        posgp(2,1)= 0.0_rp
        posgp(3,1)= 0.0_rp
        posgp(1,2)= 1.0_rp
        posgp(2,2)= 0.0_rp
        posgp(3,2)= 0.0_rp
        posgp(1,3)= 0.0_rp
        posgp(2,3)= 1.0_rp
        posgp(3,3)= 0.0_rp

        posgp(1,4)= 0.0_rp
        posgp(2,4)= 0.0_rp
        posgp(3,4)= 1.0_rp
        posgp(1,5)= 1.0_rp
        posgp(2,5)= 0.0_rp
        posgp(3,5)= 1.0_rp
        posgp(1,6)= 0.0_rp
        posgp(2,6)= 1.0_rp
        posgp(3,6)= 1.0_rp

        ! 2nd order nodes follow GiD,
        ! just in case...

        posgp(1,7)= 0.5_rp
        posgp(2,7)= 0.0_rp
        posgp(3,7)= 0.0_rp
        posgp(1,8)= 0.5_rp
        posgp(2,8)= 0.5_rp
        posgp(3,8)= 0.0_rp
        posgp(1,9)= 0.0_rp
        posgp(2,9)= 0.5_rp
        posgp(3,9)= 0.0_rp

        posgp(1,10)= 0.0_rp
        posgp(2,10)= 0.0_rp
        posgp(3,10)= 0.5_rp
        posgp(1,11)= 1.0_rp
        posgp(2,11)= 0.0_rp
        posgp(3,11)= 0.5_rp
        posgp(1,12)= 0.0_rp
        posgp(2,12)= 1.0_rp
        posgp(3,12)= 0.5_rp

        posgp(1,13)= 0.5_rp
        posgp(2,13)= 0.0_rp
        posgp(3,13)= 1.0_rp
        posgp(1,14)= 0.5_rp
        posgp(2,14)= 0.5_rp
        posgp(3,14)= 1.0_rp
        posgp(1,15)= 0.0_rp
        posgp(2,15)= 0.5_rp
        posgp(3,15)= 1.0_rp

        posgp(1,16)= 0.5_rp
        posgp(2,16)= 0.0_rp
        posgp(3,16)= 0.5_rp
        posgp(1,17)= 0.5_rp
        posgp(2,17)= 0.5_rp
        posgp(3,17)= 0.5_rp
        posgp(1,18)= 0.0_rp
        posgp(2,18)= 0.5_rp
        posgp(3,18)= 0.5_rp

        weigp(  1)= 1.0_rp/48.0_rp
        weigp(  2)= 1.0_rp/48.0_rp
        weigp(  3)= 1.0_rp/48.0_rp
        weigp(  4)= 1.0_rp/48.0_rp
        weigp(  5)= 1.0_rp/48.0_rp
        weigp(  6)= 1.0_rp/48.0_rp

        weigp(  7)= 1.0_rp/16.0_rp
        weigp(  8)= 1.0_rp/16.0_rp
        weigp(  9)= 1.0_rp/16.0_rp
        weigp( 13)= 1.0_rp/16.0_rp
        weigp( 14)= 1.0_rp/16.0_rp
        weigp( 15)= 1.0_rp/16.0_rp

        weigp( 10)= 0.0_rp
        weigp( 11)= 0.0_rp
        weigp( 12)= 0.0_rp
        weigp( 16)= 0.0_rp
        weigp( 17)= 0.0_rp
        weigp( 18)= 0.0_rp

     else
        ierro=1
     end if
  else if(ndime==2.and.ngaus==0) then
  else
     ierro=1
  end if
!
! Errors
!
!  if(ierro==1) call runend('RUPCLO: NOT AVAILABLE INTEGRATION RULE')


end subroutine rupclo
