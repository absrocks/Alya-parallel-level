subroutine ruqope(ndime,ngaus,posgp,weigp,ierro)

!-----------------------------------------------------------------------
!
!     This routine sets up the integration constants of open
!     integration rules for brick elements:
! 
!          NDIME = 1             NDIME = 2             NDIME = 3
! 
!      NGAUS  EXACT POL.     NGAUS  EXACT POL.     NGAUS  EXACT POL. 
!      -----  ----------     -----  ----------     -----  ----------
!        1      q1           1 x 1     q1          1x1x1     q1	
!        2      q3           2 x 2     q3          2x2x2     q3   
!        3      q5           3 x 3     q5          3x3x3     q5
!        4      q7           4 x 4     q7          4x4x4     q7
! 
!-----------------------------------------------------------------------
  
  use      def_kintyp
  implicit none

  integer(ip), intent(in)  :: ndime,ngaus
  integer(ip), intent(out) :: ierro
  real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
  real(rp)                 :: posgl(4),weigl(4)
  integer(ip)              :: nlocs,igaus,ilocs,jlocs,klocs
  
  ierro=0

  if(ndime==2) then
     if( ngaus > 16 ) then
        ierro = 1
        return
     end if
     nlocs=nint(sqrt(real(ngaus,rp)))
  else if(ndime==3) then
     nlocs=nint(real(ngaus,rp)**(1.0_rp/3.0_rp))
  end if
  if( ndime == 2 ) then
     if( nlocs*nlocs /= ngaus ) then
        ierro = 1 
        return
        !call runend('RUQOPE: WRONG INTEGRATION RULE')
     end if
  else if( ndime == 3 ) then
     if( nlocs*nlocs*nlocs /= ngaus ) then
        ierro = 1 
        return
        !call runend('RUQOPE: WRONG INTEGRATION RULE')
     end if     
  end if

  if(nlocs==1) then
     posgl(1)=0.0_rp
     weigl(1)=2.0_rp
  else if(nlocs==2) then
     posgl(1)=-0.577350269189626_rp
     posgl(2)= 0.577350269189626_rp
     weigl(1)= 1.0_rp
     weigl(2)= 1.0_rp
  else if(nlocs==3) then
     posgl(1)=-0.774596669241483_rp
     posgl(2)= 0.0_rp
     posgl(3)= 0.774596669241483_rp
     weigl(1)= 0.555555555555556_rp
     weigl(2)= 0.888888888888889_rp
     weigl(3)= 0.555555555555556_rp
  else if(nlocs==4)  then
     posgl(1)=-0.861136311594053_rp
     posgl(2)=-0.339981043584856_rp
     posgl(3)= 0.339981043584856_rp
     posgl(4)= 0.861136311594053_rp
     weigl(1)= 0.347854845137454_rp
     weigl(2)= 0.652145154862546_rp
     weigl(3)= 0.652145154862546_rp
     weigl(4)= 0.347854845137454_rp
  else
     ierro=1
  end if

  if(ndime==2) then
     igaus=0
     do ilocs=1,nlocs
        do jlocs=1,nlocs
           igaus=igaus+1
           weigp(  igaus)=weigl(ilocs)*weigl(jlocs)
           posgp(1,igaus)=posgl(ilocs)
           posgp(2,igaus)=posgl(jlocs)
        end do
     end do

  else if(ndime==3) then
     igaus=0
     do ilocs=1,nlocs
        do jlocs=1,nlocs
           do klocs=1,nlocs
              igaus=igaus+1
              weigp(  igaus)=weigl(ilocs)*weigl(jlocs)*weigl(klocs)
              posgp(1,igaus)=posgl(ilocs)
              posgp(2,igaus)=posgl(jlocs)
              posgp(3,igaus)=posgl(klocs)
           end do
        end do
     end do
  end if

end subroutine ruqope
