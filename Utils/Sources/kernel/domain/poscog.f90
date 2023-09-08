subroutine poscog(poscg,ndime,smplx)

!-----------------------------------------------------------------------
!
!     This routine computes the natural coordinates of the center of
!     gravity
!
!-----------------------------------------------------------------------
  use       def_kintyp
  implicit  none
  integer(ip), intent(in)  :: smplx
  integer(ip), intent(in)  :: ndime
  real(rp),    intent(out) :: poscg(ndime)

  poscg=0.0_rp

  if(smplx==1) then
     if(ndime==1) then
        poscg(1)=0.5_rp
     else if(ndime==2) then
        poscg(1)=1.0_rp/3.0_rp
        poscg(2)=1.0_rp/3.0_rp
     else if(ndime==3) then
        poscg(1)=0.25_rp
        poscg(2)=0.25_rp
        poscg(3)=0.25_rp
     end if
  end if
  
end subroutine poscog
