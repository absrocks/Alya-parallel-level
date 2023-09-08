subroutine got_elmexa(itask,pgaus,gpcod,gpvel,gprhs)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmexa
  ! NAME 
  !    got_elmexa
  ! DESCRIPTION
  !    Compute RHS of exact solution
  ! USES
  ! USED BY
  !    got_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime
  use def_gotita, only     :  kfl_exacs_got
  implicit none
  integer(ip), intent(in)  :: itask,pgaus
  real(rp),    intent(in)  :: gpcod(ndime,*),gpvel(ndime,*)
  real(rp),    intent(out) :: gprhs(ndime+1,*)
  integer(ip)              :: igaus

  if(kfl_exacs_got/=0) then
     if(itask==1) then
        call got_exacso(&
             1_ip,gpcod(1,igaus),gpvel(1,igaus),gprhs(1,igaus))
     else
        do igaus=1,pgaus
           call got_exacso(&
                2_ip,gpcod(1,igaus),gpvel(1,igaus),gprhs(1,igaus))
        end do
     end if
  end if

end subroutine got_elmexa
