subroutine elsest_inbyte(tomem,rbyte,lbyte)
  !
  ! Give unit in bytes and conversion factor
  !
  use def_elsest
  implicit none
  integer(ip),  intent(in)  :: tomem
  real(rp),     intent(out) :: rbyte
  character(2), intent(out) :: lbyte

  if(tomem>=1024*1024*1024) then
     rbyte=1024.0_rp*1024.0_rp*1024.0_rp
     lbyte='Gb'
  else if(tomem>=1024*1024) then
     rbyte=1024.0_rp*1024.0_rp
     lbyte='Mb'
  else if(tomem>=1024) then
     rbyte=1024.0_rp
     lbyte='kb'
  else
     rbyte=1.0_rp
     lbyte=' b'
  end if

end subroutine elsest_inbyte
