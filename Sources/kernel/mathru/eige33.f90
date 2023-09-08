subroutine eige33(n,a,eigva,eigve,det)
  use def_kintyp, only     :  ip,rp
  use def_parame, only     :  pi
  implicit none
  integer(ip), intent(in)  :: n
  real(rp),    intent(in)  :: a(n,n)
  real(rp),    intent(out) :: eigva(n)
  real(rp),    intent(out) :: eigve(n,n)
  real(rp),    intent(out) :: det
  integer(ip)              :: i,j
  real(rp)                 :: tra,x,y
  
  if( n == 2 ) then
     tra      = a(1,1) + a(2,2)
     det      = a(1,1) * a(2,2) - a(1,2) * a(2,1)
     x        = 0.5_rp * tra
     y        = sqrt(0.25_rp*tra*tra-det)
     eigva(1) = x + det
     eigva(2) = x - det
     if( a(2,1) /= 0.0_rp ) then
        eigve(1,1) =  eigva(1) - a(2,2)
        eigve(2,1) =  a(2,1)
        eigve(1,2) =  eigva(2) - a(2,2)
        eigve(2,2) =  a(2,1)
     else if( a(1,2) /= 0.0_rp ) then
        eigve(1,1) =  a(2,1)
        eigve(2,1) =  eigva(1) - a(1,1)
        eigve(1,2) =  a(2,1)
        eigve(2,2) =  eigva(2) - a(1,1)
     else
        eigve(1,1) =  1.0_rp
        eigve(2,1) =  0.0_rp
        eigve(1,2) =  0.0_rp
        eigve(2,2) =  1.0_rp
     end if
  else
     call runend('NOT CODED')
  end if

end subroutine eige33
