function fvecdo(n, a, b)
  use def_kintyp
  implicit none
  integer(ip) :: n,i
  real(rp)    :: a(n), b(n), fvecdo
  
  fvecdo=0.0_rp
  do i =  1, n
     fvecdo = fvecdo + a(i)*b(i)
  end do
  
end function fvecdo
