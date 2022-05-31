subroutine shape1(s,nnode,shapf,deriv,heslo,ierro)

  !-----------------------------------------------------------------------
  !
  ! This routine evaluates shape functions and their first derivates 
  ! for 1-d continuos with 2 3 & 4 nodes
  !
  !-----------------------------------------------------------------------

  use      def_kintyp
  implicit none
  integer(ip), intent(in)    :: nnode
  integer(ip), intent(inout) :: ierro
  real(rp),    intent(in)    :: s
  real(rp),    intent(out)   :: deriv(nnode), shapf(nnode), heslo(nnode)

  if(nnode==2) then     
     shapf(1) =  0.5_rp*(1.0_rp-s)
     shapf(2) =  0.5_rp*(1.0_rp+s)
     deriv(1) = -0.5_rp
     deriv(2) =  0.5_rp

  else if(nnode==3) then
     shapf(1) =  0.5_rp*s*(s-1.0_rp)
     shapf(2) =  0.5_rp*s*(s+1.0_rp)
     shapf(3) = -(s+1.0_rp)*(s-1.0_rp)
     deriv(1) =  s-0.5_rp
     deriv(2) =  s+0.5_rp
     deriv(3) = -2.0_rp*s
     if( ierro /= -1 ) then
        heslo(1) =   1.0_rp
        heslo(2) =   1.0_rp
        heslo(3) =  -2.0_rp
     end if

  else if(nnode==4) then

     shapf(1) = -9.0_rp/16.0_rp*(s+1.0_rp/3.0_rp)*(s-1.0_rp/3.0_rp)*(s-1.0_rp)
     shapf(2) =  9.0_rp/16.0_rp*(s+1.0_rp)*(s+1.0_rp/3.0_rp)*(s-1.0_rp/3.0_rp)
     shapf(3) = 27.0_rp/16.0_rp*(s+1.0_rp)*(s-1.0_rp/3.0_rp)*(s-1.0_rp)   
     shapf(4) = -27.0_rp/16.0_rp*(s+1.0_rp)*(s+1.0_rp/3.0_rp)*(s-1.0_rp)
     
     deriv(1) = -9.0_rp/16.0_rp*((s-1.0_rp/3.0_rp)*(s-1.0_rp)+(s+1.0_rp/3.0_rp)&
          *(s-1.0_rp)+(s+1.0_rp/3.0_rp)*(s-1.0_rp/3.0_rp))
     deriv(2) =  9.0_rp/16.0_rp*((s+1.0_rp/3.0_rp)*(s-1.0_rp/3.0_rp)+(s+1.0_rp)&
          *(s-1.0_rp/3.0_rp)+(s+1.0_rp)*(s+1.0_rp/3.0_rp))
     deriv(3) = 27.0_rp/16.0_rp*((s-1.0_rp/3.0_rp)*(s-1.0_rp)+ (s+1.0_rp)&
          *(s-1.0_rp)+ (s+1.0_rp)*(s-1.0_rp/3.0_rp))
     deriv(4) = -27.0_rp/16.0_rp*((s+1.0_rp/3.0_rp)*(s-1.0_rp)+(s+1.0_rp)&
          *(s-1.0_rp)+(s+1.0_rp)*(s+1.0_rp/3.0_rp)) 
    
     if( ierro /= -1 ) then
        heslo(1) =  -9.0_rp/16.0_rp*(2.0_rp*(s-1.0_rp)+2.0_rp*(s+1.0_rp/3.0_rp)+2.0_rp*(s-1.0_rp/3.0_rp))
        heslo(2) =   9.0_rp/16.0_rp*(2.0_rp*(s-1.0_rp/3.0_rp)+2.0_rp*(s+1.0_rp/3.0_rp)+2.0_rp*(s+1.0_rp))
        heslo(3) =  27.0_rp/16.0_rp*(2.0_rp*(s-1.0_rp)+2.0_rp*(s-1.0_rp/3.0_rp)+2.0_rp*(s+1.0_rp))
        heslo(4) = -27.0_rp/16.0_rp*(2.0_rp*(s-1.0_rp)+2.0_rp*(s+1.0_rp/3.0_rp)+2.0_rp*(s+1.0_rp))
     end if

  else
     ierro=1
  end if
  if( ierro == -1 ) ierro = 0

end subroutine shape1
    
