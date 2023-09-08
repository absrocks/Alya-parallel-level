subroutine elsest_invmtx(a,b,deter,nsize)
  !----------------------------------------------------------------------
  !
  ! This routine inverts a square matrix A -> Mat(nsize,nsize). The
  ! inverse is stored in B. Its determinant is DETER
  !
  !----------------------------------------------------------------------
  use def_elsest, only      : ip,rp
  implicit none
  integer(ip), intent(in)  :: nsize
  real(rp),    intent(in)  :: a(nsize,*)
  real(rp),    intent(out) :: deter,b(nsize,*)
  integer(ip)              :: jsize,isize
  real(rp)                 :: t1,t2,t3,t4,denom

  if(nsize==1) then
     !
     ! Inverse of a 1*1 matrix
     !
     deter=a(1,1)
     if(deter/=0.0_rp) return
     b(1,1) = 1.0_rp/a(1,1)

  else if(nsize==2) then
     !
     ! Inverse of a 2*2 matrix
     !
     deter=a(1,1)*a(2,2)-a(2,1)*a(1,2)
     if(deter==0.0_rp) return
     denom  = 1.0_rp/deter
     b(1,1) = a(2,2)*denom
     b(2,2) = a(1,1)*denom
     b(2,1) =-a(2,1)*denom
     b(1,2) =-a(1,2)*denom

  else if(nsize==3) then
     !
     ! Inverse of a 3*3 matrix
     !
     t1  = a(2,2)*a(3,3) - a(3,2)*a(2,3)
     t2  =-a(2,1)*a(3,3) + a(3,1)*a(2,3)
     t3  = a(2,1)*a(3,2) - a(3,1)*a(2,2)
     deter = a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3
     if(deter==0.0_rp) return
     denom  = 1.0_rp/deter
     b(1,1) = t1*denom
     b(2,1) = t2*denom
     b(3,1) = t3*denom
     b(2,2) = ( a(1,1)*a(3,3) - a(3,1)*a(1,3))*denom
     b(3,2) = (-a(1,1)*a(3,2) + a(1,2)*a(3,1))*denom
     b(3,3) = ( a(1,1)*a(2,2) - a(2,1)*a(1,2))*denom
     b(1,2) = (-a(1,2)*a(3,3) + a(3,2)*a(1,3))*denom
     b(1,3) = ( a(1,2)*a(2,3) - a(2,2)*a(1,3))*denom
     b(2,3) = (-a(1,1)*a(2,3) + a(2,1)*a(1,3))*denom

  else if(nsize==4) then
     !
     ! Inverse of a 4*4 matrix
     !
     t1=   a(2,2)*a(3,3)*a(4,4) + a(2,3)*a(3,4)*a(4,2)&
          +a(2,4)*a(3,2)*a(4,3) - a(2,3)*a(3,2)*a(4,4)&
          -a(2,2)*a(3,4)*a(4,3) - a(2,4)*a(3,3)*a(4,2)
     t2=  -a(2,1)*a(3,3)*a(4,4) - a(2,3)*a(3,4)*a(4,1)&
          -a(2,4)*a(3,1)*a(4,3) + a(2,4)*a(3,3)*a(4,1)&
          +a(2,3)*a(3,1)*a(4,4) + a(2,1)*a(3,4)*a(4,3)
     t3=   a(2,1)*a(3,2)*a(4,4) + a(2,2)*a(3,4)*a(4,1)&
          +a(2,4)*a(3,1)*a(4,2) - a(2,4)*a(3,2)*a(4,1)&
          -a(2,2)*a(3,1)*a(4,4) - a(2,1)*a(3,4)*a(4,2)
     t4=  -a(2,1)*a(3,2)*a(4,3) - a(2,2)*a(3,3)*a(4,1)&
          -a(2,3)*a(3,1)*a(4,2) + a(2,3)*a(3,2)*a(4,1)&
          +a(2,2)*a(3,1)*a(4,3) + a(2,1)*a(3,3)*a(4,2)
     deter= a(1,1)*t1 + a(1,2)*t2 + a(1,3)*t3 + a(1,4)*t4
     if(deter==0.0_rp) return
     denom = 1.0_rp/deter
     b(1,1) = t1*denom
     b(2,1) = t2*denom
     b(3,1) = t3*denom
     b(4,1) = t4*denom
     b(1,2) =(- a(1,2)*a(3,3)*a(4,4) - a(1,3)*a(3,4)*a(4,2)&
          &   - a(1,4)*a(3,2)*a(4,3) + a(1,3)*a(3,2)*a(4,4)&
          &   + a(1,2)*a(3,4)*a(4,3) + a(1,4)*a(3,3)*a(4,2))*denom
     b(2,2) =(  a(1,1)*a(3,3)*a(4,4) + a(1,3)*a(3,4)*a(4,1)&
          &   + a(1,4)*a(3,1)*a(4,3) - a(1,4)*a(3,3)*a(4,1)&
          &   - a(1,3)*a(3,1)*a(4,4) - a(1,1)*a(3,4)*a(4,3))*denom
     b(3,2) =(- a(1,1)*a(3,2)*a(4,4) - a(1,2)*a(3,4)*a(4,1)&
          &   - a(1,4)*a(3,1)*a(4,2) + a(1,4)*a(3,2)*a(4,1)&
          &   + a(1,2)*a(3,1)*a(4,4) + a(1,1)*a(3,4)*a(4,2))*denom
     b(4,2) =(  a(1,1)*a(3,2)*a(4,3) + a(1,2)*a(3,3)*a(4,1)&
          &   + a(1,3)*a(3,1)*a(4,2) - a(1,3)*a(3,2)*a(4,1)&
          &   - a(1,2)*a(3,1)*a(4,3) - a(1,1)*a(3,3)*a(4,2))*denom
     b(1,3) =(  a(1,2)*a(2,3)*a(4,4) + a(1,3)*a(2,4)*a(4,2)&
          &   + a(1,4)*a(2,2)*a(4,3) - a(1,3)*a(2,2)*a(4,4)&
          &   - a(1,2)*a(2,4)*a(4,3) - a(1,4)*a(2,3)*a(4,2))*denom
     b(2,3) =(- a(1,1)*a(2,3)*a(4,4) - a(1,3)*a(2,4)*a(4,1)&
          &   - a(1,4)*a(2,1)*a(4,3) + a(1,4)*a(2,3)*a(4,1)&
          &   + a(1,3)*a(2,1)*a(4,4) + a(1,1)*a(2,4)*a(4,3))*denom
     b(3,3) =(  a(1,1)*a(2,2)*a(4,4) + a(1,2)*a(2,4)*a(4,1)&
          &   + a(1,4)*a(2,1)*a(4,2) - a(1,4)*a(2,2)*a(4,1)&
          &   - a(1,2)*a(2,1)*a(4,4) - a(1,1)*a(2,4)*a(4,2))*denom
     b(4,3) =(- a(1,1)*a(2,2)*a(4,3) - a(1,2)*a(2,3)*a(4,1)&
          &   - a(1,3)*a(2,1)*a(4,2) + a(1,3)*a(2,2)*a(4,1)&
          &   + a(1,2)*a(2,1)*a(4,3) + a(1,1)*a(2,3)*a(4,2))*denom
     b(1,4) =(- a(1,2)*a(2,3)*a(3,4) - a(1,3)*a(2,4)*a(3,2)&
          &   - a(1,4)*a(2,2)*a(3,3) + a(1,4)*a(2,3)*a(3,2)&
          &   + a(1,3)*a(2,2)*a(3,4) + a(1,2)*a(2,4)*a(3,3))*denom
     b(2,4) =(  a(1,1)*a(2,3)*a(3,4) + a(1,3)*a(2,4)*a(3,1)&
          &   + a(1,4)*a(2,1)*a(3,3) - a(1,4)*a(2,3)*a(3,1)&
          &   - a(1,3)*a(2,1)*a(3,4) - a(1,1)*a(2,4)*a(3,3))*denom
     b(3,4) =(- a(1,1)*a(2,2)*a(3,4) - a(1,2)*a(2,4)*a(3,1)&
          &   - a(1,4)*a(2,1)*a(3,2) + a(1,4)*a(2,2)*a(3,1)&
          &   + a(1,2)*a(2,1)*a(3,4) + a(1,1)*a(2,4)*a(3,2))*denom
     b(4,4) =(  a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1)&
          &   + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)&
          &   - a(1,2)*a(2,1)*a(3,3) - a(1,1)*a(2,3)*a(3,2))*denom
  else
     !
     ! Inverse of a nsize*nsize matrix
     !
     do isize=1,nsize
        do jsize=1,nsize
           b(isize,jsize)=a(isize,jsize)
        enddo
     enddo
     call elsest_invert(b,nsize,nsize)
  end if

end subroutine elsest_invmtx
