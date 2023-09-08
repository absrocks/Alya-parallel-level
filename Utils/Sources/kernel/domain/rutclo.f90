subroutine rutclo(ndime,ngaus,posgp,weigp,ierro)

  !-----------------------------------------------------------------------
  ! 
  !     This routine sets up the integration constants of closed rules for
  !     triangles and tetrahedra 
  ! 
  !             NDIME = 2             NDIME = 3
  ! 
  !          NGAUS  EXACT POL.     NGAUS  EXACT POL. 
  !          -----  ----------     -----  ----------
  !            3       p1            4       p1
  !            4       p2            5       p2
  !            6       p1           10       p2
  !            7       p3           11       p3
  !           10       p1           15       p3
  !                                 20       p2 & x^4,...      
  ! 
  !-----------------------------------------------------------------------

  use      def_kintyp
  implicit none
  integer(ip), intent(out) :: ierro
  integer(ip), intent(in)  :: ndime,ngaus
  real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
  real(rp)                 :: w1,w2,w3
  !
  ! Line integral
  !
  ierro=0

  if(ndime==2) then
     !
     ! Area integral (triangles)
     !
     if(ngaus==3) then
        posgp(1,1)= 0.0_rp
        posgp(2,1)= 0.0_rp
        posgp(1,2)= 1.0_rp
        posgp(2,2)= 0.0_rp
        posgp(1,3)= 0.0_rp
        posgp(2,3)= 1.0_rp
        weigp(  1)= 1.0_rp/6.0_rp
        weigp(  2)= 1.0_rp/6.0_rp
        weigp(  3)= 1.0_rp/6.0_rp
     else if(ngaus==4) then
        posgp(1,1)= 0.0_rp
        posgp(2,1)= 0.0_rp
        posgp(1,2)= 1.0_rp
        posgp(2,2)= 0.0_rp
        posgp(1,3)= 0.0_rp
        posgp(2,3)= 1.0_rp
        posgp(1,4)= 1.0_rp/3.0_rp
        posgp(2,4)= 1.0_rp/3.0_rp
        weigp(  1)= 1.0_rp/24.0_rp
        weigp(  2)= 1.0_rp/24.0_rp
        weigp(  3)= 1.0_rp/24.0_rp
        weigp(  4)= 3.0_rp/8.0_rp
     else if(ngaus==6) then
        posgp(1,1)= 0.0_rp
        posgp(2,1)= 0.0_rp
        posgp(1,2)= 1.0_rp
        posgp(2,2)= 0.0_rp
        posgp(1,3)= 0.0_rp
        posgp(2,3)= 1.0_rp
        posgp(1,4)= 0.5_rp
        posgp(2,4)= 0.0_rp
        posgp(1,5)= 0.5_rp
        posgp(2,5)= 0.5_rp
        posgp(1,6)= 0.0_rp
        posgp(2,6)= 0.5_rp
        weigp(  1)= 1.0_rp/24.0_rp
        weigp(  2)= 1.0_rp/24.0_rp
        weigp(  3)= 1.0_rp/24.0_rp
        weigp(  4)= 1.0_rp/8.0_rp
        weigp(  5)= 1.0_rp/8.0_rp
        weigp(  6)= 1.0_rp/8.0_rp
     else if(ngaus==7) then
        posgp(1,1)= 0.0_rp
        posgp(2,1)= 0.0_rp
        posgp(1,2)= 1.0_rp
        posgp(2,2)= 0.0_rp
        posgp(1,3)= 0.0_rp
        posgp(2,3)= 1.0_rp
        posgp(1,4)= 0.5_rp
        posgp(2,4)= 0.0_rp
        posgp(1,5)= 0.5_rp
        posgp(2,5)= 0.5_rp
        posgp(1,6)= 0.0_rp
        posgp(2,6)= 0.5_rp
        posgp(1,7)= 1.0_rp/3.0_rp
        posgp(2,7)= 1.0_rp/3.0_rp
        weigp(  1)= 1.0_rp/40.0_rp
        weigp(  2)= 1.0_rp/40.0_rp
        weigp(  3)= 1.0_rp/40.0_rp
        weigp(  4)= 1.0_rp/15.0_rp
        weigp(  5)= 1.0_rp/15.0_rp
        weigp(  6)= 1.0_rp/15.0_rp
        weigp(  7)= 9.0_rp/40.0_rp
     else if(ngaus==10) then
        posgp(1, 1)= 0.0_rp 
        posgp(2, 1)= 0.0_rp
        posgp(1, 2)= 1.0_rp
        posgp(2, 2)= 0.0_rp
        posgp(1, 3)= 0.0_rp
        posgp(2, 3)= 1.0_rp
        posgp(1, 4)= 1.0_rp/3.0_rp
        posgp(2, 4)= 0.0_rp
        posgp(1, 5)= 2.0_rp/3.0_rp
        posgp(2, 5)= 0.0_rp
        posgp(1, 6)= 2.0_rp/3.0_rp 
        posgp(2, 6)= 1.0_rp/3.0_rp
        posgp(1, 7)= 1.0_rp/3.0_rp
        posgp(2, 7)= 2.0_rp/3.0_rp
        posgp(1, 8)= 0.0_rp
        posgp(2, 8)= 2.0_rp/3.0_rp
        posgp(1, 9)= 0.0_rp
        posgp(2, 9)= 1.0_rp/3.0_rp
        posgp(1,10)= 1.0_rp/3.0_rp
        posgp(2,10)= 1.0_rp/3.0_rp
        weigp(   1)= 1.0_rp/54.0_rp
        weigp(   2)= 1.0_rp/54.0_rp
        weigp(   3)= 1.0_rp/54.0_rp
        weigp(   4)= 1.0_rp/18.0_rp
        weigp(   5)= 1.0_rp/18.0_rp
        weigp(   6)= 1.0_rp/18.0_rp
        weigp(   7)= 1.0_rp/18.0_rp
        weigp(   8)= 1.0_rp/18.0_rp 
        weigp(   9)= 1.0_rp/18.0_rp 
        weigp(  10)= 1.0_rp/9.0_rp
     else
        ierro=1
     end if

  else if(ndime==3) then
     !
     ! Volume integral ( tetrahedra )
     !
     if(ngaus==4) then
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
        weigp(  1)= 1.0_rp/24.0_rp
        weigp(  2)= 1.0_rp/24.0_rp
        weigp(  3)= 1.0_rp/24.0_rp
        weigp(  4)= 1.0_rp/24.0_rp
     else if(ngaus==5) then
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
        posgp(1,5)= 1.0_rp/4.0_rp
        posgp(2,5)= 1.0_rp/4.0_rp
        posgp(3,5)= 1.0_rp/4.0_rp
        weigp(  1)= 1.0_rp/120.0_rp
        weigp(  2)= 1.0_rp/120.0_rp
        weigp(  3)= 1.0_rp/120.0_rp
        weigp(  4)= 1.0_rp/120.0_rp
        weigp(  5)= 2.0_rp/15.0_rp
     else if(ngaus==10) then
        posgp(1, 1)= 0.0_rp 
        posgp(2, 1)= 0.0_rp
        posgp(3, 1)= 0.0_rp
        posgp(1, 2)= 1.0_rp
        posgp(2, 2)= 0.0_rp
        posgp(3, 2)= 0.0_rp
        posgp(1, 3)= 0.0_rp
        posgp(2, 3)= 1.0_rp
        posgp(3, 3)= 0.0_rp
        posgp(1, 4)= 0.0_rp
        posgp(2, 4)= 0.0_rp
        posgp(3, 4)= 1.0_rp
        posgp(1, 5)= 0.5_rp
        posgp(2, 5)= 0.0_rp
        posgp(3, 5)= 0.0_rp
        posgp(1, 6)= 0.5_rp
        posgp(2, 6)= 0.5_rp
        posgp(3, 6)= 0.0_rp
        posgp(1, 7)= 0.0_rp
        posgp(2, 7)= 0.5_rp
        posgp(3, 7)= 0.0_rp
        posgp(1, 8)= 0.5_rp
        posgp(2, 8)= 0.0_rp
        posgp(3, 8)= 0.5_rp
        posgp(1, 9)= 0.0_rp
        posgp(2, 9)= 0.5_rp
        posgp(3, 9)= 0.5_rp
        posgp(1,10)= 0.0_rp
        posgp(2,10)= 0.0_rp
        posgp(3,10)= 0.5_rp
        weigp(   1)=-1.0_rp/120.0_rp
        weigp(   2)=-1.0_rp/120.0_rp 
        weigp(   3)=-1.0_rp/120.0_rp 
        weigp(   4)=-1.0_rp/120.0_rp 
        weigp(   5)= 1.0_rp/30.0_rp 
        weigp(   6)= 1.0_rp/30.0_rp 
        weigp(   7)= 1.0_rp/30.0_rp 
        weigp(   8)= 1.0_rp/30.0_rp 
        weigp(   9)= 1.0_rp/30.0_rp 
        weigp(  10)= 1.0_rp/30.0_rp 
     else if(ngaus==11) then
        posgp(1, 1)= 0.0_rp 
        posgp(2, 1)= 0.0_rp
        posgp(3, 1)= 0.0_rp
        posgp(1, 2)= 1.0_rp
        posgp(2, 2)= 0.0_rp
        posgp(3, 2)= 0.0_rp
        posgp(1, 3)= 0.0_rp
        posgp(2, 3)= 1.0_rp
        posgp(3, 3)= 0.0_rp
        posgp(1, 4)= 0.0_rp
        posgp(2, 4)= 0.0_rp
        posgp(3, 4)= 1.0_rp
        posgp(1, 5)= 0.5_rp
        posgp(2, 5)= 0.0_rp
        posgp(3, 5)= 0.0_rp
        posgp(1, 6)= 0.5_rp
        posgp(2, 6)= 0.5_rp
        posgp(3, 6)= 0.0_rp
        posgp(1, 7)= 0.0_rp
        posgp(2, 7)= 0.5_rp
        posgp(3, 7)= 0.0_rp
        posgp(1, 8)= 0.5_rp
        posgp(2, 8)= 0.0_rp
        posgp(3, 8)= 0.5_rp
        posgp(1, 9)= 0.0_rp
        posgp(2, 9)= 0.5_rp
        posgp(3, 9)= 0.5_rp
        posgp(1,10)= 0.0_rp
        posgp(2,10)= 0.0_rp
        posgp(3,10)= 0.5_rp
        posgp(1,11)= 1.0_rp/4.0_rp
        posgp(2,11)= 1.0_rp/4.0_rp
        posgp(3,11)= 1.0_rp/4.0_rp
        weigp(   1)= 1.0_rp/360.0_rp
        weigp(   2)= 1.0_rp/360.0_rp
        weigp(   3)= 1.0_rp/360.0_rp
        weigp(   4)= 1.0_rp/360.0_rp
        weigp(   5)= 1.0_rp/90.0_rp
        weigp(   6)= 1.0_rp/90.0_rp
        weigp(   7)= 1.0_rp/90.0_rp
        weigp(   8)= 1.0_rp/90.0_rp
        weigp(   9)= 1.0_rp/90.0_rp
        weigp(  10)= 1.0_rp/90.0_rp
        weigp(  11)= 4.0_rp/45.0_rp
     else if(ngaus==15) then
        posgp(1, 1)= 0.0_rp 
        posgp(2, 1)= 0.0_rp
        posgp(3, 1)= 0.0_rp
        posgp(1, 2)= 1.0_rp
        posgp(2, 2)= 0.0_rp
        posgp(3, 2)= 0.0_rp
        posgp(1, 3)= 0.0_rp
        posgp(2, 3)= 1.0_rp
        posgp(3, 3)= 0.0_rp
        posgp(1, 4)= 0.0_rp
        posgp(2, 4)= 0.0_rp
        posgp(3, 4)= 1.0_rp
        posgp(1, 5)= 0.5_rp
        posgp(2, 5)= 0.0_rp
        posgp(3, 5)= 0.0_rp
        posgp(1, 6)= 0.5_rp
        posgp(2, 6)= 0.5_rp
        posgp(3, 6)= 0.0_rp
        posgp(1, 7)= 0.0_rp
        posgp(2, 7)= 0.5_rp
        posgp(3, 7)= 0.0_rp
        posgp(1, 8)= 0.0_rp
        posgp(2, 8)= 0.0_rp
        posgp(3, 8)= 0.5_rp
        posgp(1, 9)= 0.5_rp
        posgp(2, 9)= 0.0_rp
        posgp(3, 9)= 0.5_rp
        posgp(1,10)= 0.0_rp
        posgp(2,10)= 0.5_rp
        posgp(3,10)= 0.5_rp
        posgp(1,11)= 1.0_rp/3.0_rp
        posgp(2,11)= 1.0_rp/3.0_rp
        posgp(3,11)= 0.0_rp
        posgp(1,12)= 1.0_rp/3.0_rp
        posgp(2,12)= 0.0_rp
        posgp(3,12)= 1.0_rp/3.0_rp
        posgp(1,13)= 1.0_rp/3.0_rp
        posgp(2,13)= 1.0_rp/3.0_rp
        posgp(3,13)= 1.0_rp/3.0_rp
        posgp(1,14)= 0.0_rp
        posgp(2,14)= 1.0_rp/3.0_rp
        posgp(3,14)= 1.0_rp/3.0_rp
        posgp(1,15)= 1.0_rp/4.0_rp
        posgp(2,15)= 1.0_rp/4.0_rp
        posgp(3,15)= 1.0_rp/4.0_rp
        weigp(   1)= 17.0_rp/5040.0_rp
        weigp(   2)= 17.0_rp/5040.0_rp
        weigp(   3)= 17.0_rp/5040.0_rp
        weigp(   4)= 17.0_rp/5040.0_rp
        weigp(   5)=  2.0_rp/315.0_rp
        weigp(   6)=  2.0_rp/315.0_rp
        weigp(   7)=  2.0_rp/315.0_rp
        weigp(   8)=  2.0_rp/315.0_rp
        weigp(   9)=  2.0_rp/315.0_rp
        weigp(  10)=  2.0_rp/315.0_rp
        weigp(  11)=  9.0_rp/840.0_rp
        weigp(  12)=  9.0_rp/840.0_rp
        weigp(  13)=  9.0_rp/840.0_rp
        weigp(  14)=  9.0_rp/840.0_rp
        weigp(  15)= 16.0_rp/315.0_rp
     else if(ngaus==20) then                           ! Integrates P2
        posgp(1, 1)= 0.0_rp                                ! and quartic mono-
        posgp(2, 1)= 0.0_rp                                ! mials to avoid zero 
        posgp(3, 1)= 0.0_rp                                ! weights at the 
        posgp(1, 2)= 1.0_rp                                ! edges
        posgp(2, 2)= 0.0_rp
        posgp(3, 2)= 0.0_rp
        posgp(1, 3)= 0.0_rp
        posgp(2, 3)= 1.0_rp
        posgp(3, 3)= 0.0_rp
        posgp(1, 4)= 0.0_rp
        posgp(2, 4)= 0.0_rp
        posgp(3, 4)= 1.0_rp
        posgp(1, 5)= 1.0_rp/3.0_rp
        posgp(2, 5)= 0.0_rp
        posgp(3, 5)= 0.0_rp
        posgp(1, 6)= 2.0_rp/3.0_rp
        posgp(2, 6)= 0.0_rp
        posgp(3, 6)= 0.0_rp
        posgp(1, 7)= 2.0_rp/3.0_rp
        posgp(2, 7)= 1.0_rp/3.0_rp
        posgp(3, 7)= 0.0_rp
        posgp(1, 8)= 1.0_rp/3.0_rp
        posgp(2, 8)= 2.0_rp/3.0_rp
        posgp(3, 8)= 0.0_rp
        posgp(1, 9)= 0.0_rp
        posgp(2, 9)= 2.0_rp/3.0_rp
        posgp(3, 9)= 0.0_rp
        posgp(1,10)= 0.0_rp
        posgp(2,10)= 1.0_rp/3.0_rp
        posgp(3,10)= 0.0_rp
        posgp(1,11)= 0.0_rp
        posgp(2,11)= 0.0_rp
        posgp(3,11)= 1.0_rp/3.0_rp
        posgp(1,12)= 2.0_rp/3.0_rp
        posgp(2,12)= 0.0_rp
        posgp(3,12)= 1.0_rp/3.0_rp
        posgp(1,13)= 0.0_rp
        posgp(2,13)= 2.0_rp/3.0_rp
        posgp(3,13)= 1.0_rp/3.0_rp
        posgp(1,14)= 0.0_rp
        posgp(2,14)= 0.0_rp
        posgp(3,14)= 2.0_rp/3.0_rp
        posgp(1,15)= 1.0_rp/3.0_rp
        posgp(2,15)= 0.0_rp
        posgp(3,15)= 2.0_rp/3.0_rp
        posgp(1,16)= 0.0_rp
        posgp(2,16)= 1.0_rp/3.0_rp
        posgp(3,16)= 2.0_rp/3.0_rp
        posgp(1,17)= 1.0_rp/3.0_rp
        posgp(2,17)= 1.0_rp/3.0_rp
        posgp(3,17)= 0.0_rp
        posgp(1,18)= 1.0_rp/3.0_rp
        posgp(2,18)= 0.0_rp
        posgp(3,18)= 1.0_rp/3.0_rp
        posgp(1,19)= 1.0_rp/3.0_rp
        posgp(2,19)= 1.0_rp/3.0_rp
        posgp(3,19)= 1.0_rp/3.0_rp
        posgp(1,20)= 0.0_rp
        posgp(2,20)= 1.0_rp/3.0_rp
        posgp(3,20)= 1.0_rp/3.0_rp
        w1 = 383.0_rp/2400.0_rp
        w2 = 1.0_rp/240.0_rp - w1
        w3 = 3.0_rp/80.0_rp  - w2
        weigp(   1)= w1
        weigp(   2)= w1
        weigp(   3)= w1
        weigp(   4)= w1
        weigp(   5)= w2
        weigp(   6)= w2
        weigp(   7)= w2
        weigp(   8)= w2
        weigp(   9)= w2
        weigp(  10)= w2
        weigp(  11)= w2
        weigp(  12)= w2
        weigp(  13)= w2
        weigp(  14)= w2
        weigp(  15)= w2
        weigp(  16)= w2
        weigp(  17)= w3
        weigp(  18)= w3
        weigp(  19)= w3
        weigp(  20)= w3
     else
        ierro=1
     end if
  end if
  !
  ! Errors
  !
  !if(ierro==1) call runend('RUTCLO: NOT AVAILABLE INTEGRATION RULE')

end subroutine rutclo
