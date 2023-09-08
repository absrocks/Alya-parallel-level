subroutine ruyope(ndime,ngaus,posgp,weigp,ierro)
  !-----------------------------------------------------------------------
  !****f* domain/ruyope
  ! NAME 
  !    ruyope
  ! DESCRIPTION
  !    This routine sets up the integration constants of open rules for
  !    the PYRA_5 element. From the following paper:
  !
  !    C.A. Felippa, A compendium of FEM integration formulas for symbolic 
  !    work, Engineering Computations 21(8), pp. 867-890 (2004). 
  !   
  ! USES
  ! USED BY
  !    rulepw
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in)    :: ndime,ngaus
  integer(ip), intent(inout) :: ierro
  real(rp),    intent(out)   :: posgp(ndime,ngaus)
  real(rp),    intent(out)   :: weigp(ngaus)
  integer(ip)                :: ii,m
  real(rp)                   :: g1,g2,g3,g4,g5,gg2(2)
  real(rp)                   :: j,k,wg9(3)
  real(rp)                   :: w1,w2,w3,w4,ww2(2)
  real(rp)                   :: jk(2,20),jk4(2,20),jk9(2,9)

  ierro    =  0

  jk(1,1)  = -1.0_rp 
  jk(2,1)  = -1.0_rp
  jk(1,2)  =  1.0_rp
  jk(2,2)  = -1.0_rp
  jk(1,3)  =  1.0_rp
  jk(2,3)  =  1.0_rp
  jk(1,4)  = -1.0_rp
  jk(2,4)  =  1.0_rp

  jk4(1,1) = -1.0_rp 
  jk4(2,1) =  0.0_rp
  jk4(1,2) =  1.0_rp
  jk4(2,2) =  0.0_rp
  jk4(1,3) =  0.0_rp
  jk4(2,3) = -1.0_rp
  jk4(1,4) =  0.0_rp
  jk4(2,4) =  1.0_rp

  jk9(1,1) = -1.0_rp 
  jk9(2,1) = -1.0_rp
  jk9(1,2) =  0.0_rp
  jk9(2,2) = -1.0_rp
  jk9(1,3) =  1.0_rp
  jk9(2,3) = -1.0_rp
  jk9(1,4) = -1.0_rp
  jk9(2,4) =  0.0_rp
  jk9(1,5) =  0.0_rp 
  jk9(2,5) =  0.0_rp
  jk9(1,6) =  1.0_rp
  jk9(2,6) =  0.0_rp
  jk9(1,7) = -1.0_rp
  jk9(2,7) =  1.0_rp
  jk9(1,8) =  0.0_rp
  jk9(2,8) =  1.0_rp
  jk9(1,9) =  1.0_rp
  jk9(2,9) =  1.0_rp

  wg9(1)   = 64.0_rp/81.0_rp
  wg9(2)   = 40.0_rp/81.0_rp
  wg9(3)   = 25.0_rp/81.0_rp

  if(       ngaus == 1 ) then
     posgp(1,1) = 0.0_rp
     posgp(2,1) = 0.0_rp
     posgp(3,1) =-0.5_rp
     weigp(  1) = 128.0_rp/27.0_rp

  else if( ngaus == 5 ) then
     g1 = 8.0_rp*sqrt(2.0_rp/15.0_rp)/5.0_rp
     do ii = 1,4
        j           = jk(1,ii)
        k           = jk(2,ii)
        posgp(1,ii) = j*g1
        posgp(2,ii) = k*g1
        posgp(3,ii) =-2.0_rp/3.0_rp
        weigp(  ii) = 81.0_rp/100.0_rp           
     end do
     posgp(1,5) = 0.0_rp
     posgp(2,5) = 0.0_rp
     posgp(3,5) = 2.0_rp/5.0_rp
     weigp(  5) = 125.0_rp/27.0_rp  
      
  else if( ngaus == 6 ) then
     g1     = sqrt(21.0_rp/35.0_rp)
     gg2(1) = 1.0_rp/6.0_rp
     gg2(2) = 1.0_rp/2.0_rp
     ww2(1) = 576.0_rp/625.0_rp
     ww2(2) = 64.0_rp/15.0_rp
     do ii=1,4
        j=jk(1,ii)
        k=jk(2,ii)
        posgp(1,ii) = j*g1
        posgp(2,ii) = k*g1
        posgp(3,ii) =-2.0_rp/3.0_rp
        weigp(  ii) = 504.0_rp/625.0_rp                    
     end do
     do ii = 5,6
        j           = jk(1,ii)
        k           = jk(2,ii)
        posgp(1,ii) = 0.0_rp
        posgp(2,ii) = 0.0_rp
        posgp(3,ii) = gg2(ii-4)
        weigp(  ii) = ww2(ii-4)                    
     end do

  else if(ngaus==8) then
     g1 = sqrt(1.0_rp/3.0_rp)
     g2 = (2.0_rp*sqrt(10.0_rp)-5.0_rp)/15.0_rp
     g3 = -2.0_rp/3.0_rp-g2
     w1 = 5.0_rp*(68.0_rp+5.0_rp*sqrt(10.0_rp))/432.0_rp
     w2 = 85.0_rp/54.0_rp-w1
     do ii=1,4
        j           = jk(1,ii)
        k           = jk(2,ii)
        posgp(1,ii) = j*g1
        posgp(2,ii) = k*g1
        posgp(3,ii) = g2
        weigp(  ii) = w1                     
     end do
     do ii=5,8
        j=jk(1,ii-4)
        k=jk(2,ii-4)
        posgp(1,ii) = j*g1
        posgp(2,ii) = k*g1
        posgp(3,ii) = g3
        weigp(  ii) = w2                     
     end do

  else if(ngaus==9) then
     g1 = 8.0_rp*sqrt( (573.0_rp+5.0_rp*sqrt(2865.0_rp))&
          &           /(109825.0_rp+969.0_rp*sqrt(2865.0_rp)) )
     g2 = sqrt( 2.0_rp*(8025.0_rp+sqrt(2865.0_rp))/35.0_rp )/37.0_rp
     g3 = -(87.0_rp+sqrt(2865.0_rp))/168.0_rp
     g4 = (-87.0_rp+sqrt(2865.0_rp))/168.0_rp
     w1 = 7.0_rp*(11472415.0_rp-70057.0_rp*sqrt(2865.0_rp))/130739500.0_rp
     w2 = 84091.0_rp/68450.0_rp-w1
     do ii=1,4
        j=jk(1,ii)
        k=jk(2,ii)
        posgp(1,ii) = j*g1
        posgp(2,ii) = k*g1
        posgp(3,ii) = g3
        weigp(  ii) = w1                     
     end do
     do ii=5,8
        j=jk(1,ii-4)
        k=jk(2,ii-4)
        posgp(1,ii) = j*g2
        posgp(2,ii) = k*g2
        posgp(3,ii) = g4
        weigp(  ii) = w2                     
     end do
     posgp(1,9) = 0.0_rp
     posgp(2,9) = 0.0_rp
     posgp(3,9) = 2.0_rp/3.0_rp
     weigp(  9) = 18.0_rp/5.0_rp

  else if(ngaus==13) then
     g1 = 7.0_rp*sqrt(35.0_rp/59.0_rp)/8.0_rp
     g2 = 224.0_rp*sqrt(336633710.0_rp/33088740423.0_rp)/37.0_rp
     g3 = sqrt(37043.0_rp/35.0_rp)/56.0_rp
     g4 = -127.0_rp/153.0_rp
     g5 = 1490761.0_rp/2842826.0_rp
     w1 = 170569.0_rp/331200_rp
     w2 = 276710106577408.0_rp/1075923777052725.0_rp
     w3 = 12827693806929.0_rp/30577384040000.0_rp
     w4 = 10663383340655070643544192.0_rp/4310170528879365193704375.0_rp
     do ii=1,4
        j=jk(1,ii)
        k=jk(2,ii)
        posgp(1,ii) = j*g1
        posgp(2,ii) = k*g1
        posgp(3,ii) = -1.0_rp/7.0_rp
        weigp(  ii) = w1         
     end do
     do ii=5,8
        j=jk4(1,ii-4)
        k=jk4(2,ii-4)
        posgp(1,ii) = j*g2
        posgp(2,ii) = k*g2
        posgp(3,ii) = -9.0_rp/28.0_rp
        weigp(  ii) = w2  
     end do
     do ii=9,12
        j=jk(1,ii-8)
        k=jk(2,ii-8)
        posgp(1,ii) = j*g3
        posgp(2,ii) = k*g3
        posgp(3,ii) = g4
        weigp(  ii) = w3  
     end do
     posgp(1,13) = 0.0_rp
     posgp(2,13) = 0.0_rp
     posgp(3,13) = g5
     weigp(  13) = w4         
        
  else if(ngaus==18) then
     g1 = sqrt(3.0_rp/5.0_rp)
     g2 = 1.0_rp-2.0_rp*(10.0_rp-sqrt(10.0_rp))/15.0_rp
     g3 = -2.0_rp/3.0_rp-g2
     w1 = 5.0_rp*(68.0_rp+5.0_rp*sqrt(10.0_rp))/432.0_rp
     w2 = 85.0_rp/54.0_rp-w1
     do ii=1,9
        j=jk9(1,ii)
        k=jk9(2,ii)
        m=abs(int(j))+abs(int(k))
        posgp(1,ii) = j*g1
        posgp(2,ii) = k*g1
        posgp(3,ii) = g2
        weigp(  ii) = w1*wg9(m+1)  
     end do
     do ii=10,18
        j=jk9(1,ii-9)
        k=jk9(2,ii-9)
        m=abs(int(j))+abs(int(k))
        posgp(1,ii) = j*g1
        posgp(2,ii) = k*g1
        posgp(3,ii) = g3
        weigp(  ii) = w2*wg9(m+1)  
     end do
  else       
     ierro=1
  end if
  !      
  ! Errors
  !
  !if(ierro==1) call runend('RUYOPE: NOT AVAILABLE QUADRATURE FOR PYRAMID ELEMENT')

end subroutine ruyope
