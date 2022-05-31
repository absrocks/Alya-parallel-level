subroutine shape2(s,t,nnode,shapf,deriv,heslo,ierro)

  !-----------------------------------------------------------------------
  !
  !    This routine evaluates shape functions and their first and
  !    second derivatives for 2-d continuos standar interpolation 
  !    elements.
  !
  !    TRIANGLES       3   6  &  10  nodes
  !    QUADRILATERALS  4   9  &  16  nodes
  !
  !-----------------------------------------------------------------------

  use      def_kintyp
  implicit none
  integer(ip), intent(in)    :: nnode
  integer(ip), intent(inout) :: ierro
  real(rp),    intent(in)    :: s,t
  real(rp),    intent(out)   :: deriv(2,nnode),shapf(nnode),heslo(3,nnode)
  integer(ip)                :: ii,jj
  real(rp)                   :: st,a1,a2,a3,ss,tt,s1,s2,s3,s4
  real(rp)                   :: t1,t2,t3,t4,s9,t9,c,a

  do ii = 1,nnode
     do jj = 1,3
        heslo(jj,ii) = 0.0_rp
     end do
  end do

  if(nnode==3) then     
     shapf(1)=1.0_rp-s-t                                
     shapf(2)=s                                  
     shapf(3)=t                                           !  3
     deriv(1,1)=-1.0_rp                                   !   
     deriv(1,2)= 1.0_rp                                   !
     deriv(1,3)= 0.0_rp                                   !
     deriv(2,1)=-1.0_rp                                   !  1       2
     deriv(2,2)= 0.0_rp
     deriv(2,3)= 1.0_rp
  else if(nnode==4) then
     st=s*t                                           
     shapf(1)=(1.0_rp-t-s+st)*0.25_rp                     !  4         3
     shapf(2)=(1.0_rp-t+s-st)*0.25_rp                     !
     shapf(3)=(1.0_rp+t+s+st)*0.25_rp                     !      
     shapf(4)=(1.0_rp+t-s-st)*0.25_rp                     !
     deriv(1,1)=(-1.0_rp+t)*0.25_rp                       !  1         2
     deriv(1,2)=(+1.0_rp-t)*0.25_rp
     deriv(1,3)=(+1.0_rp+t)*0.25_rp
     deriv(1,4)=(-1.0_rp-t)*0.25_rp
     deriv(2,1)=(-1.0_rp+s)*0.25_rp
     deriv(2,2)=(-1.0_rp-s)*0.25_rp
     deriv(2,3)=(+1.0_rp+s)*0.25_rp
     deriv(2,4)=(+1.0_rp-s)*0.25_rp
     if( ierro /= -1 ) then
        heslo(3,1)= 0.25_rp
        heslo(3,2)=-0.25_rp
        heslo(3,3)= 0.25_rp
        heslo(3,4)=-0.25_rp
     end if
  else if(nnode==6) then
     a1=1.0_rp-s-t
     a2=s                                             
     a3=t
     shapf( 1)=(2.0_rp*a1-1.0_rp)*a1                      !  3
     shapf( 2)=(2.0_rp*a2-1.0_rp)*a2                      !   
     shapf( 3)=(2.0_rp*a3-1.0_rp)*a3                      !   
     shapf( 4)=4.0_rp*a1*a2                               !  6      5
     shapf( 5)=4.0_rp*a2*a3                               !     
     shapf( 6)=4.0_rp*a1*a3                               !   
     deriv(1,1)= 1.0_rp-4.0_rp*a1                         !  1     4     2
     deriv(1,2)= 4.0_rp*a2-1.0_rp    
     deriv(1,3)= 0.0_rp           
     deriv(1,4)= 4.0_rp*(a1-a2)   
     deriv(1,5)= 4.0_rp*a3        
     deriv(1,6)=-4.0_rp*a3       
     deriv(2,1)= 1.0_rp-4.0_rp*a1    
     deriv(2,2)= 0.0_rp           
     deriv(2,3)= 4.0_rp*a3-1.0_rp    
     deriv(2,4)=-4.0_rp*a2       
     deriv(2,5)= 4.0_rp*a2        
     deriv(2,6)= 4.0_rp*(a1-a3)
     if( ierro /= -1 ) then
        heslo(1,1)= 4.0_rp
        heslo(1,2)= 4.0_rp
        heslo(1,4)=-8.0_rp
        heslo(2,1)= 4.0_rp
        heslo(2,3)= 4.0_rp
        heslo(2,6)=-8.0_rp
        heslo(3,1)= 4.0_rp
        heslo(3,4)=-4.0_rp
        heslo(3,5)= 4.0_rp
        heslo(3,6)=-4.0_rp
     end if
  else if(nnode==9) then
     ss=s*s
     st=s*t
     tt=t*t
     s1=s+1.0_rp
     t1=t+1.0_rp
     s2=s*2.0_rp
     t2=t*2.0_rp
     s9=s-1.0_rp                               
     t9=t-1.0_rp                                          !  4      7      3
     shapf( 1)=0.25_rp*s9*st*t9                           !
     shapf( 2)=0.25_rp*s1*st*t9                           !        
     shapf( 3)=0.25_rp*s1*st*t1                           !      
     shapf( 4)=0.25_rp*s9*st*t1                           !  8      9      6
     shapf( 5)=0.5_rp*(1.0_rp-ss)*t*t9                    !    
     shapf( 6)=0.5_rp*s*s1*(1.0_rp-tt)                    !   
     shapf( 7)=0.5_rp*(1.0_rp-ss)*t*t1                    !  
     shapf( 8)=0.5_rp*s*s9*(1.0_rp-tt)                    !  1      5      2
     shapf( 9)=(1.0_rp-ss)*(1.0_rp-tt)
     deriv(1,1)= 0.25_rp*t*t9*(-1.0_rp+s2)
     deriv(1,2)= 0.25_rp*(1.0_rp+s2)*t*t9
     deriv(1,3)= 0.25_rp*(1.0_rp+s2)*t*t1
     deriv(1,4)= 0.25_rp*(-1.0_rp+s2)*t*t1
     deriv(1,5)=-st*t9
     deriv(1,6)= 0.5_rp*(1.0_rp+s2)*(1.0_rp-tt)
     deriv(1,7)=-st*t1
     deriv(1,8)= 0.5_rp*(-1.0_rp+s2)*(1.0_rp-tt)
     deriv(1,9)=-s2*(1.0_rp-tt)
     deriv(2,1)= 0.25_rp*(-1.0_rp+t2)*s*s9
     deriv(2,2)= 0.25_rp*s*s1*(-1.0_rp+t2)
     deriv(2,3)= 0.25_rp*s*s1*(1.0_rp+t2)
     deriv(2,4)= 0.25_rp*s*s9*(1.0_rp+t2)
     deriv(2,5)= 0.5_rp*(1.0_rp-ss)*(-1.0_rp+t2)
     deriv(2,6)=-st*s1
     deriv(2,7)= 0.5_rp*(1.0_rp-ss)*(1.0_rp+t2)
     deriv(2,8)=-st*s9
     deriv(2,9)=-t2*(1.0_rp-ss)
     if( ierro /= -1 ) then
        heslo(1,1)= 0.5_rp*t*t9
        heslo(1,2)= 0.5_rp*t*t9
        heslo(1,3)= 0.5_rp*t*t1
        heslo(1,4)= 0.5_rp*t*t1
        heslo(1,5)=-t*t9
        heslo(1,6)= 1.0_rp-tt
        heslo(1,7)=-t*t1
        heslo(1,8)= 1.0_rp-tt
        heslo(1,9)=-2.0_rp*(1.0_rp-tt)
        heslo(2,1)= 0.5_rp*s*s9
        heslo(2,2)= 0.5_rp*s*s1
        heslo(2,3)= 0.5_rp*s*s1
        heslo(2,4)= 0.5_rp*s*s9
        heslo(2,5)= 1.0_rp-ss
        heslo(2,6)=-s*s1
        heslo(2,7)= 1.0_rp-ss
        heslo(2,8)=-s*s9
        heslo(2,9)=-2.0_rp*(1.0_rp-ss)
        heslo(3,1)= 0.25_rp*(-1.0_rp+t2)*(s9+s)
        heslo(3,2)= 0.25_rp*(-1.0_rp+t2)*(s1+s)
        heslo(3,3)= 0.25_rp*(1.0_rp+t2)*(s1+s)
        heslo(3,4)= 0.25_rp*(1.0_rp+t2)*(s9+s)
        heslo(3,5)=-s*(-1.0_rp+t2)
        heslo(3,6)=-t*s1-st
        heslo(3,7)=-s*(1.0_rp+t2)
        heslo(3,8)=-t*s9-st
        heslo(3,9)= s2*t2
     end if
  else if(nnode==10) then
     c=9.0_rp/2.0_rp
     a1=1.0_rp-s-t
     a2=2.0_rp/3.0_rp-s-t
     a3=1.0_rp/3.0_rp-s-t
     shapf( 1)=c*a1*a2*a3                                 !  3            
     shapf( 2)=c*(1.0_rp/3.0_rp-s)*(2.0_rp/3.0_rp-s)*s    !               
     shapf( 3)=c*(1.0_rp/3.0_rp-t)*(2.0_rp/3.0_rp-t)*t    !               
     shapf( 4)= 3.0_rp*c*a1*a2*s                          !  8    7  
     shapf( 5)=-3.0_rp*c*a1*(1.0_rp/3.0_rp-s)             !               
     shapf( 6)=-3.0_rp*c*(1.0_rp/3.0_rp-s)*s*t            !               
     shapf( 7)=-3.0_rp*c*s*(1.0_rp/3.0_rp-t)*t            !  9   10    6
     shapf( 8)=-3.0_rp*c*a1*(1.0_rp/3.0_rp-t)*t           !
     shapf( 9)= 3.0_rp*c*a1*a2*t                          !
     shapf(10)= 6.0_rp*c*a1*s*t                           !  1    4    5    2
     deriv(1, 1)=-c*(a1*a2+a1*a3+a2*a3)       
     deriv(1, 2)=-c*((2.0_rp/3.0_rp-s)*s&
          + (1.0_rp/3.0_rp-s)*s-(1.0_rp/3.0_rp-s)*(2.0_rp/3.0_rp-s))
     deriv(1, 3)=0.0_rp
     deriv(1, 4)= 3.0_rp*c*(a1*a2-a1*s-a2*s)
     deriv(1, 5)=-3.0_rp*c*(a1*(1.0_rp/3.0_rp-s)&
          - a1*s-(1.0_rp/3.0_rp-s)*s)
     deriv(1, 6)=-3.0_rp*c*((1.0_rp/3.0_rp-s)*t-s*t)
     deriv(1, 7)=-3.0_rp*c*((1.0_rp/3.0_rp-t)*t)
     deriv(1, 8)= 3.0_rp*c*((1.0_rp/3.0_rp-t)*t)
     deriv(1, 9)= 3.0_rp*c*(-a1*t-a2*t)
     deriv(1,10)= 6.0_rp*c*(a1*t-s*t)
     deriv(2, 1)=-c*(a1*a2+a1*a3+a2*a3)
     deriv(2, 2)= 0.0_rp
     deriv(2, 3)=-c*((2.0_rp/3.0_rp-t)*t&
          + (1.0_rp/3.0_rp-t)*t-(1.0_rp/3.0_rp-t)*(2.0_rp/3.0_rp-t))
     deriv(2, 4)= 3.0_rp*c*(-a1*s-a2*s)
     deriv(2, 5)=-3.0_rp*c*(-(1.0_rp/3.0_rp-s)*s)
     deriv(2, 6)=-3.0_rp*c*((1.0_rp/3.0_rp-s)*s)
     deriv(2, 7)=-3.0_rp*c*((1.0_rp/3.0_rp-t)*s-s*t)
     deriv(2, 8)=-3.0_rp*c*(-(1.0_rp/3.0_rp-t)*t&
          - a1*t+a1*(1.0_rp/3.0_rp-t))
     deriv(2, 9)= 3.0_rp*c*(-a1*t-a2*t+a1*a2)
     deriv(2,10)= 6.0_rp*c*(a1*s-s*t)
     if( ierro /= -1 ) then
        heslo(1, 1)= 2.0_rp*c*(a1+a2+a3) 
        heslo(1, 2)=-c*(2.0_rp-6.0_rp*s) 
        heslo(1, 3)= 0.0_rp 
        heslo(1, 4)= 3.0_rp*c*(-2.0_rp*a1-2.0_rp*a2+2.0_rp*s) 
        heslo(1, 5)=-3.0_rp*c*(-2.0_rp*a1-2.0_rp*(1.0_rp/3.0_rp-s)+2.0_rp*s) 
        heslo(1, 6)= 3.0_rp*c*2.0_rp*t 
        heslo(1, 7)= 0.0_rp 
        heslo(1, 8)= 0.0_rp  
        heslo(1, 9)= 3.0_rp*c*2.0_rp*t 
        heslo(1,10)=-6.0_rp*c*2.0_rp*t 
        heslo(2, 1)= c*(2.0_rp*a1+2.0_rp*a2+2.0_rp*a3) 
        heslo(2, 2)= 0.0_rp 
        heslo(2, 3)=-c*(2.0_rp-6.0_rp*t)
        heslo(2, 4)= 3.0_rp*c*2.0_rp*s
        heslo(2, 5)= 0.0_rp
        heslo(2, 6)= 0.0_rp
        heslo(2, 7)= 3.0_rp*c*2.0_rp*s
        heslo(2, 8)=-3.0_rp*c*(-2.0_rp*a1-2.0_rp*(1.0_rp/3.0_rp-t)+2.0_rp*t)
        heslo(2, 9)= 3.0_rp*c*(-2.0_rp*a1-2.0_rp*a2+2.0_rp*t)
        heslo(2,10)=-6.0_rp*c*2.0_rp*s
        heslo(3, 1)= 2.0_rp*c*(a1+a2+a3) 
        heslo(3, 2)= 0.0_rp  
        heslo(3, 3)= 0.0_rp 
        heslo(3, 4)= 3.0_rp*c*(-a1-a2+2.0_rp*s) 
        heslo(3, 5)=-3.0_rp*c*(-(1.0_rp/3.0_rp-s)+s) 
        heslo(3, 6)=-3.0_rp*c*(1.0_rp/3.0_rp-2.0_rp*s)
        heslo(3, 7)=-3.0_rp*c*(1.0_rp/3.0_rp-2.0_rp*t)
        heslo(3, 8)= 3.0_rp*c*(1.0_rp/3.0_rp-2.0_rp*t) 
        heslo(3, 9)= 3.0_rp*c*(-a1-a2+2.0_rp*t) 
        heslo(3,10)= 6.0_rp*c*(a1-s-t)
     end if
  else if(nnode==16) then
     a =81.0_rp/256.0_rp
     c =1.0_rp/3.0_rp
     s1=1.0_rp+s
     s2=c+s
     s3=c-s
     s4=1.0_rp-s
     t1=1.0_rp+t
     t2=c+t
     t3=c-t
     t4=1.0_rp-t
     shapf( 1) =   a*s2*s3*s4*t2*t3*t4                   ! 4    10    9    3
     shapf( 2) =   a*s1*s2*s3*t2*t3*t4                   ! 
     shapf( 3) =   a*s1*s2*s3*t1*t2*t3                   ! 
     shapf( 4) =   a*s2*s3*s4*t1*t2*t3                   ! 11   16   15    8
     shapf( 5) =-3.0_rp*a*s1*s3*s4*t2*t3*t4              !
     shapf( 6) =-3.0_rp*a*s1*s2*s4*t2*t3*t4              !
     shapf( 7) =-3.0_rp*a*s1*s2*s3*t1*t3*t4              ! 12   13   14    7
     shapf( 8) =-3.0_rp*a*s1*s2*s3*t1*t2*t4              !
     shapf( 9) =-3.0_rp*a*s1*s2*s4*t1*t2*t3              !
     shapf(10) =-3.0_rp*a*s1*s3*s4*t1*t2*t3              ! 1     5    6    2
     shapf(11) =-3.0_rp*a*s2*s3*s4*t1*t2*t4                 
     shapf(12) =-3.0_rp*a*s2*s3*s4*t1*t3*t4
     shapf(13) = 9.0_rp*a*s1*s3*s4*t1*t3*t4
     shapf(14) = 9.0_rp*a*s1*s2*s4*t1*t3*t4
     shapf(15) = 9.0_rp*a*s1*s2*s4*t1*t2*t4
     shapf(16) = 9.0_rp*a*s1*s3*s4*t1*t2*t4
     deriv(1, 1)=  a *t2*t3*t4*(-s2*s3-s2*s4+s3*s4)
     deriv(1, 2)=  a *t2*t3*t4*(-s1*s2+s1*s3+s2*s3)
     deriv(1, 3)=  a *t1*t2*t3*(-s1*s2+s1*s3+s2*s3)
     deriv(1, 4)=  a *t1*t2*t3*(-s2*s3-s2*s4+s3*s4)
     deriv(1, 5)=-3.0_rp*a*t2*t3*t4*(-s1*s3-s1*s4+s3*s4)
     deriv(1, 6)=-3.0_rp*a*t2*t3*t4*(-s1*s2+s1*s4+s2*s4)
     deriv(1, 7)=-3.0_rp*a*t1*t3*t4*(-s1*s2+s1*s3+s2*s3)
     deriv(1, 8)=-3.0_rp*a*t1*t2*t4*(-s1*s2+s1*s3+s2*s3)
     deriv(1, 9)=-3.0_rp*a*t1*t2*t3*(-s1*s2+s1*s4+s2*s4)
     deriv(1,10)=-3.0_rp*a*t1*t2*t3*(-s1*s3-s1*s4+s3*s4)
     deriv(1,11)=-3.0_rp*a*t1*t2*t4*(-s2*s3-s2*s4+s3*s4)
     deriv(1,12)=-3.0_rp*a*t1*t3*t4*(-s2*s3-s2*s4+s3*s4)
     deriv(1,13)= 9.0_rp*a*t1*t3*t4*(-s1*s3-s1*s4+s3*s4)
     deriv(1,14)= 9.0_rp*a*t1*t3*t4*(-s1*s2+s1*s4+s2*s4)
     deriv(1,15)= 9.0_rp*a*t1*t2*t4*(-s1*s2+s1*s4+s2*s4)
     deriv(1,16)= 9.0_rp*a*t1*t2*t4*(-s1*s3-s1*s4+s3*s4)
     deriv(2, 1)=  a   *s2*s3*s4*(-t2*t3-t2*t4+t3*t4)
     deriv(2, 2)=  a   *s1*s2*s3*(-t2*t3-t2*t4+t3*t4)
     deriv(2, 3)=  a   *s1*s2*s3*(-t1*t2+t1*t3+t2*t3)
     deriv(2, 4)=  a   *s2*s3*s4*(-t1*t2+t1*t3+t2*t3)
     deriv(2, 5)= -3.0_rp*a *s1*s3*s4*(-t2*t3-t2*t4+t3*t4)
     deriv(2, 6)= -3.0_rp*a *s1*s2*s4*(-t2*t3-t2*t4+t3*t4)
     deriv(2, 7)= -3.0_rp*a *s1*s2*s3*(-t1*t3-t1*t4+t3*t4)
     deriv(2, 8)= -3.0_rp*a *s1*s2*s3*(-t1*t2+t1*t4+t2*t4)
     deriv(2, 9)= -3.0_rp*a *s1*s2*s4*(-t1*t2+t1*t3+t2*t3)
     deriv(2,10)= -3.0_rp*a *s1*s3*s4*(-t1*t2+t1*t3+t2*t3)
     deriv(2,11)= -3.0_rp*a *s2*s3*s4*(-t1*t2+t1*t4+t2*t4)
     deriv(2,12)= -3.0_rp*a *s2*s3*s4*(-t1*t3-t1*t4+t3*t4)
     deriv(2,13)=  9.0_rp*a *s1*s3*s4*(-t1*t3-t1*t4+t3*t4)
     deriv(2,14)=  9.0_rp*a *s1*s2*s4*(-t1*t3-t1*t4+t3*t4)
     deriv(2,15)=  9.0_rp*a *s1*s2*s4*(-t1*t2+t1*t4+t2*t4)
     deriv(2,16)=  9.0_rp*a *s1*s3*s4*(-t1*t2+t1*t4+t2*t4)
     if( ierro /= -1 ) then
        heslo(1, 1) =&
             a *t2*t3*t4*(2.0_rp*s2-2.0_rp*s3-2.0_rp*s4)
        heslo(1, 2) =&
             a *t2*t3*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s3)
        heslo(1, 3) =&
             a *t1*t2*t3*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s3)
        heslo(1, 4) =&
             a *t1*t2*t3*(2.0_rp*s2-2.0_rp*s3-2.0_rp*s4)
        heslo(1, 5) =&
             -3.0_rp*a *t2*t3*t4*(2.0_rp*s1-2.0_rp*s3-2.0_rp*s4)
        heslo(1, 6) =&
             -3.0_rp*a *t2*t3*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s4)
        heslo(1, 7) =&
             -3.0_rp*a *t1*t3*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s3)
        heslo(1, 8) =&
             -3.0_rp*a *t1*t2*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s3)
        heslo(1, 9) =&
             -3.0_rp*a *t1*t2*t3*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s4)
        heslo(1,10) =&
             -3.0_rp*a *t1*t2*t3*(2.0_rp*s1-2.0_rp*s3-2.0_rp*s4)
        heslo(1,11) =&
             -3.0_rp*a *t1*t2*t4*(2.0_rp*s2-2.0_rp*s3-2.0_rp*s4)
        heslo(1,12) =&
             -3.0_rp*a *t1*t3*t4*(2.0_rp*s2-2.0_rp*s3-2.0_rp*s4)
        heslo(1,13) =&
             9.0_rp*a *t1*t3*t4*(2.0_rp*s1-2.0_rp*s3-2.0_rp*s4)
        heslo(1,14) =&
             9.0_rp*a *t1*t3*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s4)
        heslo(1,15) =&
             9.0_rp*a *t1*t2*t4*(-2.0_rp*s1-2.0_rp*s2+2.0_rp*s4)
        heslo(1,16) =&
             9.0_rp*a *t1*t2*t4*(2.0_rp*s1-2.0_rp*s3-2.0_rp*s4)
        heslo(2, 1) =&
             a *s2*s3*s4*(2.0_rp*t2-2.0_rp*t3-2.0_rp*t4)
        heslo(2, 2) =&
             a *s1*s2*s3*(2.0_rp*t2-2.0_rp*t3-2.0_rp*t4)
        heslo(2, 3) =&
             a *s1*s2*s3*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t3)
        heslo(2, 4) =&
             a *s2*s3*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t3)
        heslo(2, 5) =&
             -3.0_rp*a *s1*s3*s4*(2.0_rp*t2-2.0_rp*t3-2.0_rp*t4)
        heslo(2, 6) =&
             -3.0_rp*a *s1*s2*s4*(2.0_rp*t2-2.0_rp*t3-2.0_rp*t4)
        heslo(2, 7) =&
             -3.0_rp*a *s1*s2*s3*(2.0_rp*t1-2.0_rp*t3-2.0_rp*t4)
        heslo(2, 8) =&
             -3.0_rp*a *s1*s2*s3*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t4)
        heslo(2, 9) =&
             -3.0_rp*a *s1*s2*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t3)
        heslo(2,10) =&
             -3.0_rp*a *s1*s3*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t3)
        heslo(2,11) =&
             -3.0_rp*a *s2*s3*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t4)
        heslo(2,12) =&
             -3.0_rp*a *s2*s3*s4*(2.0_rp*t1-2.0_rp*t3-2.0_rp*t4)
        heslo(2,13) =&
             9.0_rp*a *s1*s3*s4*(2.0_rp*t1-2.0_rp*t3-2.0_rp*t4)
        heslo(2,14) =&
             9.0_rp*a *s1*s2*s4*(2.0_rp*t1-2.0_rp*t3-2.0_rp*t4)
        heslo(2,15) =&
             9.0_rp*a *s1*s2*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t4)
        heslo(2,16) =&
             9.0_rp*a *s1*s3*s4*(-2.0_rp*t1-2.0_rp*t2+2.0_rp*t4)
        heslo(3, 1) =&
             a*(-s2*s3-s2*s4+s3*s4)*(-t2*t3-t2*t4+t3*t4)
        heslo(3, 2) =&
             a*(-s1*s2+s1*s3+s2*s3)*(-t2*t3-t2*t4+t3*t4)       
        heslo(3, 3) =&
             a*(-s1*s2+s1*s3+s2*s3)*(-t1*t2+t1*t3+t2*t3)
        heslo(3, 4) =&
             a*(-s2*s3-s2*s4+s3*s4)*(-t1*t2+t1*t3+t2*t3)
        heslo(3, 5) =&
             -3.0_rp*a*(-s1*s3-s1*s4+s3*s4)*(-t2*t3-t2*t4+t3*t4)
        heslo(3, 6) =&
             -3.0_rp*a*(-s1*s2+s1*s4+s2*s4)*(-t2*t3-t2*t4+t3*t4)
        heslo(3, 7) =&
             -3.0_rp*a*(-s1*s2+s1*s3+s2*s3)*(-t1*t3-t1*t4+t3*t4)
        heslo(3, 8) =&
             -3.0_rp*a*(-s1*s2+s1*s3+s2*s3)*(-t1*t2+t1*t4+t2*t4)
        heslo(3, 9) =&
             -3.0_rp*a*(-s1*s2+s1*s4+s2*s4)*(-t1*t2+t1*t3+t2*t3)
        heslo(3,10) =&
             -3.0_rp*a*(-s1*s3-s1*s4+s3*s4)*(-t1*t2+t1*t3+t2*t3)
        heslo(3,11) =&
             -3.0_rp*a*(-s2*s3-s2*s4+s3*s4)*(-t1*t2+t1*t4+t2*t4)
        heslo(3,12) =&
             -3.0_rp*a*(-s2*s3-s2*s4+s3*s4)*(-t1*t3-t1*t4+t3*t4)
        heslo(3,13) =&
             9.0_rp*a*(-s1*s3-s1*s4+s3*s4)*(-t1*t3-t1*t4+t3*t4)
        heslo(3,14) =&
             9.0_rp*a*(-s1*s2+s1*s4+s2*s4)*(-t1*t3-t1*t4+t3*t4)
        heslo(3,15) =&
             9.0_rp*a*(-s1*s2+s1*s4+s2*s4)*(-t1*t2+t1*t4+t2*t4)
        heslo(3,16) =&
             9.0_rp*a*(-s1*s3-s1*s4+s3*s4)*(-t1*t2+t1*t4+t2*t4)
     end if
  else
     ierro=1
  end if
  if( ierro == -1 ) ierro = 0

end subroutine shape2
