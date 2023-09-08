subroutine ibm_dseseg(cooi1,cooi2,cooj1,cooj2,dista,proji,projj,dsegi,dsegj)
  !-----------------------------------------------------------------------
  ! NAME
  !    ibm_dseseg
  ! DESCRIPTION
  !    Minimun distance between a point and a segment
  !    INPUT
  !       cooi1,cooi2 : defines the first segment
  !       cooj1,cooj2 : defines the second segment
  !    OUPUT
  !       dista: the minimum distance
  !       proji: projection on the first segmnent
  !       projj: projection on the second segmnent
  ! USED BY
  !    ibm_dsepar
  !----------------------------------------------------------------------- 
  use def_kintyp, only           :  ip,rp 
  use def_domain, only           :  ndime
  use def_master, only           :  zeror
  implicit none
  real(rp),      intent(in)      :: cooi1(ndime),cooi2(ndime),cooj1(ndime),cooj2(ndime)
  real(rp),      intent(out)     :: dista,proji(ndime),projj(ndime)
  integer(ip)                    :: idime  
  real(rp)                       :: u(ndime),v(ndime),w(ndime),a,b,c,d,e
  real(rp)                       :: sN,tN,sD,tD,denom,dsegi,dsegj

  ! Define vectors  
  u(1)     =  cooi2(1)     - cooi1(1) 
  u(2)     =  cooi2(2)     - cooi1(2) 
  u(ndime) =  cooi2(ndime) - cooi1(ndime) 

  v(1)     =  cooj2(1)     - cooj1(1)
  v(2)     =  cooj2(2)     - cooj1(2)
  v(ndime) =  cooj2(ndime) - cooj1(ndime)

  w(1)     =  cooi1(1)     - cooj1(1)
  w(2)     =  cooi1(2)     - cooj1(2)
  w(ndime) =  cooi1(ndime) - cooj1(ndime)
  
  ! Define some useful quantities
  a = 0.0_rp
  b = 0.0_rp
  c = 0.0_rp
  d = 0.0_rp
  e = 0.0_rp
  do idime=1,ndime
     a = a + u(idime)*u(idime)
     b = b + u(idime)*v(idime)
     c = c + v(idime)*v(idime)
     d = d + u(idime)*w(idime)
     e = e + v(idime)*w(idime)
  end do

  denom = a*c - b*b       
  sD = denom
  tD = denom

  ! The lines are parallel
  if (denom < zeror) then    
     sN = 0.0_rp
     sD = 1.0_rp
     tN = e
     tD = c
  else
     sN = b*e - c*d
     tN = a*e - b*d

     if (sN < 0.0_rp) then
        sN = 0.0_rp
        tN = e
        tD = c
     elseif (sN > sD) then
        sN = sD
        tN = e + b
        tD = c
     end if
  end if

  if (tN < 0.0_rp) then
     tN = 0.0_rp
     if ( -d < 0.0_rp ) then
        sN = 0.0_rp
     elseif ( -d > a ) then
        sN = sD
     else
        sN = -d
        sD =  a
     end if
  elseif (tN > tD) then
     tN = tD
     if ( -d+b < 0.0_rp ) then
        sN = 0.0_rp
     elseif ( -d+b > a ) then
        sN = sD
     else
        sN = -d+b
        sD =  a
     end if
  end if
 
  dsegi = sN/sD
  dsegj = tN/tD
    
  proji(1)     = cooi1(1)     + dsegi*u(1) 
  proji(2)     = cooi1(2)     + dsegi*u(2) 
  proji(ndime) = cooi1(ndime) + dsegi*u(ndime) 

  projj(1)     = cooj1(1)     + dsegj*v(1) 
  projj(2)     = cooj1(2)     + dsegj*v(2) 
  projj(ndime) = cooj1(ndime) + dsegj*v(ndime) 


  dista = 0.0_rp
  do idime=1,ndime     
     dista = dista + (projj(idime) - proji(idime))*(projj(idime) - proji(idime))     
  end do

end subroutine ibm_dseseg
