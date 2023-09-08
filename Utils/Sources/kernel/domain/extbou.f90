subroutine extbou(itask,pnodb,lnodb,coord,bouno,param)
  !-----------------------------------------------------------------------
  !****f* extbou/extbou
  ! NAME
  !    extbou
  ! DESCRIPTION
  !    This routines computes the exterior normal to the IB 
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime

  implicit none
  integer(ip), intent(in)  :: itask,pnodb
  integer(ip), intent(in)  :: lnodb(pnodb)
  real(rp),    intent(in)  :: coord(ndime,*)
  real(rp),    intent(out) :: bouno(ndime)
  integer(ip), optional,intent(in)  :: param
  integer(ip)              :: p1,p2,p3
  real(rp)                 :: xfact,vec(3,3),tmp(3)

  if( pnodb == 2 ) then

     if (itask == 1) then
        p1 = lnodb(1)
        p2 = lnodb(2)
     else
        p1 = lnodb(2)
        p2 = lnodb(1)
     end if

     vec(1,1) = coord(1,p2) - coord(1,p1)
     vec(2,1) = coord(2,p2) - coord(2,p1)
     bouno(1) =  vec(2,1)
     bouno(2) = -vec(1,1)

  else if( pnodb == 3 ) then

     p1 = lnodb(1)
     if (itask == 1) then
        p2 = lnodb(2)
        p3 = lnodb(3)     
     else
        p2 = lnodb(3)
        p3 = lnodb(2)     
     end if
     call nortri(p1,p2,p3,coord,vec,ndime)
     bouno(1) = 0.5_rp*vec(1,3)
     bouno(2) = 0.5_rp*vec(2,3)
     bouno(3) = 0.5_rp*vec(3,3)

  else if( pnodb == 4 ) then

     p1 = lnodb(1)
     if (itask == 1) then
        p2 = lnodb(2)
        p3 = lnodb(3)     
     else
        p2 = lnodb(3)
        p3 = lnodb(2)     
     end if
     call nortri(p1,p2,p3,coord,vec,ndime)
     bouno(1) = 0.5_rp*vec(1,3)
     bouno(2) = 0.5_rp*vec(2,3)
     bouno(3) = 0.5_rp*vec(3,3)
     if (itask == 1) then
        p2 = lnodb(3)
        p3 = lnodb(4)     
     else
        p2 = lnodb(4)
        p3 = lnodb(3)     
     end if
     if (p3 < 0) then
        print *, 'extbou -->', lnodb(1:4)
        print *, 'param -->',param
     end if
     call nortri(p1,p2,p3,coord,vec,ndime)
     bouno(1) = bouno(1)+0.5_rp*vec(1,3)
     bouno(2) = bouno(2)+0.5_rp*vec(2,3)
     bouno(3) = bouno(3)+0.5_rp*vec(3,3)
     
  else

     call runend('EXTBOU: COULD NOT COMPUTE EXTERIOR NORMAL')

  end if
  call vecuni(ndime,bouno,xfact)

end subroutine extbou
