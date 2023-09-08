subroutine renum4(ndeg,lev,lsd,nadj,mlw,nodes,nroot,nloc)
!-----------------------------------------------------------------------
!
! Compute level structure rooted at nroot
!
!-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  logical(lg) :: back
  integer(ip) :: lsd,mlw,nodes,nroot,nloc
  integer(ip) :: lev(1:*),ndeg(1:*),nadj(1:*)
  integer(ip) :: i,kount,nlocn,il,nad,jp,nhalf,nn,iq
!
! Initialization
!
  do i = 1,nodes
     lev(i) = 0
  end do
  lev(nroot) = 1    
  ndeg(1) = nroot
  back = .false.
  lsd = 1
  kount = 1
  nloc = 1
  nlocn = 0
  mlw = 1
!
! Assign levels
!
  do while(kount<nodes)
     do il = 1,nloc
        if(.not.back) iq = ndeg(il)
        if(back) iq = ndeg(nodes+1-il)
        nad = nadj(iq)
        do while (nad>0)
           jp = nadj(nad)
           if(lev(jp)==0) then
              nlocn = nlocn + 1
              lev(jp) = lsd + 1
              kount = kount + 1
              if(back) ndeg(nlocn) = jp
              if(.not.back) ndeg(nodes+1-nlocn) = jp
           end if
           nad = nadj(nad+1)
        end do
     end do
     nloc = nlocn
     nlocn = 0
     if(nloc>mlw) mlw = nloc
     lsd = lsd + 1
     back = .not.back
  end do
  
  if(back) then
     nhalf = nodes/2
     do i = 1,nhalf
        nn = ndeg(i)
        ndeg(i) = ndeg(nodes+1-i)
        ndeg(nodes+1-i) = nn
     end do
  end if
  
end subroutine renum4
