subroutine renum2(nadj,nstart,lev,naux,nodes,ns)
!-----------------------------------------------------------------------
!
! Compute a set of posibles startings nodes
!
!-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  logical(lg) :: better
  integer(ip) :: nodes,ns
  integer(ip) :: naux(1:*),nadj(1:*),lev(1:*),nstart(1:*)
  integer(ip) :: iroot,i,nnoded,idepth,ndepth,iwidth,nwidth,lhw,lr
!
! Begin iteration
! Select initial root node arbitrarily and generate its level
! structure
!
  iroot = 1
  better = .true.
  do while (better)
     call renum4(nstart,lev,idepth,nadj,iwidth,nodes,iroot,lhw)
!
! Create a list of nodes which are at maximum distance from root
! node and store the root 
!
     ns = iroot
!
! Loop over nodes at maximum distance from root node
! Generate level structure for each node
! Set switch if a level structure of greater depth occurs
!
     better = .false.
     do i = 1,lhw
        nnoded = nstart(i)
        call renum4(naux,lev,ndepth,nadj,nwidth,nodes,nnoded,lr)
        if(ndepth>=idepth) then
           if((ndepth/=idepth).or.(nwidth<iwidth)) then
              iroot  = nnoded
              idepth = ndepth
              iwidth = nwidth
              better = .true.
           end if
        end if
     end do
  end do
  
end subroutine renum2
