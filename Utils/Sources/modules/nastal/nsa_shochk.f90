!------------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_shochk.f90
!> @date    31/08/2012
!> @author  Mariano Vazquez
!> @brief   Check if shock capturing should be applied and compute shock speed
!> @details Check if shock capturing should be applied and compute shock speed
!> @}
!------------------------------------------------------------------------
subroutine nsa_shochk(igaus,pnode,elcod,elunk,elpre,elvel,sspee)
!-----------------------------------------------------------------------

  use      def_master
  use      def_domain
  use      def_nastal

  implicit none
  integer(ip) :: igaus,pnode,inode,jnode
  real(rp)    :: &
       elcod(ndime,mnode),elunk(ndofn_nsa,mnode,ncomp_nsa),&
       elpre(mnode),elvel(ndime,mnode),cucux(ndofn_nsa,2),sspeq(ndofn_nsa),xnode,sspee
  
  
  cucux = 0.0_rp
  jnode= 0
  xnode= -10.0_rp
  do inode= 1,pnode
     if (elcod(2,inode) > 0.00001) then
        if (elcod(1,inode) > xnode) then
           xnode = elcod(1,inode)
           jnode = inode   ! jnode del nodo con x mas grande
           ! continuidad
           cucux(ndime+1,1)= elunk(3,inode,1)  ! [variable]
           cucux(ndime+1,2)= elunk(1,inode,1)  ! [flujo]
           ! energia
           cucux(ndime+2,1)= elunk(4,inode,1)
           cucux(ndime+2,2)= elvel(1,inode) *( elpre(inode) + elunk(4,inode,1))                    
           ! momento
           cucux(      1,1)= elunk(1,inode,1)
           cucux(      1,2)= elvel(1,inode) * elunk(1,inode,1) +  elpre(inode)            
        end if
     end if
  end do
  
  
  do inode= 1,pnode
     if (elcod(2,inode) > 0.00001) then
        if (inode .ne. jnode) then
           ! continuidad
           cucux(ndime+1,1)= cucux(ndime+1,1) - elunk(3,inode,1)
           cucux(ndime+1,2)= cucux(ndime+1,2) - elunk(1,inode,1)
           ! energia
           cucux(ndime+2,1)= cucux(ndime+2,1) - elunk(4,inode,1)
           cucux(ndime+2,2)= cucux(ndime+2,2) - elvel(1,inode) *( elpre(inode) + elunk(4,inode,1))                    
           ! momento
           cucux(      1,1)= cucux(      1,1) - elunk(1,inode,1)
           cucux(      1,2)= cucux(      1,2) - (elvel(1,inode) * elunk(1,inode,1) +  elpre(inode)) 
        end if
     end if
  end do

!  if ((abs(cucux(      1,1)) .gt. 5.0e-2)&
!       .and. (abs(cucux(ndime+1,1)) .gt. 5.0e-2)&
!       .and. (abs(cucux(ndime+2,1)) .gt. 5.0e-2)) then

  if ((abs(cucux(   ndime+1,1)) .gt. 5.0e-4)) then
!     sspeq(      1)= cucux(      1,2)/cucux(      1,1)
     sspeq(ndime+1)= cucux(ndime+1,2)/cucux(ndime+1,1)
!     sspeq(ndime+2)= cucux(ndime+2,2)/cucux(ndime+2,1)
  end if

  sspee = abs(sspeq(ndime+1)) 


end subroutine nsa_shochk
