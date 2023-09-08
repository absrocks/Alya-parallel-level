!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_elswar.f90
!> @date    14/10/2013
!> @author  Herbert Owen
!> @brief   Send elmat * elswa(2) to rhs
!> @details Send elmat * elswa(2) to rhs
!> @} 
!------------------------------------------------------------------------
subroutine por_elswar(pnode,ncomp_por,elrhs,elmat,elswa)
  use def_parame

  implicit none

  integer(ip), intent(in)  :: pnode,ncomp_por
  real(rp), intent(inout)  :: elrhs(pnode) 
  real(rp), intent(in)     :: elmat(pnode,pnode)
  real(rp), intent(in)     :: elswa(pnode,ncomp_por)

  integer(ip)              :: inode,jnode

  do inode = 1,pnode 
     do jnode = 1,pnode
        elrhs(inode) = elrhs(inode) - elmat(inode,jnode) * elswa(jnode,2)
     end do
  end do

end subroutine por_elswar
