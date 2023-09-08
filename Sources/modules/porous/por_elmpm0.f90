!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_elmpm0.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Compute elemental matrix for pressure inicialization - very similar por_solite but only lapl & grav.
!> @details Compute elemental matrix for pressure inicialization - very similar por_solite but only lapl & grav.
!> @} 
!------------------------------------------------------------------------
subroutine por_elmpm0(&
     pnode,pgaus,gpdif,gpgra,&
     gpsha,gpcar,gpvol,&
     elmat,elrhs)
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  implicit none
  integer(ip), intent(in)    :: pnode,pgaus
  real(rp),    intent(in)    :: gpdif(ndime,pgaus)     ! Note that we include anisotropic difusion 
  real(rp),    intent(in)    :: gpgra(ndime,pgaus)     ! Gravity terms  (g/Bw)(rho_o*k_ro/mu_o+rho_w*k_rw/mu_w)K(idime)*gravi_por(idime)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(out)   :: elmat(pnode,pnode),elrhs(pnode)

  integer(ip)                :: inode,jnode,idime,igaus
  real(rp)                   :: fact2,gpper(pnode),elrh0(pnode),xmuit
  !
  ! Initialization
  !
  do inode = 1,pnode
     elrhs(inode) = 0.0_rp
     do jnode = 1,pnode
        elmat(jnode,inode) = 0.0_rp
     end do
  end do

  do igaus = 1,pgaus
     !
     ! calculus of residual resid and perturbation function gpper
     !
     do inode = 1,pnode
        gpper(inode) = gpsha(inode,igaus)         
     end do
     ! 
     !   Diffusion Term & gravity term
     !
     fact2 = gpvol(igaus)
     do inode = 1, pnode
        do jnode= 1,inode-1  ! Off diagonal terms
           xmuit = 0.0_rp
           do idime= 1,ndime
              xmuit = xmuit + gpcar(idime,jnode,igaus)*gpcar(idime,inode,igaus)*gpdif(idime,igaus)
           end do
           xmuit = xmuit*fact2
           elmat(inode,jnode)= elmat(inode, jnode) + xmuit
           elmat(jnode,inode)= elmat(jnode, inode) + xmuit
        end do
        xmuit = 0.0_rp       ! diagonal terms & gravity term
        elrh0(inode) =  0.0_rp
        do idime=1, ndime
           xmuit = xmuit + gpcar(idime,inode,igaus)*gpcar(idime,inode,igaus)*gpdif(idime,igaus)
           elrh0(inode) = elrh0 (inode) + gpcar(idime,inode,igaus) * gpgra(idime,igaus) * fact2
        end do
        elmat(inode,inode)= elmat(inode,inode) + xmuit*fact2
     end do
     !
     ! Assembly of the matrix add
     !
     do inode =1, pnode
        elrhs(inode) = elrhs(inode) + elrh0(inode)
     end do
  end do

end subroutine por_elmpm0

