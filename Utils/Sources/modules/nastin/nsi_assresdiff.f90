subroutine nsi_assresdiff(&
                         itask,pnode,pevat,lnods,elrbu,elrbp,sens_mesh,idime_dof,ipoin_dof,elvel,elpre)
  !-----------------------------------------------------------------------
  !****f* mathru/nsi_assresdiff
  ! NAME 
  !    nsi_assresdiff
  ! DESCRIPTION
  !    Assembly an elemental matrix ELMAT in global matrix AMATR
  ! USES
  ! USED BY
  !    ***_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  nzdom,ndime,r_dom,c_dom,r_sol,c_sol,&
       &                        r_sym,c_sym,nzsym,mnodb,npoin
  use def_master, only       :  solve
  use def_kermod, only       :  kfl_ndvars_opt
  implicit none
  integer(ip), intent(in)    :: itask,pnode,pevat
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: elrbu(ndime,pnode)
  real(rp),    intent(in)    :: elrbp(pnode)
!   real(rp),    intent(inout) :: resdiff(kfl_ndvars_opt,*)
  real(rp),    intent(inout) :: sens_mesh(ndime,*)
  real(rp),    intent(in)    :: elvel(ndime, pnode,*)
  real(rp),    intent(in)    :: elpre(pnode)
!   integer(ip), intent(in)    :: idesvar
  integer(ip), intent(in)    :: idime_dof,ipoin_dof
  integer(ip)                :: ndofn,inode,jnode,iposi,jposi,izsym 
  integer(ip)                :: idime,jdime,izsol,jpoin,ipoin,jcolu,ind

  select case(itask)
  
!   case(1_ip)
!      !
!      ! Momentum and continuity RHS only
!      !
!      do inode = 1,pnode
!         ipoin = lnods(inode)
!         do idime = 1,ndime
!            ind = (ipoin-1)*ndime + idime
!            !$OMP ATOMIC
!            resdiff(idesvar,ind) = resdiff(idesvar,ind) + elrbu(idime,inode)
!         end do
!         !$OMP ATOMIC
!         resdiff(idesvar,npoin*ndime + ipoin) = resdiff(idesvar,npoin*ndime + ipoin) + elrbp(inode)
!      end do
     
  case(1_ip)
     !
     ! Momentum and continuity RHS only
     !
     do inode = 1,pnode
       do idime = 1,ndime
         ! dR_m/dX Lambda_m
         !$OMP ATOMIC
         sens_mesh(idime_dof,ipoin_dof) = sens_mesh(idime_dof,ipoin_dof) + elrbu(idime,inode)*elvel(idime,inode,2)
       enddo
       ! dR_c/dX Lambda_c
       !$OMP ATOMIC
       sens_mesh(idime_dof,ipoin_dof) = sens_mesh(idime_dof,ipoin_dof) + elrbp(inode)*elpre(inode)
     end do

  end select

end subroutine nsi_assresdiff
