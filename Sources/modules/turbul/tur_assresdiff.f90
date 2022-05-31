subroutine tur_assresdiff(&
                         itask,pnode,elres_diff,sens_mesh,idime_dof,ipoin_dof,eltur)
  !-----------------------------------------------------------------------
  !****f* mathru/tur_assresdiff
  ! NAME 
  !    tur_assresdiff
  ! DESCRIPTION
  !    
  ! USES
  ! USED BY
  !    ***_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  nzdom,ndime,r_dom,c_dom,r_sol,c_sol,&
       &                        r_sym,c_sym,nzsym,mnodb,npoin
  use def_turbul, only       :  nturb_tur
  use def_master, only       :  solve
  use def_kermod, only       :  kfl_ndvars_opt
  
  implicit none
  integer(ip), intent(in)    :: itask,pnode
  real(rp),    intent(in)    :: elres_diff(pnode)
  real(rp),    intent(inout) :: sens_mesh(ndime,*)
  real(rp),    intent(in)    :: eltur(nturb_tur,pnode,*)
  integer(ip), intent(in)    :: idime_dof,ipoin_dof
  integer(ip)                :: ndofn,inode 
  integer(ip)                :: idime

  select case(itask)
        
  case(1_ip)
     
     do inode = 1,pnode
       ! dR_t/dX Lambda_t
       sens_mesh(idime_dof,ipoin_dof) = sens_mesh(idime_dof,ipoin_dof) + elres_diff(inode)*eltur(1,inode,3)
     end do
     
  end select

end subroutine tur_assresdiff
