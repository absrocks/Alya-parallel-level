subroutine lev_assmat(&
     elmat,amatr,pnode,mnode,nequa,lnods,lplev,algso,kfl_timet_lev)
  !-----------------------------------------------------------------------
  !****f* Levels/lev_assmat
  ! NAME 
  !    lev_assmat
  ! DESCRIPTION
  !    This routine performs the assembly of the big matrix
  ! USES
  !    skyase
  !    csrase
  ! USED BY
  !    lev_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_solver
  implicit none
  integer(ip), intent(in)    :: pnode,nequa,mnode,kfl_timet_lev
  integer(ip), intent(in)    :: algso
  integer(ip), intent(in)    :: lplev(nequa)
  integer(ip), intent(in)    :: lnods(pnode)
  real(rp),    intent(in)    :: elmat(mnode,mnode)
  real(rp),    intent(inout) :: amatr(*)

  if(algso==0.and.kfl_timet_lev==2) then
     call skyase(elmat,amatr,pnode,1_ip,1_ip,mnode,nequa,lnods,lplev,2_ip)
  else
     call csrase(elmat,amatr,1_ip,pnode,mnode,lnods,2_ip)
  end if
  
end subroutine lev_assmat
