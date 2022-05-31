
subroutine tem_elmdcost_all(&
		    pnode,pgaus,gpsha,gpvol,elrhs,eltem_forw)
  !------------------------------------------------------------------------
  !****f* Temper/tem_elmdcost_all
  ! NAME 
  !    tem_elmdcost_all
  ! DESCRIPTION
  !    This subroutine calculates elemental contributions of the dcost .. d(f)/d(U)
  !    and add this to elrhs as ELRHS = ELRHS + ELMDCOST
  !    Deponds on the couplement may ELMDCOST = 0 
  ! USES
  ! USED BY
  !    tem_elmope
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use def_kermod, only       :  gasco,kfl_cost_type
!   use def_master, only       :  tempe_forw
  use def_temper, only       :  kfl_regim_tem
  use def_master
  
  implicit none

  integer(ip),    intent(in)        :: pnode,pgaus 
  real(rp),       intent(in)        :: eltem_forw(pnode)
  real(rp),       intent(in)        :: gpvol(pgaus)                          ! |J|*w
  real(rp),       intent(in)        :: gpsha(pnode,pgaus)
  real(rp),       intent(inout)     :: elrhs(pnode)
  real(rp)                          :: elmdcost(pnode)                       ! dF/dT
  real(rp)                          :: ipoin
  real(rp)                          :: elunk(pgaus)
  integer(ip)                       :: inode,igaus

  if (kfl_cost_type == 1) then
    !
    ! Initialization
    !
    elmdcost = 0.0_rp
    elunk = 0.0_rp
    !
    ! calculate dF/dT where F is sum[T-100] all over the domain
    ! 
    do igaus=1,pgaus
	do inode =1,pnode
	  elunk(igaus) = elunk(igaus) + gpsha(inode,igaus)*eltem_forw(inode)
	end do
	do inode =1,pnode
	  elmdcost(inode) = elmdcost(inode) + 2*gpvol(igaus)*gpsha(inode,igaus)*(elunk(igaus)-100.0_rp)
	end do
    end do 
    !
    ! elrhs = elrhs - elmdcost
    !
    do inode =1,pnode
      elrhs(inode) = elrhs(inode) - elmdcost(inode)
    end do
  
  endif
  
end subroutine tem_elmdcost_all
