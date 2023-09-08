
subroutine tem_elmcost_all(&
		      eltem,pnode,pgaus,gpsha,gpvol,lnods)

  !------------------------------------------------------------------------
  !****f* Temper/tem_elmcost_all
  ! NAME 
  !    tem_elmcost_all
  ! DESCRIPTION
  !    This subroutine calculates elemental contributions of the cost .. f
  ! USES
  ! USED BY
  !    tem_elmope
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use def_kermod, only       :  kfl_adj_prob,costf,kfl_cost_type
  use def_master, only       :  prthe,tempe_forw
  use def_temper, only       :  kfl_regim_tem

  implicit none
  
  integer(ip), intent(in)    :: pnode,pgaus
  
  real(rp)    :: elunk(pgaus),eltem_pert(pnode)
  integer(ip) :: inode,igaus,ipoin
  
  real(rp),       intent(in)        :: eltem(pnode),gpvol(pgaus)                       ! |J|*w
  real(rp),       intent(in)        :: gpsha(pnode,pgaus)
  integer(ip),    intent(in)        :: lnods(pnode)
  real(rp)                          :: eltem_forw(pnode)

  if (kfl_cost_type == 1) then
    !
    ! Initialization
    !
    elunk = 0.0_rp
    !
    ! read the temperature from the forward values
    !
    if (kfl_adj_prob == 1_ip) then
      do inode=1,pnode
	ipoin=lnods(inode)
	eltem_forw(inode) = tempe_forw(ipoin,1)
      end do
    else
      do inode=1,pnode
	eltem_forw(inode) = eltem(inode)
      enddo
    endif
    !
    ! calculate F as the sum[T-100] in all the domain
    !
    do igaus=1,pgaus
	do inode =1,pnode
	  elunk(igaus) = elunk(igaus) + gpsha(inode,igaus)*eltem_forw(inode)
	end do
	costf = costf + gpvol(igaus)*(elunk(igaus)-100.0_rp)*(elunk(igaus)-100.0_rp)
    end do
    
  endif
  
end subroutine tem_elmcost_all
