
subroutine nsi_elmcost_all(&
		      elvel,pnode,pgaus,gpsha,gpvol,lnods,costf)

  !------------------------------------------------------------------------
  !****f* Nastin/nsi_elmcost_all
  ! NAME 
  !    nsi_elmcost_all
  ! DESCRIPTION
  !    This subroutine calculates elemental contributions of the cost .. f
  ! USES
  ! USED BY
  !------------------------------------------------------------------------

  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use def_kermod, only       :  kfl_adj_prob,kfl_cost_type

  implicit none
  
  integer(ip), intent(in)    :: pnode,pgaus
  
  real(rp)    :: gpunk(pgaus),elvel_pert(pnode)
  integer(ip) :: inode,igaus,ipoin
  
  real(rp),       intent(in)        :: gpvol(pgaus)                      ! |J|*w
  real(rp),       intent(in)        :: elvel(ndime, pnode),gpsha(pnode,pgaus)
  real(rp),       intent(inout)     :: costf
  integer(ip),    intent(in)        :: lnods(pnode)

  if (kfl_cost_type == 2) then  
    !
    !  Initialization
    !
    do igaus=1,pgaus
	  gpunk(igaus) = 0.0_rp
    end do
    !
    ! read the veloc from the forward values
    !    
    do igaus=1,pgaus
	do inode =1,pnode
	  gpunk(igaus) = gpunk(igaus) + gpsha(inode,igaus)*elvel(1,inode)
	end do
        !$OMP ATOMIC
	costf = costf + gpvol(igaus)*gpunk(igaus)**2
    end do
  
  endif

end subroutine nsi_elmcost_all
