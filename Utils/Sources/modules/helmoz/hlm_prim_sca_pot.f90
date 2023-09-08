subroutine hlm_prim_sca_pot(pnode,elcod,pscapo)

!-------------------------------------------------------------------------------------------
! Sources/modules/helmoz/hlm_prim_sca_pot.f90
! NAME
!    hlm_prim_sca_pot
! DESCRIPTION
!    This routine calculates values of primary electric scalar potential in element nodes.
!    Description of emmet_hlm values is given in hlm_prim_vec_pot.f90
!
! INPUT ARGUMENTS
!    PNODE .... Number of nodes of an element
!    ELCOD .... Cartesian coordinates of element nodes
! OUTPUT ARGUMENTS
!    PSCAPO .... Values of primary scalar potential in elemet nodes 
! USED BY
!    hlm_elmope
!-------------------------------------------------------------------------------------------

 use def_kintyp, only              :  ip,rp
 use def_domain, only              :  ndime
 use def_helmoz

 implicit none

 !Dummy arguments
 integer(ip), intent(in)          :: pnode
 real(rp)   , intent(in)          :: elcod(ndime,pnode)
 complex(rp), intent(out)         :: pscapo(pnode)

 !Local variables
 integer(ip)                      :: inode

 select case ( emmet_hlm )

 case ( 1_ip )
	do inode = 1,pnode
 		pscapo(inode) = (0.0_rp,0.0_rp)
	end do 
 case ( 2_ip )
	do inode = 1,pnode
		pscapo(inode) = (0.0_rp,0.0_rp)
	end do 
 case ( 3_ip )
	do inode = 1,pnode
		pscapo(inode) = (0.0_rp,0.0_rp)
	end do 
 case ( 4_ip )
	do inode = 1,pnode
		pscapo(inode) = (0.0_rp,0.0_rp)
	end do       
 case ( 5_ip )
	do inode = 1,pnode
		pscapo(inode) = (0.0_rp,0.0_rp)
	end do       
 end select

end subroutine hlm_prim_sca_pot

