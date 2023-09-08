!> hlm_elm_equs.f90
!> @file hlm_elm_equs.f90 
!> @fn hlm_elm_equs 
!> Text automatically added for recognition of doxygen parser
!>
subroutine hlm_elm_equsdiff(pnode,pgaus,gprea,dgprea,gpvol,gprhs,gpsha,gpcar,elmat,elrhs)

  !--------------------------------------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_elm_equs.f90
  ! NAME
  !	hlm_elm_equs
  ! DESCRIPTION
  !     This routine creates element equations, i.e. element matrix and element RHS, for an element in a mesh.
  ! INPUT ARGUMENTS
  !     PNODE .... Number of element nodes
  !     PGAUS .... Number of element Gauss points
  !     GPREA .... Reaction term: i * omega * mu * sigma(i)
  !     DGPREA ... Reaction term: i * omega * mu * dsigma(i)
  !     GPVOL .... w * |J|(gauss point)
  !     GPRHS .... f (interior, volume, loads)
  !     GPSHA .... Element trial functions
  !     GPCAR .... Cartesian derivatives of element trial functions
  ! OUTPUT ARGUMENTS
  !     ELMAT .... Element matrix
  !     ELRHS .... Element RHS 
  !--------------------------------------------------------------------------------------------------------

  !Declaration statements
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_helmoz, only       :  ncond_hlm
  
  implicit none

  !Dummy arguments
  integer(ip), intent(in)    :: pnode,pgaus
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  complex(rp), intent(in)    :: gprea(ncond_hlm,pgaus)
  complex(rp), intent(in)    :: dgprea(ncond_hlm,pgaus)
  complex(rp), intent(in)    :: gprhs(4,pnode)
  complex(rp), intent(out)   :: elmat(4*pnode,4*pnode)
  complex(rp), intent(out)   :: elrhs(4*pnode)

  !Local variables
  integer(ip)                :: tnode
  integer(ip)                :: inode,jnode,igaus
  integer(ip)                :: ivar,jvar,var

  tnode = 4_ip * pnode       !pnode * 4 unknown parameters in a node (Asx, Asy, Asz, Psis)

  !Initialization: ELMAT and ELRHS

  !Initialization without loop unrolling: ELMAT and ELRHS

  !do jnode=1,tnode
  !	elrhs(jnode) = (0.0_rp,0.0_rp)
  !	do inode=1,tnode
  !      	elmat(inode, jnode) = (0.0_rp,0.0_rp)
  !	end do
  !end do

  !Initialization with loop unrolling: ELMAT and ELRHS

  do jnode = 1,tnode,4

	elrhs(jnode)   = (0.0_rp,0.0_rp)
	elrhs(jnode+1) = (0.0_rp,0.0_rp)
	elrhs(jnode+2) = (0.0_rp,0.0_rp)
	elrhs(jnode+3) = (0.0_rp,0.0_rp)
	do inode = 1,tnode,4
		elmat(inode,   jnode) = (0.0_rp,0.0_rp)
		elmat(inode+1, jnode) = (0.0_rp,0.0_rp)
		elmat(inode+2, jnode) = (0.0_rp,0.0_rp)
		elmat(inode+3, jnode) = (0.0_rp,0.0_rp)
	end do
	do inode = 1,tnode,4
		elmat(inode,   jnode+1) = (0.0_rp,0.0_rp)
		elmat(inode+1, jnode+1) = (0.0_rp,0.0_rp)
		elmat(inode+2, jnode+1) = (0.0_rp,0.0_rp)
		elmat(inode+3, jnode+1) = (0.0_rp,0.0_rp)
	end do
	do inode = 1,tnode,4
		elmat(inode,   jnode+2) = (0.0_rp,0.0_rp)
		elmat(inode+1, jnode+2) = (0.0_rp,0.0_rp)
		elmat(inode+2, jnode+2) = (0.0_rp,0.0_rp)
		elmat(inode+3, jnode+2) = (0.0_rp,0.0_rp)
	end do
	do inode = 1,tnode,4
		elmat(inode,   jnode+3) = (0.0_rp,0.0_rp)
		elmat(inode+1, jnode+3) = (0.0_rp,0.0_rp)
		elmat(inode+2, jnode+3) = (0.0_rp,0.0_rp)
		elmat(inode+3, jnode+3) = (0.0_rp,0.0_rp)
	end do

  end do

  !Create element RHS: ELRHS
  do igaus = 1,pgaus
	do inode = 1,pnode
		var = 4_ip * (inode - 1_ip)
		do jnode = 1,pnode
			elrhs(var+1) = elrhs(var+1) - gpvol(igaus) * dgprea(1,igaus) * &
						   gpsha(inode,igaus) * (gprhs(1,jnode) * gpsha(jnode,igaus))
			elrhs(var+2) = elrhs(var+2) - gpvol(igaus) * dgprea(2,igaus) * &
						   gpsha(inode,igaus) * (gprhs(2,jnode) * gpsha(jnode,igaus))
			elrhs(var+3) = elrhs(var+3) - gpvol(igaus) * dgprea(3,igaus) * &
						   gpsha(inode,igaus) * (gprhs(3,jnode) * gpsha(jnode,igaus))
			elrhs(var+4) = elrhs(var+4) + gpvol(igaus) * (&       !Paper version   !!! gprea originally !!!
						   dgprea(1,igaus) * gprhs(1,jnode) * gpsha(jnode,igaus) * gpcar(1,inode,igaus)&
						 + dgprea(2,igaus) * gprhs(2,jnode) * gpsha(jnode,igaus) * gpcar(2,inode,igaus)&
						 + dgprea(3,igaus) * gprhs(3,jnode) * gpsha(jnode,igaus) * gpcar(3,inode,igaus))
		end do
	end do
  end do

  !Create element matrix: ELMAT
  do igaus = 1,pgaus
	do jnode = 1,pnode

		jvar = 4_ip * (jnode - 1_ip)
		!Create 4 x 4 blocks (submatrices) that form strictly upper triangular part of element matrix
		do inode = 1,jnode-1
			ivar = 4_ip * (inode - 1_ip)
			elmat(ivar+1,jvar+1) = elmat(jvar+1,ivar+1)
			elmat(ivar+2,jvar+1) = (0.0_rp,0.0_rp)
			elmat(ivar+3,jvar+1) = (0.0_rp,0.0_rp)
			elmat(ivar+4,jvar+1) = elmat(ivar+4,jvar+1) + gpvol(igaus) * &
								   gprea(1,igaus) * gpsha(inode,igaus) * gpcar(1,jnode,igaus)
		end do
		do inode = 1,jnode-1
			ivar = 4_ip * (inode - 1_ip)
			elmat(ivar+1,jvar+2) = (0.0_rp,0.0_rp)
			elmat(ivar+2,jvar+2) = elmat(jvar+2,ivar+2)
			elmat(ivar+3,jvar+2) = (0.0_rp,0.0_rp)
			elmat(ivar+4,jvar+2) = elmat(ivar+4,jvar+2) + gpvol(igaus) * &
								   gprea(2,igaus) * gpsha(inode,igaus) * gpcar(2,jnode,igaus)
		end do
		do inode = 1,jnode-1
			ivar = 4_ip * (inode - 1_ip)
			elmat(ivar+1,jvar+3) = (0.0_rp,0.0_rp)
			elmat(ivar+2,jvar+3) = (0.0_rp,0.0_rp)
			elmat(ivar+3,jvar+3) = elmat(jvar+3,ivar+3)
			elmat(ivar+4,jvar+3) = elmat(ivar+4,jvar+3) + gpvol(igaus) * &
								   gprea(3,igaus) * gpsha(inode,igaus) * gpcar(3,jnode,igaus)
		end do
		do inode = 1,jnode-1
			ivar = 4_ip * (inode - 1_ip)
			elmat(ivar+1,jvar+4) = elmat(ivar+4,jvar+1)
			elmat(ivar+2,jvar+4) = elmat(ivar+4,jvar+2)
			elmat(ivar+3,jvar+4) = elmat(ivar+4,jvar+3)
			elmat(ivar+4,jvar+4) = elmat(jvar+4,ivar+4)
		end do

		!Create 4 x 4 blocks (submatrices) that form lower triangular part of element matrix
		do inode = jnode,pnode 
			ivar = 4_ip * (inode - 1_ip)
			elmat(ivar+1,jvar+1) = elmat(ivar+1,jvar+1) - gpvol(igaus) * (&       !or +
								!   gpcar(1,inode,igaus) * gpcar(1,jnode,igaus)&	
								! + gpcar(2,inode,igaus) * gpcar(2,jnode,igaus)&
								! + gpcar(3,inode,igaus) * gpcar(3,jnode,igaus)&
								 - gprea(1,igaus) * gpsha(inode,igaus) * gpsha(jnode,igaus))       !or +
			elmat(ivar+2,jvar+1) = (0.0_rp,0.0_rp)
			elmat(ivar+3,jvar+1) = (0.0_rp,0.0_rp)
			elmat(ivar+4,jvar+1) = elmat(ivar+4,jvar+1) + gpvol(igaus) * &
								   gprea(1,igaus) * gpsha(inode,igaus) * gpcar(1,jnode,igaus)
		end do
		do inode = jnode,pnode
			ivar = 4_ip * (inode - 1_ip)
			elmat(ivar+1,jvar+2) = (0.0_rp,0.0_rp)
			elmat(ivar+2,jvar+2) = elmat(ivar+2,jvar+2) - gpvol(igaus) * (&       !or +
								!   gpcar(1,inode,igaus) * gpcar(1,jnode,igaus)&	
								! + gpcar(2,inode,igaus) * gpcar(2,jnode,igaus)&
								! + gpcar(3,inode,igaus) * gpcar(3,jnode,igaus)&
								 - gprea(2,igaus) * gpsha(inode,igaus) * gpsha(jnode,igaus))       !or +
			elmat(ivar+3,jvar+2) = (0.0_rp,0.0_rp)
			elmat(ivar+4,jvar+2) = elmat(ivar+4,jvar+2) + gpvol(igaus) * &
								   gprea(2,igaus) * gpsha(inode,igaus) * gpcar(2,jnode,igaus)
		end do
		do inode = jnode,pnode
			ivar = 4_ip * (inode - 1_ip)
			elmat(ivar+1,jvar+3) = (0.0_rp,0.0_rp)
			elmat(ivar+2,jvar+3) = (0.0_rp,0.0_rp)
			elmat(ivar+3,jvar+3) = elmat(ivar+3,jvar+3) - gpvol(igaus) * (&       !or +
								!   gpcar(1,inode,igaus) * gpcar(1,jnode,igaus)&	
								! + gpcar(2,inode,igaus) * gpcar(2,jnode,igaus)&
								! + gpcar(3,inode,igaus) * gpcar(3,jnode,igaus)&
								 - gprea(3,igaus) * gpsha(inode,igaus) * gpsha(jnode,igaus))       !or +
			elmat(ivar+4,jvar+3) = elmat(ivar+4,jvar+3) + gpvol(igaus) * &
								   gprea(3,igaus) * gpsha(inode,igaus) * gpcar(3,jnode,igaus)
		end do
		do inode = jnode,pnode
			ivar = 4_ip * (inode - 1_ip)
			elmat(ivar+1,jvar+4) = elmat(ivar+4,jvar+1)
			elmat(ivar+2,jvar+4) = elmat(ivar+4,jvar+2)
			elmat(ivar+3,jvar+4) = elmat(ivar+4,jvar+3)
			elmat(ivar+4,jvar+4) = elmat(ivar+4,jvar+4) - gpvol(igaus) * (&
								   gprea(1,igaus) * gpcar(1,inode,igaus) * gpcar(1,jnode,igaus)&	
								 + gprea(2,igaus) * gpcar(2,inode,igaus) * gpcar(2,jnode,igaus)&
								 + gprea(3,igaus) * gpcar(3,inode,igaus) * gpcar(3,jnode,igaus))
		end do

	end do
  end do
  
end subroutine hlm_elm_equsdiff

