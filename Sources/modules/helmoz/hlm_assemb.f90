subroutine hlm_assemb(itask,pnode,lnods,elmat,elrhs)

  !------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_assemb.f90
  ! NAME 
  !    hlm_assemb
  ! DESCRIPTION
  !    This routine assembles elemet equations into the system equations.
  ! INPUT ARGUMENTS
  !    ITASK
  !    PNODE .... Number of element nodes
  !    LNODS .... List of global numbers for element nodes
  !    ELMAT .... Element matrix
  !    ELRHS .... Element RHS 
  ! OUTPUT ARGUMENTS
  !    There are no explicit output arguments. The routine creates the system 
  !    matrix AMATX in BCSR format and the system right-hand-side RHSIX. Both 
  !    of them are declared as pointers in def_master module.
  ! USES
  ! USED BY
  !    hlm_elmope
  !------------------------------------------------------------------------

  use def_kintyp, only     :  ip,rp
  use def_master, only     :  amatx,rhsix,kfl_paral
  use def_domain, only     :  r_sol,c_sol

  implicit none

  integer(ip), intent(in)  :: itask,pnode
  integer(ip), intent(in)  :: lnods(pnode)
  complex(rp), intent(in)  :: elmat(4*pnode,4*pnode),elrhs(4*pnode)

  integer(ip)              :: inode,jnode,izsol,jcolu,ipoin,jpoin
  integer(ip)              :: node,poin,zsol,brsol,ivar,jvar

  !Create system RHS: RHSIX
  do inode = 1,pnode
	node  = 4_ip * (inode - 1_ip)
	ipoin = lnods(inode)     
	poin  = 4_ip * (ipoin - 1_ip)
	rhsix(poin+1) = rhsix(poin+1) + elrhs(node+1)
	rhsix(poin+2) = rhsix(poin+2) + elrhs(node+2)
	rhsix(poin+3) = rhsix(poin+3) + elrhs(node+3)
	rhsix(poin+4) = rhsix(poin+4) + elrhs(node+4) 
  end do

  !Create system matrix: AMATX  
  if ( itask == 1_ip ) then

	do inode = 1,pnode
		ivar = 4_ip * (inode - 1_ip)
		ipoin = lnods(inode)

		do jnode = 1,pnode
			jvar = 4_ip * (jnode - 1_ip)
			jpoin = lnods(jnode)
			izsol = r_sol(ipoin)
			brsol = r_sol(ipoin+1)
			jcolu = c_sol(izsol)

			do while( jcolu /= jpoin .and. izsol < brsol)
				izsol = izsol + 1_ip
				jcolu = c_sol(izsol)
			end do

			if ( izsol < brsol ) then
				zsol = 16_ip * (izsol - 1_ip) 

				amatx(zsol+1) = amatx(zsol+1) + elmat(ivar+1,jvar+1)
				amatx(zsol+2) = (0.0_rp,0.0_rp)
				amatx(zsol+3) = (0.0_rp,0.0_rp)
				amatx(zsol+4) = amatx(zsol+4) + elmat(ivar+1,jvar+4)

				amatx(zsol+5) = (0.0_rp,0.0_rp)
				amatx(zsol+6) = amatx(zsol+6) + elmat(ivar+2,jvar+2)
				amatx(zsol+7) = (0.0_rp,0.0_rp)
				amatx(zsol+8) = amatx(zsol+8) + elmat(ivar+2,jvar+4)

				amatx(zsol+9)  = (0.0_rp,0.0_rp)
				amatx(zsol+10) = (0.0_rp,0.0_rp)
				amatx(zsol+11) = amatx(zsol+11) + elmat(ivar+3,jvar+3)
				amatx(zsol+12) = amatx(zsol+12) + elmat(ivar+3,jvar+4)

				amatx(zsol+13) = amatx(zsol+13) + elmat(ivar+4,jvar+1)
				amatx(zsol+14) = amatx(zsol+14) + elmat(ivar+4,jvar+2)
				amatx(zsol+15) = amatx(zsol+15) + elmat(ivar+4,jvar+3)
				amatx(zsol+16) = amatx(zsol+16) + elmat(ivar+4,jvar+4)
			else
				write (*,*) 'Elementary matrix not assembled! ', ipoin, jpoin
			end if
		end do

	end do

  end if

end subroutine hlm_assemb
