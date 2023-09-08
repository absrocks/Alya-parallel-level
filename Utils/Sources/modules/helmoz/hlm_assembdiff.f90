!> hlm_assembdiff.f90
!> @file hlm_assembdiff.f90 
!> @fn hlm_assembdiff 
!> Text automatically added for recognition of doxygen parser
!>
subroutine hlm_assembdiff(itask,pnode,lnods,elmat,elrhs,indvars)

  !------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_assembdiff.f90
  ! NAME 
  !    hlm_assembdiff
  ! DESCRIPTION
  !    This routine assembles elemet equations into the system equations.
  ! INPUT ARGUMENTS
  !    ITASK
  !    PNODE .... Number of element nodes
  !    LNODS .... List of global numbers for element nodes
  !    ELMAT .... Element matrix
  !    ELRHS .... Element RHS 
  ! OUTPUT ARGUMENTS
  !    There are no explicit output arguments. The routine creates the
  !    derivatives of system matrix DAMAT in BCSR format and the system 
  !    right-hand-side DRHSI. Both of them are declared as pointers in def_master module.
  ! USES
  ! USED BY
  !    hlm_elmope
  !------------------------------------------------------------------------

  use def_kintyp, only     :  ip,rp
  use def_master, only     :  damatx,drhsix,kfl_paral
  use def_domain, only     :  r_sol,c_sol

  implicit none

  integer(ip), intent(in)  :: itask,pnode
  integer(ip), intent(in)  :: lnods(pnode)
  complex(rp), intent(in)  :: elmat(4*pnode,4*pnode),elrhs(4*pnode)
  integer(ip), intent(in)  :: indvars

  integer(ip)              :: inode,jnode,izsol,jcolu,ipoin,jpoin
  integer(ip)              :: node,poin,zsol,brsol,ivar,jvar

  !Create system RHS: DRHSI
  do inode = 1,pnode
	node  = 4_ip * (inode - 1_ip)
	ipoin = lnods(inode)     
	poin  = 4_ip * (ipoin - 1_ip)
	!drhsi(poin+1,indvars) = drhsi(poin+1,indvars) + elrhs(node+1)
	!drhsi(poin+2,indvars) = drhsi(poin+2,indvars) + elrhs(node+2)
	!drhsi(poin+3,indvars) = drhsi(poin+3,indvars) + elrhs(node+3)
	!drhsi(poin+4,indvars) = drhsi(poin+4,indvars) + elrhs(node+4) 
	drhsix(poin+1) = drhsix(poin+1) + elrhs(node+1)
	drhsix(poin+2) = drhsix(poin+2) + elrhs(node+2)
	drhsix(poin+3) = drhsix(poin+3) + elrhs(node+3)
	drhsix(poin+4) = drhsix(poin+4) + elrhs(node+4) 
!!!!For testing
	!drhsi(poin+1) = elrhs(node+1)
	!drhsi(poin+2) = elrhs(node+2)
	!drhsi(poin+3) = elrhs(node+3)
	!drhsi(poin+4) = elrhs(node+4) 
!!!!For testing
  end do

  !Create system matrix: DAMAT  
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

				damatx(zsol+1) = damatx(zsol+1) + elmat(ivar+1,jvar+1)
				damatx(zsol+2) = (0.0_rp,0.0_rp)
				damatx(zsol+3) = (0.0_rp,0.0_rp)
				damatx(zsol+4) = damatx(zsol+4) + elmat(ivar+4,jvar+1)

				damatx(zsol+5) = (0.0_rp,0.0_rp)
				damatx(zsol+6) = damatx(zsol+6) + elmat(ivar+2,jvar+2)
				damatx(zsol+7) = (0.0_rp,0.0_rp)
				damatx(zsol+8) = damatx(zsol+8) + elmat(ivar+4,jvar+2)

				damatx(zsol+9)  = (0.0_rp,0.0_rp)
				damatx(zsol+10) = (0.0_rp,0.0_rp)
				damatx(zsol+11) = damatx(zsol+11) + elmat(ivar+3,jvar+3)
				damatx(zsol+12) = damatx(zsol+12) + elmat(ivar+4,jvar+3)

				damatx(zsol+13) = damatx(zsol+13) + elmat(ivar+1,jvar+4)
				damatx(zsol+14) = damatx(zsol+14) + elmat(ivar+2,jvar+4)
				damatx(zsol+15) = damatx(zsol+15) + elmat(ivar+3,jvar+4)
				damatx(zsol+16) = damatx(zsol+16) + elmat(ivar+4,jvar+4)
			end if
		end do

	end do

  end if

!  !Create system matrix: DAMAT  
!  if ( itask == 1_ip ) then
!
!	do inode = 1,pnode
!		ivar = 4_ip * (inode - 1_ip)
!		ipoin = lnods(inode)
!
!		do jnode = 1,pnode
!			jvar = 4_ip * (jnode - 1_ip)
!			jpoin = lnods(jnode)
!			izsol = r_sol(ipoin)
!			brsol = r_sol(ipoin+1)
!			jcolu = c_sol(izsol)
!
!			do while( jcolu /= jpoin .and. izsol < brsol)
!				izsol = izsol + 1_ip
!				jcolu = c_sol(izsol)
!			end do
!
!			if ( izsol < brsol ) then
!				zsol = 16_ip * (izsol - 1_ip) 
!
!				damat(zsol+1,indvars) = damat(zsol+1,indvars) + elmat(ivar+1,jvar+1)
!				damat(zsol+2,indvars) = (0.0_rp,0.0_rp)
!				damat(zsol+3,indvars) = (0.0_rp,0.0_rp)
!				damat(zsol+4,indvars) = damat(zsol+4,indvars) + elmat(ivar+4,jvar+1)
!
!				damat(zsol+5,indvars) = (0.0_rp,0.0_rp)
!				damat(zsol+6,indvars) = damat(zsol+6,indvars) + elmat(ivar+2,jvar+2)
!				damat(zsol+7,indvars) = (0.0_rp,0.0_rp)
!				damat(zsol+8,indvars) = damat(zsol+8,indvars) + elmat(ivar+4,jvar+2)
!
!				damat(zsol+9,indvars)  = (0.0_rp,0.0_rp)
!				damat(zsol+10,indvars) = (0.0_rp,0.0_rp)
!				damat(zsol+11,indvars) = damat(zsol+11,indvars) + elmat(ivar+3,jvar+3)
!				damat(zsol+12,indvars) = damat(zsol+12,indvars) + elmat(ivar+4,jvar+3)
!
!				damat(zsol+13,indvars) = damat(zsol+13,indvars) + elmat(ivar+1,jvar+4)
!				damat(zsol+14,indvars) = damat(zsol+14,indvars) + elmat(ivar+2,jvar+4)
!				damat(zsol+15,indvars) = damat(zsol+15,indvars) + elmat(ivar+3,jvar+4)
!				damat(zsol+16,indvars) = damat(zsol+16,indvars) + elmat(ivar+4,jvar+4)
!			end if
!		end do
!
!	end do
!
!  end if
end subroutine hlm_assembdiff

