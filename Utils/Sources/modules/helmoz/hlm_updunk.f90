subroutine hlm_updunk(itask)

  !--------------------------------------------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_updunk.f90
  ! NAME 
  !    hlm_updunk
  ! DESCRIPTION
  !    This routine performs several updates.
  ! USED BY
  !    hlm_solite (itask=1): obtains the initial guess for vector of unknowns UNKNX
  !    hlm_solite (itask=2): vector of unknowns updates smgvr_hlm (secondary magnetic vector potential) and 
  !                                                     selsp_hlm (secondary electric scalar potential) vector
  !--------------------------------------------------------------------------------------------------------------

  use def_parame
  use def_master
  use def_domain
  use def_helmoz

  implicit none

  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,poin,icomp
  
  if ( INOTMASTER ) then

	select case ( itask ) 
	case ( 1_ip )
	!Obtain the initial guess for unknx       
	do icomp = 1,nzrhx
		unknx(icomp) = cmplx(0.0_rp,0.0_rp,kind=rp)
	end do

	case ( 2_ip )
		select case ( ppout_hlm ) 
		case ( 1_ip )
			!write(*,*) 'Output: Primary potentials'
			!Only primary potentials (for testing)
			smgvp_hlm => pmgvp_hlm
			selsp_hlm => pelsp_hlm
		case ( 2_ip )
			!write(*,*) 'Output: Secondary potentials'
			!Secondary potentials only
			do ipoin = 1,npoin
				poin = 4_ip * (ipoin - 1_ip)
				smgvp_hlm(1,ipoin) = unknx(poin+1)
				smgvp_hlm(2,ipoin) = unknx(poin+2)
				smgvp_hlm(3,ipoin) = unknx(poin+3)
				selsp_hlm(  ipoin) = unknx(poin+4)

				!smgvp_hlm(1,ipoin) = rhsix(poin+1)    !For testing
				!smgvp_hlm(2,ipoin) = rhsix(poin+2)
				!smgvp_hlm(3,ipoin) = rhsix(poin+3)
				!selsp_hlm(  ipoin) = rhsix(poin+4)
			end do
		case ( 3_ip )
			!write(*,*) 'Output: Total potentials'
			!Primary and secondary potentials together
			do ipoin = 1,npoin
				poin = 4_ip * (ipoin - 1_ip)
				smgvp_hlm(1,ipoin) = unknx(poin+1) + pmgvp_hlm(1,ipoin)
				smgvp_hlm(2,ipoin) = unknx(poin+2) + pmgvp_hlm(2,ipoin)
				smgvp_hlm(3,ipoin) = unknx(poin+3) + pmgvp_hlm(3,ipoin)
				selsp_hlm(  ipoin) = unknx(poin+4) + pelsp_hlm(  ipoin)
			end do
		end select
	end select

  end if

end subroutine hlm_updunk

