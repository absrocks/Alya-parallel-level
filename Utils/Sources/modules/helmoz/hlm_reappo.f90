subroutine hlm_reappo()

  !-----------------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_reappo.f90
  ! NAME 
  !    hlm_reappo
  ! DESCRIPTION
  !    This routine reads data from the file
  !    1) 'pvecpo.dat' (interpolation mode)
  !       First, it reads values of z and r cylindrical coordinates of nodes of 
  !       a reference mesh, of size nz * nr, in the z-r plane.
  !       Then, it reads values of primary vector potential in these nodes.
  !    2) Or from the file 'pripot.dat' (direct mode).
  !       It reads values of primary vector potential in the mesh nodes
  ! USES
  ! USED BY
  !    hlm_reaphy
  !-----------------------------------------------------------------------------------	

  use def_parame
  use def_inpout
  use def_master
  use def_helmoz
  use def_domain

  implicit none

  integer(ip) :: i,iz,ir,ierror
  real(rp)    :: realp1,imagp1,realp2,imagp2,realp3,imagp3
  character(len=16) :: filename

  select case (ppcod_hlm)
  case (1_ip)
  	filename = 'pvecpo.dat'
  case (2_ip)
  	filename = 'pripot.dat'
  case default
  	write (*,*) 'ERROR in hlm_reappo: invalid value of ppcod_hlm!'
  end select

  open(unit=1,status='old',file=filename,iostat=ierror)

  if ( ierror == 0_ip ) then

	select case (ppcod_hlm)
	case (1_ip)

	  	read(1,*) nz_hlm       !Number of points in the z direction, nz
	  	read(1,*) nr_hlm       !Number of points in the r direction, nr
	
	  	call hlm_memphy(2_ip)
	
	  	do iz = 1,nz_hlm
			read(1,*) z_hlm(iz)       !Create an array of z cylindrical coordinates of nodes of a reference mesh in the z-r plane
	  	end do
	
	  	do ir = 1,nr_hlm
			read(1,*) r_hlm(ir)       !Create an array of r cylindrical coordinates of nodes of a reference mesh in the z-r plane
	  	end do
	
	  	do iz = 1,nz_hlm
			do ir = 1,nr_hlm
				read(1,*) realp1,imagp1
				!Create a matrix of values of primary vector potential in nodes of a reference mesh in the z-r plane
				pvepo_hlm(iz,ir) = cmplx(realp1,imagp1,kind=rp)
			end do
	  	end do

	case (2_ip)

!		write (*,*) 'STARTS'
	  	call hlm_memphy(2_ip)
	  	do i = 1,npoin
			read(1,*) realp1,imagp1,realp2,imagp2,realp3,imagp3
			!Read a value of primary vector potential in the mesh node
			pmgvp_hlm(1,i) = cmplx(realp1,imagp1,kind=rp)
			pmgvp_hlm(2,i) = cmplx(realp2,imagp2,kind=rp)
			pmgvp_hlm(3,i) = cmplx(realp3,imagp3,kind=rp)
! Complex precision test
if (i == 1) then
	write (*,*) 'Doubles:  ', realp1, imagp1
	write (*,*) 'Complex: ', pmgvp_hlm(1,i)
!	write (*,*) 'DComplex:', cmplx(realp1,imagp1,kind=rp)
end if
	  	end do
		write (*,*) npoin, ' lines of primary vector potential read'

	end select
        
	close (unit=1)

  else
	write(*,*) 'ERROR OPENING FILE: ',filename
  end if

end subroutine hlm_reappo

