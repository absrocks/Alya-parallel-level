subroutine hlm_wls(tpoin,nn,tp,clnod1,clnod2,coefs1,coefs2,coefs3,coefs4)

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_wls.f90
  ! NAME
  !    hlm_wls
  ! DESCRIPTION
  !    This routine performs the weighted, local least squares approximation
  !    of function values from scattered data - 
  !    Weighted Least Squares (WLS) method
  ! INPUT ARGUMENTS
  !    TPOIN ... Cartesian coordinates of a test point (xt, yt, zt)
  !    NN ...... Number of closest mesh nodes to a test point, N
  !    TP ...... Number of a test point
  !    CLNOD1 ... Global numbers of N closest mesh nodes to a test point
  !    CLNOD2 ... Distances between each of the N closest mesh nodes 
  !               and a test point
  ! OUTPUT ARGUMENTS
  !    COEFS1 ... Unknown coefficients to be minimized for interpolation of 
  !               x component of vector secondary potential
  !    COEFS2 ... Unknown coefficients to be minimized for interpolation of 
  !               y component of vector secondary potential
  !    COEFS3 ... Unknown coefficients to be minimized for interpolation of 
  !               z component of vector secondary potential
  !    COEFS4 ... Unknown coefficients to be minimized for interpolation of 
  !               scalar secondary potential 
  ! USES
  !    hlm_findinv
  ! USED BY
  !    hlm_mlsint
  !-----------------------------------------------------------------------

  use def_master
  use def_domain
  use def_helmoz
  use def_parame

  implicit none

  real(rp),    intent(in)  :: tpoin(3)
  integer(ip), intent(in)  :: nn,tp
  integer(ip), intent(in)  :: clnod1(nn)
  real(rp),    intent(in)  :: clnod2(nn)
  complex(rp), intent(out) :: coefs1(4),coefs2(4),coefs3(4),coefs4(4)

  integer(ip) :: ii,jj,kk,error1,error2, airpresent
  real(rp)    :: momat(4,4),invmomat(4,4), momat2(4,4),invmomat2(4,4)
  complex(rp) :: bb1(4),bb2(4),bb3(4),bb4(4)       
  real(rp)    :: w(nn),xx(nn),yy(nn),zz(nn),wx,wy,wz,h
  complex(rp) :: fvx(nn),fvy(nn),fvz(nn),fs(nn)

  momat(:,:) = 0.0_rp
  bb1(:) = (0.0_rp,0.0_rp)
  bb2(:) = (0.0_rp,0.0_rp)
  bb3(:) = (0.0_rp,0.0_rp)
  bb4(:) = (0.0_rp,0.0_rp)
  coefs1(:) = (0.0_rp,0.0_rp)
  coefs2(:) = (0.0_rp,0.0_rp)
  coefs3(:) = (0.0_rp,0.0_rp)
  coefs4(:) = (0.0_rp,0.0_rp)
  error1 = 0
  error2 = 0

  airpresent = 0
  kk = nn * (tp - 1_ip)
  do ii = 1,nn
	xx(ii) = clcoor_hlm(1,kk+ii)
	yy(ii) = clcoor_hlm(2,kk+ii)
	zz(ii) = clcoor_hlm(3,kk+ii)

	if ( IMASTER ) then
		fvx(ii) = smgvp_hlm(1,kk+ii)
		fvy(ii) = smgvp_hlm(2,kk+ii)
		fvz(ii) = smgvp_hlm(3,kk+ii)
		fs(ii)  = selsp_hlm(kk+ii)
	else  ! == if sequential
		fvx(ii) = smgvp_hlm(1,clnod1(ii))
		fvy(ii) = smgvp_hlm(2,clnod1(ii))
		fvz(ii) = smgvp_hlm(3,clnod1(ii))
		fs(ii)  = selsp_hlm(clnod1(ii))
	end if

	if (airpl_hlm /= 1.0e30_rp .and. ((zz(ii)<airpl_hlm .and. tpoin(3)>=airpl_hlm) .or. (zz(ii)>airpl_hlm .and. tpoin(3)<airpl_hlm))) then
		airpresent = airpresent + 1        !Points on the other side of the air-plane are present, need to treat scalar potential in a different way
	endif
  end do

  !Mean of tetra's sides
  h = 0.0_rp
  do ii = 1,4
	do jj = ii+1,4
		h = h + sqrt((xx(ii)-xx(jj)) * (xx(ii)-xx(jj)) + (yy(ii)-yy(jj)) * (yy(ii)-yy(jj)) + (zz(ii)-zz(jj)) * (zz(ii)-zz(jj)))
	end do
  end do
  h = h / 6.0_rp

  !Calculate values of a weighting function
  do ii = 1,nn
	w(ii) = exp(-clnod2(ii)*clnod2(ii)/(h*h))
!	w(ii) = (1-clnod2(ii)/h)**4*(4*clnod2(ii)/h+1)
  end do

  !Creation of the Moment Matrix
  do ii = 1,nn
	wx = w(ii) * xx(ii)
	momat(1,1) = momat(1,1) + wx * xx(ii)
	momat(1,2) = momat(1,2) + wx * yy(ii)
	momat(1,3) = momat(1,3) + wx * zz(ii)
	momat(1,4) = momat(1,4) + wx

	wy = w(ii) * yy(ii) 
	momat(2,1) = momat(2,1) + wy * xx(ii)
	momat(2,2) = momat(2,2) + wy * yy(ii)
	momat(2,3) = momat(2,3) + wy * zz(ii)
	momat(2,4) = momat(2,4) + wy

	wz = w(ii) * zz(ii)
	momat(3,1) = momat(3,1) + wz * xx(ii)
	momat(3,2) = momat(3,2) + wz * yy(ii)
	momat(3,3) = momat(3,3) + wz * zz(ii)
	momat(3,4) = momat(3,4) + wz

	momat(4,1) = momat(4,1) + wx
	momat(4,2) = momat(4,2) + wy
	momat(4,3) = momat(4,3) + wz
	momat(4,4) = momat(4,4) + w(ii)

	bb1(1) = bb1(1) + wx * fvx(ii)
	bb1(2) = bb1(2) + wy * fvx(ii)
	bb1(3) = bb1(3) + wz * fvx(ii)
	bb1(4) = bb1(4) + w(ii) * fvx(ii)

	bb2(1) = bb2(1) + wx * fvy(ii)
	bb2(2) = bb2(2) + wy * fvy(ii)
	bb2(3) = bb2(3) + wz * fvy(ii)
	bb2(4) = bb2(4) + w(ii) * fvy(ii)

	bb3(1) = bb3(1) + wx * fvz(ii)
	bb3(2) = bb3(2) + wy * fvz(ii)
	bb3(3) = bb3(3) + wz * fvz(ii)
	bb3(4) = bb3(4) + w(ii) * fvz(ii)

	bb4(1) = bb4(1) + wx * fs(ii)
	bb4(2) = bb4(2) + wy * fs(ii)
	bb4(3) = bb4(3) + wz * fs(ii)
	bb4(4) = bb4(4) + w(ii) * fs(ii)
  end do

  if (airpresent >= 1) then         ! If the number of air points is bigger than 10%
	momat2(:,:) = 0.0_rp
	bb4(:) = (0.0_rp,0.0_rp)
	do ii = 1,nn
		if (.not.((zz(ii)<airpl_hlm .and. tpoin(3)>=airpl_hlm) .or. (zz(ii)>airpl_hlm .and. tpoin(3)<airpl_hlm))) then
			wx = w(ii) * xx(ii)
			momat2(1,1) = momat2(1,1) + wx * xx(ii)
			momat2(1,2) = momat2(1,2) + wx * yy(ii)
			momat2(1,3) = momat2(1,3) + wx * zz(ii)
			momat2(1,4) = momat2(1,4) + wx

			wy = w(ii) * yy(ii) 
			momat2(2,1) = momat2(2,1) + wy * xx(ii)
			momat2(2,2) = momat2(2,2) + wy * yy(ii)
			momat2(2,3) = momat2(2,3) + wy * zz(ii)
			momat2(2,4) = momat2(2,4) + wy

			wz = w(ii) * zz(ii)
			momat2(3,1) = momat2(3,1) + wz * xx(ii)
			momat2(3,2) = momat2(3,2) + wz * yy(ii)
			momat2(3,3) = momat2(3,3) + wz * zz(ii)
			momat2(3,4) = momat2(3,4) + wz

			momat2(4,1) = momat2(4,1) + wx
			momat2(4,2) = momat2(4,2) + wy
			momat2(4,3) = momat2(4,3) + wz
			momat2(4,4) = momat2(4,4) + w(ii)

			bb4(1) = bb4(1) + wx * fs(ii)
			bb4(2) = bb4(2) + wy * fs(ii)
			bb4(3) = bb4(3) + wz * fs(ii)
			bb4(4) = bb4(4) + w(ii) * fs(ii)
		end if
	end do

	call hlm_findinv(momat2,invmomat2,4,error2)
  endif

  !Find the inverse of the Moment Matrix
  call hlm_findinv(momat,invmomat,4,error1)

  if ( error1 == 0_ip .and. error2 == 0_ip ) then
	!Calculate unknown coefficients
	do ii = 1,4
		do jj = 1,4
			coefs1(ii) = coefs1(ii) + invmomat(ii,jj) * bb1(jj)
			coefs2(ii) = coefs2(ii) + invmomat(ii,jj) * bb2(jj)
			coefs3(ii) = coefs3(ii) + invmomat(ii,jj) * bb3(jj)
			if (airpresent >= 1) then
				coefs4(ii) = coefs4(ii) + invmomat2(ii,jj) * bb4(jj)
			else
				coefs4(ii) = coefs4(ii) + invmomat(ii,jj) * bb4(jj)
			end if
		end do
	end do
  else
	write(*,*)'NOT ABLE TO INTERPOLATE FUNCTION - INVERSE ERROR: ', airpresent, tpoin(1), tpoin(2), tpoin(3)
  end if

end subroutine hlm_wls
