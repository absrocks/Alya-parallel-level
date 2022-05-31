subroutine hlm_prim_vec_pot(pnode,elcod,lnods,pvecpo)

!----------------------------------------------------------------------------------------------------
! Sources/modules/helmoz/hlm_prim_vec_pot.f90
! NAME
!    hlm_prim_vec_pot
! DESCRIPTION
!    This routine calculates values of primary magnetic vector potential in element nodes.
!
!    emmet_hlm = 1 .... CSEM with current loop 
!                   Induction logging problem formulated using cylindrical geometry and
!                   primary vector potential from a horizontal loop source located in a homogeneous
!                   formation described by a uniform electrical conductivity. 
!    In this case, the routine calculates values of primary vector potential in element nodes   
!    using values of primary vector potential in nodes of a reference mesh in the z-r plane,
!    in the following way:
!    1. transforms Cartesian coordinates of element nodes into cylindrical ones.
!    2. for each element node, it finds the node in a reference mesh in the z-r plane whose
!    r and z coordinates are the closest to the corresponding coordinates of the element node.
!    3. assigns the value of primary vector potential in the found reference node to 
!    the corresponding element node. This is the value of primary potential in the azimuthal
!    direction.
!    4. calculates the values of primary potential in Cartesian coordinates.
!
!    emmet_hlm = 2 .... CSEM with HED
!    emmet_hlm = 3 .... CSEM with VED
!    
!    emmet_hlm = 4 or 5 .... Magnetotellurics is an electromagnetic geophysical method of imaging 
!                                 the earth's subsurface by measuring natural variations of electrical 
!                                 and magnetic fields at the Earth's surface. 
!    The Far Field solution for an electric dipole in the infinite is a plane wave, so we can consider 
!    the natural source as a big electric dipole very far away. In order to get a complete impedance tensor,
!    we need to run the program twice for two orthogonal orientations (polarizations) of the dipole. 
!    For the 1st FWD run (1st polarization), we need to orient the dipole in the NS direction (emmet_hlm = 4) and 
!    for the 2nd FWD run (2nd polarization), in the WE direction (emmet_hlm = 5).
! INPUT ARGUMENTS
!    PNODE .... Number of nodes of an element
!    ELCOD .... Cartesian coordinates of element nodes
! OUTPUT ARGUMENTS
!    PVECPO .... Values of primary vector potential in elemet nodes in the directions of
!                Cartesian coordinates
! USED BY
!    hlm_elmope
!----------------------------------------------------------------------------------------------------

 use def_kintyp, only              :  ip,rp
 use def_domain, only              :  ndime, npoin
 use def_parame 
 use def_helmoz

 implicit none

 !Dummy arguments
 integer(ip), intent(in)          :: pnode
 real(rp)   , intent(in)          :: elcod(ndime,pnode)
 integer(ip), intent(in)          :: lnods(pnode)
 complex(rp), intent(out)         :: pvecpo(ndime,pnode) 

 !Local variables
 integer(ip)                      :: inode,idime,rmini,rmaxi,rind,zmini,zmaxi,zind,ii
 real(rp)                         :: x,y,r,theta,z,rr,kr,Rd,absz
 complex(rp)                      :: imag,ids4,rr5,imso(ncond_hlm),gam(ncond_hlm),expgr(ncond_hlm),skobka(ncond_hlm) !Temp dipole variables
 real(rp)                         :: cycod(ndime,pnode)                     !Cylindrical coordinates, (r,theta,z), of element nodes
 complex(rp)                      :: rc,k                                   !Wave number
 complex(rp)                      :: cypvpo(pnode)                          !Value of Ap(r,z) in element nodes in the azimuthal direction 
 real(rp)                         :: ir1,ir2,iz1,iz2,delta,d1,d2,d3         !Interpolation variables
 integer(ip)                      :: r1ind, z2ind, mtind1, mtind2

!write (*,*) 'ppcod_hlm: ', ppcod_hlm, ', dipole: ', elcur_hlm, length_hlm, ', offsets: ', xoffs_hlm, yoffs_hlm, zoffs_hlm

 select case ( emmet_hlm )

 case ( 1_ip : 3_ip )

  if (ppcod_hlm == 1_ip .or. (ppcod_hlm == 2_ip .and. emmet_hlm == 1)) then

	!Transformation of Cartesian coordinates (x,y,z) of element nodes into cylindrical coordinates (r,theta,z)
	do inode = 1,pnode
		x = elcod(1,inode) - xoffs_hlm
		y = elcod(2,inode) - yoffs_hlm
		cycod(1,inode) = sqrt(x*x+y*y)            !r = sqrt(x^2 + y^2)
		if ( x == 0.0_rp .and. y == 0.0_rp ) then
			cycod(2,inode) = 100.0_rp
		else
			cycod(2,inode) = atan2(y,x)       !theta = arctg(y/x), -pi < theta < pi
		end if
		cycod(3,inode) = elcod(3,inode)           !z = z
	end do

  end if

  if (ppcod_hlm == 1_ip) then

	!Calculation of values of primary magnetic potential, Ap(r,z), in element nodes 
	do inode = 1,pnode
		r = cycod(1,inode)
		z = cycod(3,inode)

		rmini = 1_ip
		rmaxi = nr_hlm
		call rec_bin_srch(r_hlm,r,rmini,rmaxi,rind)       !Search in the array of r coordinates in z-r plane 
														  !for the closest value to the r coordinate of the element node 
		
		zmini = 1_ip
		zmaxi = nz_hlm
		call rec_bin_srch(z_hlm,z,zmini,zmaxi,zind)       !Search in the array of z coordinates in z-r plane
														  !for the closest value to the |z| coordinate of the element node

        	!Barycentric interpolation
        	if ( r < r_hlm(rind) ) then                       !'rind' is the index of the closest value to r 
            	r1ind = rind - 1_ip
        	else
            	r1ind = rind + 1_ip
        	end if
        	if ( z < z_hlm(zind) ) then                       !'zind' is the index of the closest value to z
            	z2ind = zind - 1_ip
        	else
            	z2ind = zind + 1_ip
        	end if
        	ir1 = r_hlm(r1ind)
        	iz1 = z_hlm(zind)
        	ir2 = r_hlm(rind)
        	iz2 = z_hlm(z2ind)
         
			!Barycentric coordinates of a triangle, sum of which is always 1      
			delta = (ir1 - r_hlm(rind)) * (iz2 - z_hlm(zind)) - (ir2 - r_hlm(rind)) * (iz1 - z_hlm(zind))
			d1 = (ir1 * iz2 - ir2 * iz1 + r * (iz1 - iz2) + z * (ir2 - ir1)) / delta
			d2 = (ir2 * z_hlm(zind) - r_hlm(rind) * iz2 + r * (iz2 - z_hlm(zind)) + z * (r_hlm(rind) - ir2)) / delta
			d3 = 1.0_rp - d1 - d2                                           
                                                            
			cypvpo(inode) = d1 * pvepo_hlm(zind,rind) + d2 * pvepo_hlm(zind,r1ind) + d3 * pvepo_hlm(z2ind,rind)

			!cypvpo(inode) = pvepo_hlm(zind,rind)              !No interpolation
 	enddo

 	!Calculation of values of primary magnetic potential in the directions of Cartesian coordinates
 	do inode = 1,pnode
		if (emmet_hlm == 1_ip ) then                       !cypvpo(inode) is in the azimuthal direction
			if (cycod(2,inode) == 100.0_rp) then
	 			pvecpo(1,inode) = (0.0_rp,0.0_rp)          !Aptheta is zero on the z axis
				pvecpo(2,inode) = (0.0_rp,0.0_rp)
			else 
				theta = cycod(2,inode)
		 		pvecpo(1,inode) = -cypvpo(inode) * sin(theta)       !Apx = Aptheta * (-sin(theta))
				pvecpo(2,inode) =  cypvpo(inode) * cos(theta)       !Apy = Aptheta * (cos(theta))
			end if
			pvecpo(3,inode) = (0.0_rp,0.0_rp)                           !Apz = 0 + i0 
		else if (emmet_hlm == 2_ip ) then                  !cypvpo(inode) is in the x direction
			pvecpo(1,inode) = cypvpo(inode)
			pvecpo(2,inode) = (0.0_rp,0.0_rp)
			pvecpo(3,inode) = (0.0_rp,0.0_rp)
		else                                               !cypvpo(inode) is in the z direction
			pvecpo(1,inode) = (0.0_rp,0.0_rp)
			pvecpo(2,inode) = (0.0_rp,0.0_rp)
			pvecpo(3,inode) = cypvpo(inode)
		end if
 	end do

  else if (ppcod_hlm == 2_ip) then

	if (emmet_hlm == 1_ip ) then
		write (*,*) 'DIRECT (ppcod_hlm = 2) is currently not supported for Current Loop (emmet_hlm = 1)'
	else
		do inode = 1,pnode
			pvecpo(1,inode) = pmgvp_hlm(1,lnods(inode))
			pvecpo(2,inode) = pmgvp_hlm(2,lnods(inode))
			pvecpo(3,inode) = pmgvp_hlm(3,lnods(inode))
			if (lnods(inode)<10) write (*,*) pvecpo(1,inode),pvecpo(2,inode),pvecpo(3,inode)
		enddo
	end if

  else if (ppcod_hlm == 3_ip) then                         !Generation of primary potentials here

	if (emmet_hlm == 1_ip ) then
		write (*,*) 'MYSELF (ppcod_hlm = 3) is currently not supported for Current Loop (emmet_hlm = 1)'
	else
		do inode = 1,pnode
			x = elcod(1,inode) - xoffs_hlm
			y = elcod(2,inode) - yoffs_hlm
			z = elcod(3,inode)
			absz = abs(z-zoffs_hlm)
			rr = sqrt(x*x+y*y+absz*absz)
			if ( rr < 1.0e-5_rp ) then
				write (*,*) 'Changed rr:', x, y, absz
				rr = 1.0e-5_rp
			end if

			imag = cmplx(0.0_rp,1.0_rp,kind=rp)
            rr5 = rr*rr*rr*rr*rr
			ids4 = cmplx(0.0_rp,-elcur_hlm*length_hlm/(4.0_rp*pi*anguf_hlm),kind=rp)
			do ii=1,ncond_hlm
				imso(ii) = cmplx(0.0_rp,anguf_hlm*perma_hlm(1)*bckco_hlm(ii),kind=rp)
	            gam(ii) = sqrt(imso(ii))
	            expgr(ii) = exp(imag*gam(ii)*rr)
	            skobka(ii) = (imso(ii)*rr*rr+3.0_rp*imag*gam(ii)*rr-3.0_rp)
			enddo

			if (emmet_hlm == 2_ip ) then                    !HED
				pvecpo(1,inode) = 1/(4.0_rp*pi)*perma_hlm(1)*elcur_hlm*length_hlm*(expgr(1)/rr) + ids4/bckco_hlm(1)*(-expgr(1)/rr5*(x*x*skobka(1)+rr*rr*(1.0_rp-imag*gam(1)*rr)))
				pvecpo(2,inode) = ids4/bckco_hlm(2)*(-x*y*expgr(2)/rr5*skobka(2))
				pvecpo(3,inode) = ids4/bckco_hlm(3)*(x*absz*expgr(3)/rr5*skobka(3))
				if (z-zoffs_hlm>0) pvecpo(3,inode) = -pvecpo(3,inode)
			else                                            !VED
				pvecpo(3,inode) = 1/(4.0_rp*pi)*perma_hlm(1)*elcur_hlm*length_hlm*(expgr(3)/rr) + ids4/bckco_hlm(3)*(-expgr(3)/rr5*((z-zoffs_hlm)*(z-zoffs_hlm)*skobka(3)+rr*rr*(1.0_rp-imag*gam(3)*rr)))
				pvecpo(2,inode) = ids4/bckco_hlm(2)*(-(z-zoffs_hlm)*y*expgr(2)/rr5*skobka(2))
				pvecpo(1,inode) = ids4/bckco_hlm(1)*(abs(x)*(z-zoffs_hlm)*expgr(1)/rr5*skobka(1))
	            if (x>0) pvecpo(1,inode) = -pvecpo(1,inode)
			end if

!			if (lnods(inode)<10) write (*,*) lnods(inode), pvecpo(1,inode),pvecpo(2,inode),pvecpo(3,inode)
!			if (lnods(inode)==1) write (*,*) imso, ids4, rr, skobka

			pmgvp_hlm(1,lnods(inode)) = pvecpo(1,inode)
			pmgvp_hlm(2,lnods(inode)) = pvecpo(2,inode)
			pmgvp_hlm(3,lnods(inode)) = pvecpo(3,inode)
		enddo
	end if

  end if

	! Old approach: formula from the Green Book (Nabighian, 1988)
	!I = 1000.0_rp
	!dS = 100.0_rp
	!Rd = -25.0_rp
	!kr = sqrt(anguf_hlm*perma_hlm(1)*bckco_hlm/2)
	!k  = cmplx(kr,-kr,kind=rp)
 	!do inode = 1,pnode
		!x = elcod(1,inode)
		!y = elcod(2,inode)
		!z = elcod(3,inode)
		!rr = sqrt(x*x+y*y+(z-Rd)*(z-Rd))
		!if ( rr < 1.0e-10_rp ) rr = 1.0e-10_rp
		!rc = cmplx(0.0_rp,-rr,kind=rp)
		!pvecpo(1,inode) = I*dS/(4_ip*pi*rr)*exp(k*rc) 
		!pvecpo(2,inode) = (0.0_rp,0.0_rp)
		!pvecpo(3,inode) = (0.0_rp,0.0_rp)
	!end do

 case ( 4_ip : 5_ip )

	!elcur_hlm = 50.0_rp
	!length_hlm = 1.0_rp
	Rd = 5.0_rp

	kr = sqrt(anguf_hlm*perma_hlm(1)*bckco_hlm(1)/2)
	k  = cmplx(kr,-kr,kind=rp)

	if (emmet_hlm == 4_ip ) then
		mtind1 = 1;
		mtind2 = 2;
	else
		mtind1 = 2;
		mtind2 = 1;
	end if

 	do inode = 1,pnode
		x = elcod(1,inode)
		y = elcod(2,inode)
		z = elcod(3,inode)
   
		rr = sqrt(x*x+(Rd-y)*(Rd-y)+z*z)
		if ( rr < 1.0e-10_rp ) rr = 1.0e-10_rp
		rc = cmplx(0.0_rp,-rr,kind=rp)

		pvecpo(mtind1,inode) = elcur_hlm*length_hlm/(4_ip*pi*rr)*exp(k*rc) 
		pvecpo(mtind2,inode) = (0.0_rp,0.0_rp)
		pvecpo(3,inode) = (0.0_rp,0.0_rp)
	end do

 end select

end subroutine hlm_prim_vec_pot


recursive subroutine rec_bin_srch(xx,val,mini,maxi,ind)

!-----------------------------------------------------------------------------------------
! NAME
!    rec_bin_srch
! DESCRIPTION
!    This routine performs binary search of a sorted array 'xx' and it returns the index 
!    of an element whose value is the closest to a value of 'val'.
! INPUT ARGUMENTS
!     XX   .... Array of real values
!     VAL  .... Real value that is compared to values of array elements
!     MINI .... Index of the first element in array
!     MAXI .... Index of the last element in array
! OUTPUT ARGUMENTS
!     IND  .... Index of an element whose value is the closest to the value of 'val'
!-----------------------------------------------------------------------------------------

 use def_kintyp, only       :  ip,rp

 implicit none

 real(rp), intent(in)      :: xx(*)
 real(rp), intent(in)      :: val
 integer(ip), intent(in)   :: mini,maxi
 integer(ip), intent(out)  :: ind

 integer(ip)               :: mid

 mid = (mini + maxi) / 2

 	if ( xx(mid) > val ) then
		if ( xx(mid-1) <= val ) then
			if ( (val - xx(mid-1)) < (xx(mid) - val)) then
				ind = mid-1
			else
				ind = mid
			end if
		else
			call rec_bin_srch(xx,val,mini,mid-1,ind)
		end if
	else if ( xx(mid) < val ) then
		if ( xx(mid+1) >= val ) then
			if ( (val - xx(mid)) < (xx(mid+1) - val)) then
				ind = mid
			else
				ind = mid+1
			end if
		else
			call rec_bin_srch(xx,val,mid+1,maxi,ind)
		end if
	else 
		ind = mid
	end if

end subroutine rec_bin_srch

