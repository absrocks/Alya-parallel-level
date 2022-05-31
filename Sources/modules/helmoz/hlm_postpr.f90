subroutine hlm_postpr(ivari)

  !-----------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_postpr.f90
  ! NAME 
  !    hlm_postpr
  ! DESCRIPTION
  !    This routine gets values of FE-computed potentials that 
  !    are needed for computation of field vectors.
  ! USES
  ! USED BY
  !    hlm_outvar
  !-----------------------------------------------------------------

  use def_parame
  use def_master
  use def_domain
  use def_helmoz

  implicit none

  integer(ip), intent(in) :: ivari

  integer(ip)             :: nn,ii,jj,kk

  nn = nmlsi_hlm       !Number of closest mesh nodes to a site needed for the MLSI

  select case (ivari)  
    case(1_ip)
      do ii = 1,nsite_hlm
  	    kk = nn * (ii - 1_ip)
	      do jj = 1,nn
				  smgvp_hlm(1,kk+jj) = parx2(1,clnod1_hlm(kk+jj))
          smgvp_hlm(2,kk+jj) = parx2(2,clnod1_hlm(kk+jj))
          smgvp_hlm(3,kk+jj) = parx2(3,clnod1_hlm(kk+jj))
        enddo
      enddo  
     
    case(2_ip)
      do ii = 1,nsite_hlm
  	    kk = nn * (ii - 1_ip)
	      do jj = 1,nn
				  selsp_hlm(kk+jj) = parx1(clnod1_hlm(kk+jj))
        enddo
      enddo  
      
  end select

end subroutine hlm_postpr
