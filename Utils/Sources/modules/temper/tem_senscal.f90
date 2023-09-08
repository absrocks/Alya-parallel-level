!> tem_senscal.f90
!> @file tem_senscal.f90 
!> @fn tem_senscal 
!> This subroutine calculates the sensitivities using discrete adjoint
!>

subroutine tem_senscal

  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod, only : kfl_ndvars_opt,sens,kfl_dvar_type
  use def_domain
  use def_temper
  implicit none
  integer(ip) :: idesvar,ipoin
  real(rp)    :: resdiff_des(npoin)
  real(rp)    :: inner,dummr
  
  if (kfl_dvar_type == 1) then
    ! initialization
    inner = 0.0_rp
    do ipoin=1,npoin
      resdiff_des(ipoin)=0.0_rp
    end do
	
    ! start sensitivities calculations
    do idesvar = 1,kfl_ndvars_opt
	do ipoin=1,npoin
	  resdiff_des(ipoin)=resdiff_tem(idesvar,ipoin)
	end do
	do ipoin=1,npoin
	  sens(idesvar) = sens(idesvar) + resdiff_des(ipoin)*tempe(ipoin,1)
	enddo
	call pararr('SUM',0_ip,1,sens(idesvar))
    end do

    print*,'sens tem', sens
  endif
  
end subroutine tem_senscal
