subroutine rad_output()
  !------------------------------------------------------------------------
  !****f* Radiat/rad_output
  ! NAME 
  !    rad_output
  ! DESCRIPTION
  !    Output and postprocess of solution
  ! USES
  ! USED BY
  !    rad_iniunk
  !    rad_endite
  !    rad_endste
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  use mod_iofile
  implicit none
  integer(ip) :: ivari,ivarp
  !
  ! Initial solution, end of a time step and and of run
  !
  do ivarp = 1,nvarp
     ivari = ivarp
     call posdef(11_ip,ivari)
     call rad_outvar(ivari)
  end do

  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! Boundary conditions
     !
     if( npp_bound_rad == 1 ) then 
        npp_bound_rad = 0
        call rad_outbcs()
     end if
     !
     ! Calculations on sets
     !
     call rad_outset()
     !
     ! Calculations on witness points
     !
     call rad_outwit()
     !
     ! Error w/r exact solution
     !
     if( kfl_exacs_rad /= 0 .and. ISEQUEN ) then
        call rad_exaerr(1_ip)
     end if

  else if( ittyp == ITASK_ENDRUN ) then
     !
     ! End of the run
     !
     if( kfl_splot_rad == 1 .and. ISEQUEN ) then
        call suplot(1,radav_rad(1:npoin,1),lun_splot_rad)
        close(lun_splot_rad)
     end if

     if( kfl_psmat_rad > 0 .and. INOTMASTER ) then
        kfl_psmat_rad = 0
        call pspltm(&
             npoin,npoin,1_ip,0_ip,c_dom,r_dom,amatr,&
             trim(title)//': '//namod(modul),0_ip,18.0_rp,'cm',&
             0_ip,0_ip,2_ip,lun_psmat_rad)
        call rad_openfi(4_ip)
     end if
     !
     ! Residual
     !
     !if( postp(1) % npp_stepi(8) > 0 ) call rad_outvar(8_ip)

  end if

end subroutine rad_output
