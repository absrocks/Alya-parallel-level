subroutine rad_iniunk()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_iniunk
  ! NAME 
  !    rad_iniunk
  ! DESCRIPTION
  !    This routine sets up the initial condition for the radiation.
  ! USED BY
  !    rad_begste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  implicit none
  integer(ip) :: ipoin,icomp

  if( kfl_rstar == 0 ) then  

     if( INOTMASTER ) then
        !
        !  IF we are a test case, prepare mock-up temperature and density
        !
        call rad_initst()
        !
        ! Prepare initial conditions 
        !
        call rad_inicnd()
        !
        ! Smooth boundary conditions
        !
        call rad_smobcs()  !!F Maybe I should smoothe the initial conditions to obey the BCS
        !
        ! When solving an exact solution, impose bc
        !
        call rad_exabcs(1_ip)  !!F Do some exact solutions
        !
        ! Load initial conditions
        !
        if( kfl_rstar == 0 ) then        
           icomp = min(3_ip,ncomp_rad)
           do ipoin = 1,npoin
              if( kfl_fixno_rad(1,ipoin) >= 0 ) then
!!$                 radav_rad(ipoin,icomp) = bvess_rad(ipoin,1)
                 radav_rad(ipoin,1)     = radav_rad(ipoin,icomp)
              end if
           end do
        end if
        !
        ! Latex output format
        !
        call rad_outlat(one)

     end if
  else
     !
     ! Read restart file
     !
     call rad_restar(1_ip)
  end if

end subroutine rad_iniunk
