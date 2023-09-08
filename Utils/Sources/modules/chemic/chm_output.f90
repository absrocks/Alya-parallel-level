subroutine chm_output()
  !------------------------------------------------------------------------
  !****f* Temper/chm_output
  ! NAME 
  !    chm_output
  ! DESCRIPTION
  !    End of a time step 
  ! 
  ! ITASK = -1 ... Initial solution
  ! ITASK =  0 ... When timemarching is true. There is output or
  !                post-process of results if required.
  ! ITASK =  1 ... When timemarching is false. Output and/or 
  !                post-process of results is forced if they have not 
  !                been written previously. 
  ! USES
  ! USED BY
  !    chm_iniunk
  !    chm_endite
  !    chm_endste
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  implicit none
  integer(ip) :: ivari,ivarp
  !
  ! Initial solution, end of a time step and and of run
  !
  do ivarp = 1,nvarp
     ivari = ivarp
     call posdef(11_ip,ivari)
     call chm_outvar(ivari)
  end do

  if( postp(1) % npp_stepi(30) /= 0 ) then
     if( mod(ittim, postp(1) % npp_stepi(30) ) == 0 ) then   ! AVERAGE frequency
        avtim_chm = cutim                                    ! Update reference time for time-averaging
     endif
  endif


  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then 
     !
     ! Calculations on sets
     !
     call chm_outset()
     !
     ! Calculations on witness points
     !
     call chm_outwit()
     !
     ! user element intergrals
     !
     call chm_elmusr()

  end if

end subroutine chm_output
