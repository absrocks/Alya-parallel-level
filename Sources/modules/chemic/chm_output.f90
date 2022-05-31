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
  use mod_output_postprocess
  implicit none
  external    :: chm_outvar
  integer(ip) :: ivari,ivarp
  !
  ! Initial solution, end of a time step and and of run
  !
  call output_postprocess_variables(chm_outvar)
  !
  ! Update reference time for time-averaging
  ! 
  do ivarp = 1,nvarp
     if( postp(1) % npp_stepi(ivarp,0) /= 0 ) then
        if( mod(ittim, postp(1) % npp_stepi(ivarp,0) ) == 0 ) then   
           avtim_chm = cutim                                    
        endif
     endif
  end do


  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then 
     !
     ! Calculations on sets
     !
     call chm_outset()
     !
     ! Calculations on witness points
     !
     call chm_outwit()
     
  end if

end subroutine chm_output
