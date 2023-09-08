subroutine got_output()
  !------------------------------------------------------------------------
  !****f* Gotita/got_output
  ! NAME 
  !    got_output
  ! DESCRIPTION
  !    End of a GOTITA time step 
  !    ITASK=0 ... When timemarching is true. There is output or
  !                post-process of results if required.
  !    ITASK=1 ... When timemarching is false. Output and/or post-process
  !                of results is forced if they have not been written 
  !                previously
  ! USES
  !    got_outvar
  !    postpr
  ! USED BY
  !    got_endste (itask=1)
  !    got_turnof (itask=2)
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_gotita
  use mod_postpr
  use mod_iofile
  implicit none
  integer(ip) :: ivari,ivarp
  !
  ! Initial solution, end of a time step and and of run
  !
  do ivarp = 1,nvarp
     ivari = ivarp
     call posdef(11_ip,ivari)
     call got_outvar(ivari)
  end do

  if( ittyp == ITASK_ENDRUN ) then

     if(kfl_exacs_got/=0.and.kfl_paral==-1) then
        call got_exaerr()
     end if

  end if

end subroutine got_output
