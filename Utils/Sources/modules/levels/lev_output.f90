subroutine lev_output()
  !-----------------------------------------------------------------------
  !****f* Levels/lev_output
  ! NAME 
  !    lev_output
  ! DESCRIPTION
  !    End of a time step 
  !    ITASK=0 ... When timemarching is true. There is output or post-process
  !                of results if required.
  !    ITASK=1 ... When timemarching is false. Output and/or post-process of
  !                results is forced if they have not been written previously.
  ! USED BY
  !    lev_iniunk
  !    lev_endste
  !***
  !-----------------------------------------------------------------------
  use  def_parame
  use  def_master
  use  def_domain
  use  def_levels
  implicit none
  integer(ip) :: ivari,ivarp

  do ivarp = 1,nvarp
     ivari = ivarp
     call posdef(11_ip,ivari)
!print*,'a=',ivari
     if( ivari == 4 ) then
        call lev_openfi(5_ip)
        call lev_outint()           
        call lev_openfi(8_ip)
     else
        call lev_outvar(ivari)
     end if
!print*,'b=',ivari

  end do
  !
  ! Deallocate redistantiation
  !
  call lev_memall(10_ip)

  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! Calculations on sets
     !
     call lev_outset()

  end if
!print*,'c=',ivari

end subroutine lev_output
