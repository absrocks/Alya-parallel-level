subroutine ker_output()
  !------------------------------------------------------------------------
  !****f* Kermod/ker_output
  ! NAME 
  !    ker_output
  ! DESCRIPTION
  !    Output module variables
  ! USES
  !    ker_outvar
  ! USED BY
  !    ker_endste (itask=1)
  !    ker_turnof (itask=2)
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  implicit none
  integer(ip) :: ivari,ivarp
  !
  ! Postprocess variables
  !
  do ivarp = 1,nvarp
     ivari = ivarp  
     call posdef(11_ip,ivari)
     call ker_outvar(ivari)
  end do

end subroutine ker_output
