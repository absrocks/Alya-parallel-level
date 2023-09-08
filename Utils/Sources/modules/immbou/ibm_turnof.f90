subroutine ibm_turnof()
  !-----------------------------------------------------------------------
  !****f* Immbou/ibm_turnof
  ! NAME 
  !    ibm_turnof
  ! DESCRIPTION
  !    End a run
  ! USES
  ! USED BY
  !    Immbou
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  implicit none
  integer(ip) :: iimbo

  print *,"1"
  do iimbo = 1,nimbo
     write(90+kfl_paral,'(10(1x,e12.6))') imbou(iimbo)%posil(1:3,1),imbou(iimbo)%posia(1:3,1)
  end do
  print *,"2"

end subroutine ibm_turnof
