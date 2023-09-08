subroutine qua_turnof()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_turnof
  ! NAME 
  !    qua_turnof
  ! DESCRIPTION
  !    This routine closes the run for the Schrodinger equation
  ! USES
  !    qua_outcpu
  !    qua_output
  !    outfor
  !    outsol
  ! USED BY
  !    Quanty
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_quanty
  use mod_outfor, only : outfor
  implicit none
  !
  ! Output results
  !
  call qua_output(1_ip)
  !
  ! Output solver statistics
  !
  solve_sol => solve_qua(1:)
  call outfor(37_ip,lun_outpu_qua,' ')
  !
  ! Write tail for formatted files
  !
  call outfor(6_ip,lun_outpu_qua,' ')
  !
  ! Output latex file
  !
  call qua_outlat(2_ip)
  !
  ! Close used files
  !
  close(lun_outpu_qua)
  close(lun_conve_qua)
  close(lun_solve_qua)
  close(lun_rstar_qua)
  close(lun_witne_qua)
  close(lun_ppseu_qua)

end subroutine qua_turnof

