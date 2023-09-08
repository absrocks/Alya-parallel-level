subroutine qua_reanut()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_reanut
  ! NAME 
  !    qua_reanut
  ! DSCRIPTION
  !    This routine reads the numerical treatment for quanty
  ! INPUT
  ! OUTPUT

  ! USES
  !    ecoute
  ! USED BY
  !    qua_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_quanty
  use def_domain
  use def_solver
  use def_inpout
  use mod_ecoute, only :  ecoute
  implicit none

  if( INOTSLAVE ) then
     ! 
     !  Initializations (defaults)
     !
     miinn_qua     =  1                               ! One internal iteration
     kfl_tiacc_qua =  1                               ! First order time integ.
     kfl_tisch_qua =  1                               ! Trapezoidal rule
     cotol_qua     =  0.1_rp                          ! Internal tolerance
     !solve_sol     => solve_qua                       ! Algebraic solver type
     solve_sol     => solve(1:) 
     eigen_sol     => eigen_qua                       ! Eigenvalue solver type
     !
     ! Reach the section
     !
     call ecoute('qua_reanut')
     do while(words(1)/='NUMER')
        call ecoute('qua_reanut')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDNU')
        call ecoute('qua_reanut')

        if(words(1)=='MAXIM') then
           miinn_qua = int(param(2))

        else if(words(1)=='CONVE') then
           cotol_qua = param(2)

        else if(words(1)=='MEZCL') then
           mezcla = param(2)
        else if(words(1)=='OUTPU') then
           noutput = param(2)

        else if(words(1)=='EIGEN') then
           call reaeig(1_ip)

        else if(words(1)=='ALGEB') then
           call reasol(1_ip)

        else if(words(1)=='PRECO') then 
           call reasol(2_ip)

        end if
     end do

  end if

end subroutine qua_reanut
 
