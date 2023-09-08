subroutine qua_doiter()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_doiter
  ! NAME 
  !    qua_doiter
  ! DESCRIPTION
  !    This routine controls the internal loop of the Schrodinger equation.
  ! USES
  !    qua_begite
  !    qua_solite
  !    qua_endite
  ! USED BY
  !    Quanty
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_quanty
  implicit none
  integer(ip) :: ipoin

  if( kfl_stead_qua == 0 ) then
     !
     ! lee los PP, genera la funcion potencial para el caso DFT 
     ! genera la densidad inicial para el caso DFT o ALL_electron
     !
     call qua_begite()   
     !
     ! tengo rhoon(npoin,1), v_pot_ps(npoin)
     !
     do while( kfl_goite_qua == 1 )
        !
        ! Resuelvo la ecuacion de Poisson
        !
        if( kfl_dftgs_qua == 1 .or. kfl_alele_qua == 1 ) then
           call qua_solPoi()         
           write(6,*) 'resuelvo poisson=',kfl_paral
        endif
        !
        ! resulvo ecuacion de K-S or Shrodinger
        !
        call qua_solite()    

        write(6,*) 'resuelvo K-S=',kfl_paral
        !
        ! actualizo rhoon(uno)
        !
        call qua_endite(one)

     end do

     if( kfl_dftgs_qua == 1 .or. kfl_alele_qua == 1 ) then
        !
        ! calculo energias varias
        !
        call qua_integr(2_ip)
     end if

     call qua_endite(two)

  end if

end subroutine qua_doiter


