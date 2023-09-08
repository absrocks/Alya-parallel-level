subroutine neu_parall(order)
  !-----------------------------------------------------------------------
  !****f* Parall/neu_parall
  ! NAME
  !    neu_parall
  ! DESCRIPTION
  !    This routine exchange data 
  ! USES
  ! USED BY
  !    Reapro
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_neutro
  use def_domain
  use def_inpout
  use def_solver
  use mod_memory
  use mod_opebcs
  use mod_ADR,   only : ADR_parallel_data
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: dummi

  if( ISEQUEN ) return

  strre = 'neu_paral'
  strin = 'neu_paral'
  strch = 'neu_paral'
  nullify(parin)
  nullify(parre)

  select case ( order )

  case ( 1_ip )

     !-------------------------------------------------------------------
     !
     ! Broadcast physical data read in *.neu.dat file
     !
     !-------------------------------------------------------------------

     do parii = 1,2 
        npari = 0 ; nparr = 0 ; nparc = 0
        !
        ! Physical problem
        !           
        call iexcha(num_energies_neu) 
        call iexcha(num_directions_neu)
        call iexcha(kfl_icosa_neu)
        call iexcha(kfl_snord_neu)
        call rexcha(aniso_neu)
        !
        ! ADR data
        !
        call ADR_parallel_data(ADR_neu)
        !
        ! Send/receive
        ! 
        if( parii == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'PARIN','neu_parall',parin,npari,'DO_NOT_INITIALIZE')
           call memory_alloca(mem_modul(1:2,modul),'PARRE','neu_parall',parre,nparr,'DO_NOT_INITIALIZE')
           if( ISLAVE ) call Parall(2_ip)
        end if
     end do
     if( IMASTER ) call Parall(2_ip)
     call memory_deallo(mem_modul(1:2,modul),'PARRE','neu_parall',parre)
     call memory_deallo(mem_modul(1:2,modul),'PARIN','neu_parall',parin)

  case ( 2_ip ) 

     !-------------------------------------------------------------------
     !
     ! Other data
     !
     !-------------------------------------------------------------------

     call spnbcs(tncod_neu)
     call spbbcs(tbcod_neu)

     do parii = 1,2 
        npari = 0 ; nparr = 0 ; nparc = 0
        !
        ! Numerical treatment
        !
        call iexcha(miinn_neu) 
        call iexcha(kfl_smobo_neu)
        call rexcha(cotol_neu) 
        call rexcha(relax_neu) 
        call rexcha(nitsche_neu)

        solve_sol => solve(1:1)
        call soldef(1_ip)
        !
        ! Postprocess
        !
        call posdef(1_ip,dummi)
        !
        ! Boundary conditions
        !
        !
        ! ADR data
        !
        call ADR_parallel_data(ADR_neu)
        !
        ! Send/receive
        ! 
        if( parii == 1 ) then
           call memory_alloca(mem_modul(1:2,modul),'PARIN','neu_parall',parin,npari,'DO_NOT_INITIALIZE')
           call memory_alloca(mem_modul(1:2,modul),'PARRE','neu_parall',parre,nparr,'DO_NOT_INITIALIZE')
           if( ISLAVE ) call Parall(2_ip)
        end if
     end do
     if( IMASTER ) call Parall(2_ip)
     call memory_deallo(mem_modul(1:2,modul),'PARRE','neu_parall',parre)
     call memory_deallo(mem_modul(1:2,modul),'PARIN','neu_parall',parin)

  end select

  npari = 0
  nparr = 0
  nparc = 0
  nullify(parin)
  nullify(parre)

end subroutine neu_parall

