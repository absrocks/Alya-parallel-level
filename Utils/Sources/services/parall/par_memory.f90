subroutine par_memory(itask)
  !-------------------------------------------------------------------------------
  !****f* parall/par_memory
  ! NAME
  !    par_memory
  ! DESCRIPTION
  !    Allocate memory for partition dimensions and arrays
  ! INPUT
  ! OUTPUT
  ! USED BY
  !    par_create_graph_arrays
  !***
  !-------------------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  use def_parall 
  use mod_memory
  implicit none
  integer(ip), intent(in) :: itask

  select case(itask)

  case(1_ip)
     !
     ! Allocate memory for graph partition arrays
     ! 
     call memory_alloca(mem_servi(1:2,servi),'LEPAR_PAR','par_memory' , lepar_par ,  nelem )
     call memory_alloca(mem_servi(1:2,servi),'LEPER_PAR','par_memory' , leper_par ,  nelem )
     call memory_alloca(mem_servi(1:2,servi),'LEINV_PAR','par_memory' , leinv_par ,  nelem )
     call memory_alloca(mem_servi(1:2,servi),'LBPAR_PAR','par_memory' , lbpar_par ,  nboun )
     call memory_alloca(mem_servi(1:2,servi),'LBPER_PAR','par_memory' , lbper_par ,  nboun )
     call memory_alloca(mem_servi(1:2,servi),'LBINV_PAR','par_memory' , lbinv_par ,  nboun )
     call memory_alloca(mem_servi(1:2,servi),'LNPAR_PAR','par_memory' , lnpar_par ,  npoin )
     call memory_alloca(mem_servi(1:2,servi),'LNPER_PAR','par_memory' , lnper_par ,  npoin )
     call memory_alloca(mem_servi(1:2,servi),'LNINV_PAR','par_memory' , lninv_par ,  npoin )

  case(2_ip)

     call memory_alloca(mem_servi(1:2,servi),'LNINV_LOC','par_memory' , lninv_loc , npoin_total )
     call memory_alloca(mem_servi(1:2,servi),'XLNIN_LOC','par_memory' , xlnin_loc , npart_par+1_ip )

  case(3_ip)

     call memory_alloca(mem_servi(1:2,servi),'BADJ' ,'par_memory' , badj  , gnb+1_ip )
     call memory_alloca(mem_servi(1:2,servi),'BDOM' ,'par_memory' , bdom  , npoin_total-gni )
     call memory_alloca(mem_servi(1:2,servi),'BPOIN','par_memory' , bpoin , npoin_total-gni )

  case(4_ip)

     if( .not. READ_AND_RUN() ) then

        call memory_deallo(mem_servi(1:2,servi),'XLNIN_LOC','par_memory' , xlnin_loc )
        call memory_deallo(mem_servi(1:2,servi),'LNINV_LOC','par_memory' , lninv_loc )
        call memory_deallo(mem_servi(1:2,servi),'LNINV_PAR','par_memory' , lninv_par )
        call memory_deallo(mem_servi(1:2,servi),'LNPER_PAR','par_memory' , lnper_par )
        call memory_deallo(mem_servi(1:2,servi),'LNPAR_PAR','par_memory' , lnpar_par )
        call memory_deallo(mem_servi(1:2,servi),'LBINV_PAR','par_memory' , lbinv_par )
        call memory_deallo(mem_servi(1:2,servi),'LBPER_PAR','par_memory' , lbper_par )
        call memory_deallo(mem_servi(1:2,servi),'LBPAR_PAR','par_memory' , lbpar_par )
        call memory_deallo(mem_servi(1:2,servi),'LEINV_PAR','par_memory' , leinv_par )
        call memory_deallo(mem_servi(1:2,servi),'LEPER_PAR','par_memory' , leper_par )
        call memory_deallo(mem_servi(1:2,servi),'LEPAR_PAR','par_memory' , lepar_par )

     end if
     !
     ! XLNIN_LOC is used for restart: do not deallocate
     !
     call memory_deallo(mem_servi(1:2,servi),'GINDE_PAR','par_memory' , ginde_par )
     call memory_deallo(mem_servi(1:2,servi),'LNEIG_PAR','par_memory' , lneig_par )

     call memory_deallo(mem_servi(1:2,servi),'NELEW_PAR','par_memory' , nelew_par )

  case(5_ip) 

     call memory_deallo(mem_servi(1:2,servi),'BADJ'      ,'par_memory' , badj      )
     call memory_deallo(mem_servi(1:2,servi),'BDOM'      ,'par_memory' , bdom      )
     call memory_deallo(mem_servi(1:2,servi),'BPOIN'     ,'par_memory' , bpoin     )
     call memory_deallo(mem_servi(1:2,servi),'NEIGHDOM'  ,'par_memory' , neighDom  )
     call memory_deallo(mem_servi(1:2,servi),'LCOMM_PAR' ,'par_memory' , lcomm_par )

  case(7_ip)

     call memory_deallo(mem_servi(1:2,servi) ,'IADUAL'    ,'par_memory' , iaDual     )
     call memory_deallo(mem_servi(1:2,servi) ,'JADUAL'    ,'par_memory' , jaDual     )
     call memory_deallo(mem_servi(1:2,servi) ,'TRANSLDUAL','par_memory' , translDual )
     call memory_deallo(mem_servi(1:2,servi) ,'COLOURS'   ,'par_memory' , colours    )

  case(8_ip)

     call memory_alloca(mem_servi(1:2,servi),'GINDE_PAR','par_memory' , ginde_par , 4_ip , npart_par+1_ip )
     call memory_alloca(mem_servi(1:2,servi),'LNEIG_PAR','par_memory' , lneig_par , npart_par      )
     call memory_alloca(mem_servi(1:2,servi),'LEIND_PAR','par_memory' , leind_par , npart_par+1_ip )
     call memory_alloca(mem_servi(1:2,servi),'LBIND_PAR','par_memory' , lbind_par , npart_par+1_ip )
     call memory_alloca(mem_servi(1:2,servi),'LSUBZ_PAR','par_memory' , lsubz_par , npart_par+1_ip )
     !
     ! Allocate memory for dimensions
     !
     call memory_alloca(mem_servi(1:2,servi),'NPOIN_PAR','par_memory' , npoin_par , npart_par )
     call memory_alloca(mem_servi(1:2,servi),'NELEM_PAR','par_memory' , nelem_par , npart_par )
     call memory_alloca(mem_servi(1:2,servi),'NBOUN_PAR','par_memory' , nboun_par , npart_par )
     call memory_alloca(mem_servi(1:2,servi),'NSKEW_PAR','par_memory' , nskew_par , npart_par )
     call memory_alloca(mem_servi(1:2,servi),'SLFBO_PAR','par_memory' , slfbo_par , npart_par )
     call memory_alloca(mem_servi(1:2,servi),'NELEW_PAR','par_memory' , nelew_par , npart_par )

  case(11_ip)

     call memory_alloca(mem_servi(1:2,servi),'SOLVE_SOL(1)%DISPL','par_memory' , solve_sol(1) % displ , npart_par+1_ip )
     call memory_alloca(mem_servi(1:2,servi),'SOLVE_SOL(1)%LCOUN','par_memory' , solve_sol(1) % lcoun , npart_par+1_ip )

  case(-11_ip)

     call memory_deallo(mem_servi(1:2,servi),'SOLVE_SOL(1)%DISPL','par_memory' , solve_sol(1) % displ )
     call memory_deallo(mem_servi(1:2,servi),'SOLVE_SOL(1)%LCOUN','par_memory' , solve_sol(1) % lcoun )

  case(12_ip)

    call memory_alloca(mem_servi(1:2,servi),'SOLVE_SOL(1)%XBIG'  ,'par_memory' , solve_sol(1) % xbig   , &
                       solve_sol(1) % nbig             )
    call memory_alloca(mem_servi(1:2,servi),'SOLVE_SOL(1)%LBIG'  ,'par_memory' , solve_sol(1) % lbig   , &
                       solve_sol(1) % nbig             )
    call memory_alloca(mem_servi(1:2,servi),'SOLVE_SOL(1)%XSMALL','par_memory' , solve_sol(1) % xsmall , &
                       max(1_ip,solve_sol(1) % nsmall) )

  case(13_ip)

     call memory_alloca(mem_servi(1:2,servi),'SOLVE_SOL(1)%DISP4'  ,'par_memory' , solve_sol(1) % disp4 , npart_par+1_ip )
     call memory_alloca(mem_servi(1:2,servi),'SOLVE_SOL(1)%LCOU4'  ,'par_memory' , solve_sol(1) % lcou4 , npart_par+1_ip )

  case(-15_ip)

     call memory_deallo(mem_servi(1:2,servi),'SOLVE_SOL(1)%LBIG','par_memory' , solve_sol(1) % lbig )

  case( 16_ip)

     call memory_alloca(mem_servi(1:2,servi),'XADJDOM','par_memory' , xadjDom , npart_par+1_ip )

  case(-16_ip)

     call memory_deallo(mem_servi(1:2,servi),'XADJDOM','par_memory' , xadjDom )
     call memory_deallo(mem_servi(1:2,servi),'ADJDOM' ,'par_memory' , adjDom  )

  case( 17_ip)

     call memory_alloca(mem_servi(1:2,servi),' ADJDOM','par_memory' , adjDom , xadjDom(npart_par+1_ip)-1_ip )

  case(18_ip)

     call memory_alloca(mem_servi(1:2,servi),'NELEM_TOT','par_memory' , nelem_tot , npart_par+1_ip )
     call memory_alloca(mem_servi(1:2,servi),'NPOIN_TOT','par_memory' , npoin_tot , npart_par+1_ip )
     call memory_alloca(mem_servi(1:2,servi),'NPOIA_TOT','par_memory' , npoia_tot , npart_par+1_ip )
     call memory_alloca(mem_servi(1:2,servi),'NBOUN_TOT','par_memory' , nboun_tot , npart_par+1_ip )

  case(20_ip)

     !
     ! Required for REREAD option to store/read lbper_par array from partition files (GROWSMARTER)
     !
     call memory_alloca(mem_servi(1:2,servi),'LBPER_PAR','par_memory' , lbper_par , nboun_total )

  end select

end subroutine par_memory
