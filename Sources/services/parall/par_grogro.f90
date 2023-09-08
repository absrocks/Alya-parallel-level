subroutine par_grogro()
  !-------------------------------------------------------------------------------
  !****f* parall/par_grogro
  ! NAME
  !    par_grogro
  ! DESCRIPTION
  !    Compute the communication strategy between groups.
  !    All arrays are embedded in structure COMLE(ICOML)%COMMD:
  !    bound_dim .............................. Size of the communication arrays
  !    neights(ii) ............................ Subdomain kk
  !    bound_perm(:) .......................... Permutation
  !    bound_size(ii) -> bound_size(ii+1)-1 ... Communication arrays with kk
  !    Others are temporary arrays.
  ! INPUT
  ! OUTPUT
  ! USED BY
  !***
  !-------------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_solver
  use def_parall
  implicit none
  integer(ip) :: inter,ndual_par
  !
  ! Pointers and level update
  !
  comle(icoml)%ngrou =  solve_sol(1)%ngrou
  comle(icoml)%lgrou => solve_sol(1)%lgrou
  solve_sol(1)%icoml =  icoml

  if( IMASTER .and. kfl_ptask /= 2 ) then
     !
     ! Domain graph: NEIGHDOM, LNEIG_PAR, ADJDOM, XADJDOM
     !
     call par_domgro(inter)         
     !
     ! Dual graph: NBDUAL, TRANSLDUAL, TRANSL, IADUAL, JADUAL
     !                
     call par_duagro(ndual_par)        
     !
     ! Color graph: NBCOLOURS, COLOUR
     !
     call par_colgro(ndual_par) 
     !
     ! Communication scheduling: COMMUNSORT, LCOMM_PAR
     !
     call par_comgro(ndual_par)
     !
     ! Arrays for communications: BADJ, BDOM, BPOIN
     !
     call par_arrgro(inter)
     !
     ! Send communcation strategy
     !
     call par_sengro()
     !
     ! Deallocate memory
     !
     call par_memgro(10_ip)

  else if( ISLAVE ) then
     !
     ! Receive communcation strategy
     !
     call par_sengro()

  end if

  icoml = icoml + 1

end subroutine par_grogro 
