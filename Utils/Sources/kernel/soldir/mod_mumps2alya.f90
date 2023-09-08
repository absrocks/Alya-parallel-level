!>------------------------------------------------------------------------
!> @addtogroup Direct_Solver
!> @{
!> @file    mod_mumps2alya.f90
!> @author  Guillaume Houzeaux
!> @date    19/03/2018
!> @brief   ToolBox for MUMPS
!> @details ToolBox for MUMPS
!>
!------------------------------------------------------------------------

module mod_mumps2alya

  use def_kintyp,             only : ip,rp,lg
  use def_kintyp,             only : direct_solver_typ
  use def_master,             only : INOTMASTER,IMASTER
  use def_master,             only : IPARALL
  use def_kermod,             only : ndivi
  use def_domain,             only : meshe 
  use mod_matrix,             only : matrix_output_dense_format
  use mod_matrix,             only : matrix_full_to_half
  use mod_renumbering,        only : renumbering_lexical_order_type
  use mod_memory,             only : memory_alloca
  use mod_memory,             only : memory_deallo
  use mod_graphs,             only : graphs_csr_to_coo
  use mod_parall,             only : PAR_COMM_MY_CODE_WM4
  use mod_parall,             only : PAR_CODE_SIZE
  use mod_communications,     only : PAR_SUM
  use mod_communications,     only : PAR_GATHER
  use mod_communications,     only : PAR_GATHERV
  use mod_communications,     only : PAR_SCATTERV
  use mod_communications,     only : PAR_INTERFACE_NODE_EXCHANGE 
  use mod_alya2agmg,          only : mathal2matfinal,calc_lexplode_ndofn

  implicit none 
  private

#ifdef MUMPS
  INTERFACE
     SUBROUTINE DMUMPS( id )
       !!!DEC$ ATTRIBUTES C,REFERENCE,NOMIXED_STR_LEN_ARG, ALIAS:'_DMUMPS'   ::
       !!!DMUMPS
       INCLUDE 'dmumps_struc.h'
       TYPE (DMUMPS_STRUC) :: id
     END SUBROUTINE DMUMPS
  END INTERFACE
#endif
  
  public :: direct_solver_initialization_mumps  ! Symbolic factorization
  public :: direct_solver_factorization_mumps   ! Initialization
  public :: direct_solver_solution_mumps        ! Solution of system Ax=b
  public :: direct_solver_cleaning_mumps        ! Solution of system Ax=b

  integer(ip), pointer   :: ia_ndofn(:)
  integer(ip), pointer   :: ja_ndofn(:)
  integer(ip), pointer   :: lndofn_id(:)
  integer(ip), pointer   :: lndofn_jd(:)
  integer(ip), pointer   :: lndofn_iz(:)  
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-03-19
  !> @brief   Initializaiton of MUMPS
  !> @details Initializaiton of MUMPS
  !> 
  !-----------------------------------------------------------------------

  subroutine direct_solver_initialization_mumps(direct_solver)

    use def_master,         only : kfl_paral
    use mod_communications, only : PAR_COMM_SPLIT
    
    type(direct_solver_typ), intent(inout) :: direct_solver 
    integer(ip)                            :: nzso2,nnnn2
    integer(ip)                            :: nn_own,ndof,nz_loc,ipoin,idime
    integer(ip),             pointer       :: lninv_lex(:)
    integer(ip),             pointer       :: lninv_lex_sca(:)
    integer(4)                             :: PAR_COMM_FINAL4,my_new_rank4,kfl_paral4
    logical(lg)                            :: in_parallel

#ifdef MUMPS
    !
    ! Should we run in parallel
    !
    if( direct_solver % kfl_paral > 0 .and. IPARALL ) then
       in_parallel = .true.
    else
       in_parallel = .false.
    end if
    !
    ! Lexical numbering
    !
    if( in_parallel ) then
       nullify(lninv_lex)
       nullify(lninv_lex_sca)
       !       call renumbering_lexical_order_type(meshe(ndivi),lninv_lex,'WITHOUT HALOS')
       ! extrañamente estaba como la linea de arriba (sin ndofn que es optional pero supuestamenet todo esta preparado para que funcione con el )
       ! mirando renumbering_lexical_order_type y renumbering_lexical_order_node tood parec indicar que estan listas paar trabajr con ndof!!!
       ! Al intentar correr peta -- esta mal paar ndof porque renumbering_lexical_order_node permr(:)    y tendría que ser permr(:,:)  y el size del
       ! primer componente ser ndof . Tla como esta ahora falla en call PAR_INTERFACE_NODE_EXCHANGE(permr,'MAX','IN MY CODE') !!!!!
       ! Guillaume em dijo que a aprtir del miercoles esta libre y lo mira
       ! Por ahora usare un lninv_lex_sca  y obtener lninv_lex a partir de _sca
!       call renumbering_lexical_order_type(meshe(ndivi),lninv_lex,'WITHOUT HALOS',direct_solver % ndof)
       call renumbering_lexical_order_type(meshe(ndivi),lninv_lex_sca,'WITHOUT HALOS',1_ip)

       ndof   = direct_solver % ndof
       if( .not. associated(lninv_lex) .and. INOTMASTER ) then
          call memory_alloca(direct_solver % memor,'lninv_lex','direct_solver_initialization_mumps',lninv_lex, meshe(ndivi) % npoin * ndof)
       end if
       do ipoin = 1,meshe(ndivi) % npoin
          do idime = 1,ndof
             lninv_lex( idime + (ipoin-1) * ndof) = idime + (lninv_lex_sca(ipoin)-1) * ndof
          end do
       end do
    end if
    
    if( INOTMASTER ) then
       !
       ! Initialize MUMPS
       !
       direct_solver % mumps_par % JOB = -1_4
       direct_solver % mumps_par % SYM =  0_4
       if( in_parallel ) then
          direct_solver % mumps_par % PAR =  1
          direct_solver % mumps_par % COMM = PAR_COMM_MY_CODE_WM4
       else
          kfl_paral4 = int(kfl_paral,4)
          if( IPARALL ) then
             call PAR_COMM_SPLIT(kfl_paral4,PAR_COMM_FINAL4,my_new_rank4,'IN MY CODE WITHOUT MASTER')
          end if
          direct_solver % mumps_par % PAR =  1
          direct_solver % mumps_par % COMM = PAR_COMM_FINAL4
       end if
       
       call DMUMPS(direct_solver % mumps_par)
       
       nullify(direct_solver % mumps_par % IRN_LOC)
       nullify(direct_solver % mumps_par % JCN_LOC)
       nullify(direct_solver % mumps_par % A      )
       nullify(direct_solver % mumps_par % A_LOC  )
       nullify(direct_solver % mumps_par % RHS    )  
       !
       ! Graph in COO format (DIRECT_SOLVER % MUMPS_PAR % IRN_LOC,DIRECT_SOLVER % MUMPS_PAR % CJN_LOC)
       !
       if( direct_solver % mumps_par % SYM >= 1 ) then
          if ( direct_solver % ndof > 1_ip ) call runend('MUMPS sym not ready for ndof > 1')
          nz_loc = ( direct_solver % nz - direct_solver % nn ) / 2 + direct_solver % nn
       else
          nz_loc = direct_solver % nz * ( direct_solver % ndof ** 2 )
       end if
       call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % MUMPS_PAR % IRN_LOC','solver_direct_solver_initialize',direct_solver % mumps_par % IRN_LOC,nz_loc)
       call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % MUMPS_PAR % JCN_LOC','solver_direct_solver_initialize',direct_solver % mumps_par % JCN_LOC,nz_loc)

       if ( direct_solver % ndof > 1_ip ) then
          nullify(ia_ndofn)
          nullify(ja_ndofn)
          nullify(lndofn_id)
          nullify(lndofn_jd)
          nullify(lndofn_iz) 
          call memory_alloca(direct_solver % memor,'ia_ndofn','solver_direct_solver_initialize',ia_ndofn,   meshe(ndivi) % npoin * direct_solver % ndof + 1 )
          call memory_alloca(direct_solver % memor,'ja_ndofn','solver_direct_solver_initialize',ja_ndofn,   nz_loc)
          call memory_alloca(direct_solver % memor,'lndofn_id','solver_direct_solver_initialize',lndofn_id, nz_loc)
          call memory_alloca(direct_solver % memor,'lndofn_jd','solver_direct_solver_initialize',lndofn_jd, nz_loc)
          call memory_alloca(direct_solver % memor,'lndofn_iz','solver_direct_solver_initialize',lndofn_iz, nz_loc)
          !
          ! Obtain values just allocated
          !
          ! ver si esto funciona para el caso r_dom_own c_dom_own -- parecerái que si -- probar
          ! Ojo en graphs_csr_to_coo necesito hasta nz y npoin voy a pasar meshe(ndivi) % npoin  en lugar de  meshe(ndivi) % npoin_own
          !
!          call calc_lexplode_ndofn(direct_solver % ndof, meshe(ndivi) % npoin_own, direct_solver % nz , meshe(ndivi) % r_dom_own, meshe(ndivi) % c_dom_own, &   
          !               ia_ndofn,  ja_ndofn, lndofn_id, lndofn_jd, lndofn_iz)


!!!return   !!!!  hhhhhhhhhhhhhhhh    con este anda ok
          ! calc_lexplode_ndofn esta creando un pete horrible con le ruened anterior y uno en direct_solver_factorization_mumps  pasa bien el creado de la matriz
          ! si dejo el de factorization pero pongo este runeend debajo de calc_lexplode_ndofn  peta en nsi_elmope_all



          
          call calc_lexplode_ndofn(direct_solver % ndof, meshe(ndivi) % npoin, direct_solver % nz ,  direct_solver % ia,  direct_solver % ja, &   
               ia_ndofn,  ja_ndofn, lndofn_id, lndofn_jd, lndofn_iz)    ! creo que va andar bien aunque calc_lexplode_ndofn  esta escriat con npoin_own no veo que haya problema
                                                                        ! npoin y los vecrtores de tamaño corecto

          if( direct_solver % mumps_par % SYM >= 1 ) then
             call runend('SYMMETRIC not ready')
          else
             call graphs_csr_to_coo(direct_solver % nn * direct_solver % ndof, 1_ip, ia_ndofn, ja_ndofn, nz_loc, direct_solver % mumps_par % IRN_LOC, direct_solver % mumps_par % JCN_LOC)
          end if
       else  ! ndof ==1
          if( direct_solver % mumps_par % SYM >= 1 ) then
             call graphs_csr_to_coo(direct_solver % nn, 1_ip, direct_solver % ia, direct_solver % ja, nz_loc, direct_solver % mumps_par % IRN_LOC, direct_solver % mumps_par % JCN_LOC,'SYMMETRIC')
          else
             call graphs_csr_to_coo(direct_solver % nn, 1_ip, direct_solver % ia, direct_solver % ja, nz_loc, direct_solver % mumps_par % IRN_LOC, direct_solver % mumps_par % JCN_LOC)
          end if
       end if
       !
       ! Global lexical numbering
       !
       ! hhh esto solo sirve para ndof==1   pero no seria dificil crear lninv_lex   es más tal ves crear un lninv_lex_tmp para el escalar y luego que lninv_lex sea el vectorial
       ! mentira al volver a mirar renumbering_lexical_order_type vi que pasandole ndof (optinal) ya te saca en teoria lninv_lex para ndof
       if( in_parallel ) then
          direct_solver % mumps_par % IRN_LOC = lninv_lex(direct_solver % mumps_par % IRN_LOC)
          direct_solver % mumps_par % JCN_LOC = lninv_lex(direct_solver % mumps_par % JCN_LOC)
       end if       
       !
       ! Dimensions required by MUMPS
       !
       nn_own = direct_solver % nn_own
       ndof   = direct_solver % ndof
       nnnn2  = direct_solver % nn_own
       nzso2  = direct_solver % nz_own

       ! hh ok ndof no sym -- ya puse runened más arriba
       if( in_parallel ) then
          call PAR_SUM(nnnn2,'IN MY CODE WITHOUT MASTER') 
          call PAR_SUM(nzso2,'IN MY CODE WITHOUT MASTER')
          direct_solver % mumps_par % NZ_LOC = nz_loc        
          if( direct_solver % mumps_par % MYID == 0 ) then          
             direct_solver % mumps_par % N  = nnnn2 *   direct_solver % ndof
             if( direct_solver % mumps_par % SYM >= 1 ) then
                direct_solver % mumps_par % NZ = ( ( nzso2 - nnnn2 ) / 2 + nnnn2 ) * ( direct_solver % ndof ** 2 )
             else
                direct_solver % mumps_par % NZ = nzso2 * ( direct_solver % ndof ** 2 )             
             end if
          end if
       else
          direct_solver % mumps_par % N      = direct_solver % nn
          direct_solver % mumps_par % NZ     = nz_loc
          direct_solver % mumps_par % NZ_LOC = nz_loc
       end if
       !
       ! Gathering and scattering array
       !
       if( in_parallel ) then
          if( direct_solver % mumps_par % MYID == 0 ) then
             call memory_alloca(direct_solver % memor,'gatsca_mumps','direct_solver_initialization_mumps',direct_solver % gatsca_mumps,PAR_CODE_SIZE-1_ip,lboun=0_ip)
          else
             call memory_alloca(direct_solver % memor,'gatsca_mumps','direct_solver_initialization_mumps',direct_solver % gatsca_mumps,1_ip)
          end if
          call PAR_GATHER(ndof * nn_own,direct_solver % gatsca_mumps,'IN MY CODE WITHOUT MASTER')
       end if
       
    end if

    if( in_parallel ) call renumbering_lexical_order_type(meshe(ndivi),lninv_lex,'DEALLOCATE')  ! Deallocates lninv_lex

    !--------------------------------------------------------------------
    !
    ! Options
    !
    !--------------------------------------------------------------------

    direct_solver % mumps_par % ICNTL(1)  = 6  ! Output stream for error messages
    direct_solver % mumps_par % ICNTL(2)  = 0  ! Output stream for diagnostic printing, statistics, and warning messages
    direct_solver % mumps_par % ICNTL(3)  = 6  ! Output stream for global information, collected on the host
    direct_solver % mumps_par % ICNTL(4)  = 2  ! Level of printing for error, warning, and diagnostic messages

    direct_solver % mumps_par % ICNTL(5)  = 0  ! Matrix in COO format
    direct_solver % mumps_par % ICNTL(11) = 1  ! Compute statistics
    direct_solver % mumps_par % ICNTL(18) = 3  ! Distributed matrix
    
#endif

  end subroutine direct_solver_initialization_mumps

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-03-20
  !> @brief   MUMPS factorization
  !> @details Analysis and factorization with MUMPS
  !> 
  !-----------------------------------------------------------------------

  subroutine direct_solver_factorization_mumps(direct_solver,aa)

    use def_master, only : kfl_paral
    use def_domain, only : ndime,coord
    use mod_renumbering, only : renumbering_sfc
    
    type(direct_solver_typ), intent(inout) :: direct_solver    
    real(rp), target,        intent(in)    :: aa(direct_solver % ndof ** 2 * direct_solver % nz )
    integer(ip)                            :: nn,ndof,nz_loc,n,nn_own
    logical(lg)                            :: in_parallel
    
    real(rp),                pointer       :: coord_gat(:,:)

#ifdef MUMPS
    
    if( INOTMASTER ) then
       !
       ! Should we run in parallel
       !
       if( direct_solver % kfl_paral > 0 .and. IPARALL ) then
          in_parallel = .true.
       else
          in_parallel = .false.
       end if
       !
       ! Input parameters
       !
       ndof   = direct_solver % ndof
       nn     = direct_solver % nn
       nz_loc = direct_solver % mumps_par % NZ_LOC

       if(associated(direct_solver % mumps_par % A_LOC)) then
          call memory_deallo(direct_solver % memor,'A_LOC','direct_solver_factorization_mumps',direct_solver % mumps_par % A_LOC)
       end if
       call memory_alloca(direct_solver % memor,'A_LOC','direct_solver_factorization_mumps',direct_solver % mumps_par % A_LOC,nz_loc)
       if( direct_solver % mumps_par % SYM >= 1 ) then
          call matrix_full_to_half(nn,ndof,direct_solver % ia,direct_solver % ja,aa,direct_solver % mumps_par % A_LOC)
       else
          if ( ndof == 1 ) then   !
             direct_solver % mumps_par % A_LOC(1:nz_loc) = aa(1:nz_loc)
          else
             
             !               subroutine mathal2matfinal  (nz_own,nz2,ndofn,lndofn_id,lndofn_jd,lndofn_iz,amatr_w_halo,amatr_final)
             ! OJO el segundo argumento es nz_2  -- pero en realidad creo que ya no es neecsario con nz_own dberái ser suficiente.
             ! lo que voy a hacer es pasarle nz_own a nz_2  si funciona lo más correcto luego seria quitar el argumeneto nz_2  -  habrái que ver que lo de agmg
             ! siga andando ok por eso no lo hago ahora
             call mathal2matfinal  (nz_loc/(ndof**2), nz_loc/(ndof**2), ndof, lndofn_id, lndofn_jd, &
                  lndofn_iz, aa(1:nz_loc), direct_solver % mumps_par % A_LOC(1:nz_loc))
             ! SUPONGO aca se podria deeallocar lndofn_id lndofn_jd  -- pero revisar antes
          end if
       end if
       !
       ! Output
       !
       !if( direct_solver % mumps_par % MYID == 0 ) then
       !   write(90+kfl_paral,'(a,10(1x,i10))') 'N=',direct_solver % mumps_par % N
       !   write(90+kfl_paral,'(a,10(1x,i10))') 'NZ=',direct_solver % mumps_par % NZ
       !end if
       !write(90+kfl_paral,'(a,10(1x,i10))')    'NZ_LOC=',direct_solver % mumps_par % NZ_LOC
       !write(90+kfl_paral,'(a,10(1x,es10.3))') 'A_LOC=',aa(1:direct_solver % nz)
       !write(90+kfl_paral,'(a,10(1x,i10))')    'IRN_LOC=',direct_solver % mumps_par % IRN_LOC(1:direct_solver % nz)
       !write(90+kfl_paral,'(a,10(1x,i10))')    'JCN_LOC=',direct_solver % mumps_par % JCN_LOC(1:direct_solver % nz)
       !call matrix_output_dense_format(direct_solver % nn,1_ip,direct_solver % ia,direct_solver % ja,aa)
       !
       ! Metis
       !
       if( 1 == 2 ) then ! BEWARE if you change this you need to resolve !!!call PAR_GATHERV   bellow
          direct_solver % mumps_par % ICNTL(7) = 1
          nullify(coord_gat)
          if( direct_solver % mumps_par % MYID == 0 ) then
             n = direct_solver % mumps_par % N
             call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % MUMPS_PAR % PERM_IN','solver_direct_solver_initialize',direct_solver % mumps_par % PERM_IN,n,'IDENTITY')
             allocate(coord_gat(ndime,n))
          else
             call memory_alloca(direct_solver % memor,'DIRECT_SOLVER % MUMPS_PAR % PERM_IN','solver_direct_solver_initialize',direct_solver % mumps_par % PERM_IN,1_ip,'IDENTITY')
             allocate(coord_gat(1,1))
          end if
          nn_own = direct_solver % nn_own
          direct_solver % gatsca_mumps = direct_solver % gatsca_mumps * ndime
          if( iparall ) then
             ! BEWARE I commented out this line because it complaints:
             ! There is no matching specific subroutine for this generic subroutine call.   [PAR_GATHERV]
             ! in any case it is inside a 1==2 
!!!call PAR_GATHERV(coord(1:ndime,1:nn_own),coord_gat,direct_solver % gatsca_mumps,'IN MY CODE WITHOUT MASTER')
          else
             coord_gat(1:ndime,1:nn_own) = coord(1:ndime,1:nn_own)
          end if
          direct_solver % gatsca_mumps = direct_solver % gatsca_mumps / ndime
          if(direct_solver % mumps_par % MYID == 0 ) then
             call renumbering_sfc(256_ip,coord_gat,direct_solver % mumps_par % PERM_IN)
          end if
          deallocate(coord_gat)
       end if
       ! 
       ! Analysis
       !
       direct_solver % mumps_par % JOB = 1_4
       call DMUMPS(direct_solver % mumps_par)
       !
       ! Factorization
       !
       direct_solver % mumps_par % JOB = 2_4     
       call DMUMPS(direct_solver % mumps_par)

    end if

#endif

  end subroutine direct_solver_factorization_mumps

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-03-20
  !> @brief   MUMPS factorization
  !> @details Analysis and factorization with MUMPS
  !> 
  !-----------------------------------------------------------------------

  subroutine direct_solver_solution_mumps(direct_solver,rhsid,unkno)

    use def_master, only : kfl_paral
    
    type(direct_solver_typ), intent(inout) :: direct_solver    
    real(rp), target,        intent(in)    :: rhsid(direct_solver % ndof * direct_solver % nn)     
    real(rp), target,        intent(out)   :: unkno(direct_solver % ndof * direct_solver % nn)
    integer(ip)                            :: nn,nn_own,nn_tot,ndof
    integer(ip)                            :: nunkn
    logical(lg)                            :: in_parallel

#ifdef MUMPS
    
    if( INOTMASTER ) then
       !
       ! Should we run in parallel
       !
       if( direct_solver % kfl_paral > 0 .and. IPARALL ) then
          in_parallel = .true.
       else
          in_parallel = .false.
       end if
       !
       ! Gather RHS
       !
       ndof   = direct_solver % ndof
       nn     = direct_solver % nn
       nunkn  = ndof * nn
       nn_own = direct_solver % nn_own
       nn_tot = int(direct_solver % mumps_par % N,ip)
       !
       ! Scatter RHS
       !       
       if( direct_solver % mumps_par % MYID == 0 ) then
          call memory_alloca(direct_solver % memor,'direct_solver % mumps_par % rhs','direct_solver_factorization_mumps',direct_solver % mumps_par % RHS,ndof * nn_tot)
       else
          call memory_alloca(direct_solver % memor,'direct_solver % mumps_par % rhs','direct_solver_factorization_mumps',direct_solver % mumps_par % RHS,1_ip)
       end if 
       if( in_parallel ) then
          call PAR_GATHERV(rhsid,direct_solver % mumps_par % RHS,ndof * nn_own,direct_solver % gatsca_mumps,'IN MY CODE WITHOUT MASTER')
       else if( direct_solver % mumps_par % MYID == 0 ) then
          direct_solver % mumps_par % RHS(1:nunkn) = rhsid(1:nunkn)
       end if
       !
       ! Solution
       !
       direct_solver % mumps_par % JOB = 3_4 
       call DMUMPS(direct_solver % mumps_par)
       !
       ! Gather solution
       !
       unkno(1:nunkn) = 0.0_rp
       if( in_parallel ) then
          call PAR_SCATTERV(direct_solver % mumps_par % RHS,unkno,direct_solver % gatsca_mumps,ndof * nn_own,'IN MY CODE WITHOUT MASTER')
          call PAR_INTERFACE_NODE_EXCHANGE(ndof,unkno,'SUM','IN MY CODE')
       else
          unkno(1:nunkn) = direct_solver % mumps_par % RHS(1:nunkn)
       end if
       !
       ! Deallocate
       !
       call memory_deallo(direct_solver % memor,'direct_solver % mumps_par_rhs','direct_solver_factorization_mumps',direct_solver % mumps_par % RHS)
       !
       ! Deallocate
       !
       !if( direct_solver % mumps_par % SYM == 1 ) then
       !   call memory_deallo(direct_solver % memor,'A_LOC','direct_solver_factorization_mumps',direct_solver % mumps_par % A_LOC)
       !else
       !   nullify( direct_solver % mumps_par % A_LOC )
       !end if

    end if

#endif

  end subroutine direct_solver_solution_mumps

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-03-20
  !> @brief   MUMPS factorization
  !> @details Analysis and factorization with MUMPS
  !> 
  !-----------------------------------------------------------------------

  subroutine direct_solver_cleaning_mumps(direct_solver)
    
    type(direct_solver_typ), intent(inout) :: direct_solver    

#ifdef MUMPS
    call memory_deallo(direct_solver % memor,'DIRECT_SOLVER % MUMPS_PAR % IRN_LOC','direct_solver_cleaning_mumps',direct_solver % mumps_par % IRN_LOC)
    call memory_deallo(direct_solver % memor,'DIRECT_SOLVER % MUMPS_PAR % JCN_LOC','direct_solver_cleaning_mumps',direct_solver % mumps_par % JCN_LOC)
    direct_solver % mumps_par % JOB = -2_4     
    call DMUMPS(direct_solver % mumps_par)    
#endif    
    call memory_deallo(direct_solver % memor,'gatsca_mumps'                       ,'direct_solver_cleaning_mumps',direct_solver % gatsca_mumps)
    
  end subroutine direct_solver_cleaning_mumps

end module mod_mumps2alya
!> @}
 
