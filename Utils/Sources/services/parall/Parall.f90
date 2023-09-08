subroutine Parall(order)
  !-----------------------------------------------------------------------
  !****f* parall/Parall
  ! NAME
  !    Parall
  ! DESCRIPTION
  !    This routine is the bridge for Parall service
  ! USES
  !
  ! USED BY
  !
  !***
  !-----------------------------------------------------------------------
  
  use def_parame
  use def_master
  use def_parall
  use def_domain
  use mod_postpr
  use mod_par_virfil
  use mod_parall,             only : commd
  use mod_parall,             only : PAR_INTEGER,PAR_COMM_WORLD
  use mod_parall,             only : lun_outpu_par
  use mod_outfor,             only : outfor
  implicit none
  
  integer(ip), intent(in) :: order
  integer(ip)             :: ichec,dummi
  real(rp)                :: time0,time1,time2,dummr

  servi=5

  call cputim(time0)

  if( order == -1 ) then

     call cputim(time1)
     call par_initia()                                      ! Initialize MPI
     call par_errors(1_ip)
     call cputim(time2)
     cpu_paral(1)=time2-time1

  else if( order == 0 ) then

     call cputim(time1)
     call par_reapro()                                   ! Read service data
     call par_inidat()                                   ! Send initial data
     call par_sendat(one)                                ! Send/receive data from Reapro
     call par_openfi(zero)
     call par_openfi(one)
     if(kfl_paral==0) call outfor(27_ip,lun_outpu_par,' ')
     call cputim(time2)
     cpu_paral(2)=time2-time1

  else if( order == -2 ) then

     if(nproc_par>1) then                                ! Checkpoint for communication
        call par_livinf(11_ip,' ',ichec)
        call par_chkpoi(ichec)
        if(ichec==0) then
           call par_livinf(12_ip,' ',ichec)
        else
           call runend('MPI IS NOT WORKING WELL')
        end if
     end if

  else if( order == 2000 ) then
     !
     ! Reorder a graph using METIS
     !
     call runend('PARALLL: 2000 IS OBSOLETE')     
     !call par_metis(&
     !     2_ip      , dummi     , dummi     , dummi  , dummi      , dummi      , &
     !     dummi     , dummi     , dummi     , dummi  , nnren_par  , iaren_par  , &
     !     jaren_par , invpr_par , permr_par , dummi  , dummi      , dummi      , &
     !     mem_servi(1:2,servi)  )

  else if( IPARALL .or. kfl_ptask == 0 ) then

     select case ( order )

     case(  1)
        !
        ! Partition mesh
        !
        call par_prepro()                                   ! Partition graph
        call par_errors(2_ip)
        call par_openfi(two)
        !
        ! Output partition info
        !
        if( kfl_virfi_par == 1 ) then                       ! Virtual file
           call par_dumbuf(-1_ip)
           kfl_virfi_par = 0
           !call par_inibuf(-1_ip)
        end if
        !
        ! Info and possibly end the run
        !
        if( PART_AND_WRITE() ) then                         ! Partition only: end of the run
           call outmem()
           call outcpu()
           call outlat(2_ip)
           call outlat(3_ip)
           call cputim(time2)
           cpu_servi(1,servi) =  cpu_servi(1,servi) + time2-time0
           call par_turnof()
           call runend('O.K.!')
        end if
        if( IMASTER .and. .not. READ_AND_RUN() ) then
           if(kfl_freme==0)then
              call memgeo(3_ip)                              ! Deallocate Master geometry memory, except LTYPE, LNODS, COORDS, LELPO, PELPO
           else
              call memgeo(-2_ip)                             ! Deallocate Master geometry memory
           end if
        end if
        if( IMASTER ) then
           call par_memory(4_ip)                            ! Deallocate memory of partition arrays
        end if
        kfl_ptask = 1                                       ! Switch to normal execution
        call vocabu(-1_ip,0_ip,0_ip)

     case(  2)
        call par_broadc()                                   ! Broadcast arrays
     case(  3)
        call par_sendin()                                   ! Send arrays
     case(  4)
        call par_receiv()                                   ! Receive arrays
     case(  5)
        call par_operat(one)                                ! Compute minimum
     case(  6)
        call par_finali(1_ip)                                   ! Finalize MPI - Now 6 is only for the OK case
     case(  7)
        call par_turnof()                                   ! Finalize service
     case(  9)
        call par_operat(three)                              ! Compute sum
     case( 10)
        call par_operat(two)                                ! Compute maximum
     case (11)
        call par_sendat(2_ip)                               ! Parallel partitioning
     case(400)
        call par_slexch()                                   ! Exchange data between slaves
     case( 12)
        call par_sengeo(two)                                ! Send/receive LBOEL
     case( 14)
        call par_openfi(five)                               ! Compute min max and average element volume
     case( 15)
        call par_comset(1_ip)                               ! Communication of sets information
     case( 16)
        call par_solpls(one)                                ! Bridge to Solpls service
     case( 17)
        call par_solpls(two)                                ! Bridge to Solpls service
     case( 18)
        kfl_desti_par=0                                     ! Receive from Master
        call par_receiv()
     case( 19)
        if( IMASTER ) then                                  ! Receives from Slave 1
           kfl_desti_par=1
           call par_receiv()
        else
           kfl_desti_par=0                                  ! Slave 1 sends to master
           call par_sendin()
        end if
     case( 20)
        call par_barrie()                                   ! MPI Barrier
     case( 21)
        pard1=nproc_par-1                                   ! Bridge
        pard2=nsire_par
     case( 22)
        pard1=npart_par                                     ! Bridge
        npart=npart_par
     case( 23)
        kfl_desti_par=pard2                                 ! Bridge
     case( 24)
        kfl_desti_par=0                                     ! Send to master
        call par_sendin()
     case( 25)
        !call par_locnum(pard1,pard2)
        call runend('PARALL: OBSOLETE OPTION 25')
     case( 26)
        call par_cregro()
     case( 27)
        call par_livinf(15_ip,' ',0_ip)
     case( 28)

     case( 29)
        call par_livinf(16_ip,' ',0_ip)
     case( 30)
        call par_livinf(17_ip,' ',0_ip)
     case( 31)
        call par_barrie()
     case( 32)
        call par_slegro()                                   ! Group send/receive
     case( 33)
        call par_grogro()                                   ! Group communcation strategy
     case( 34)
        call par_gatsen()                                   ! Group communcation strategy
     case( 35)
        call par_gatgro()                                   ! Group communcation strategy
     case( 36)
        call par_minxch()                                   ! Compute minimum between slaves
     case( 37)
        call par_brotyp(1_ip)                               ! Broadcast types
     case( 38)
        call par_sendat(6_ip)                               ! Send/receive postprocess data
     case( 39)
        call par_comset(2_ip)                               ! Communication of witness point information
     case( 42)
        call par_sen2ma()                                   ! Send arrays to master
     case( 43)
        call par_sleskc()                                   ! SKCOS
     case( 44)
        call par_comset(3_ip)                               ! Define who owns the witness point
     case( 45)
        call par_comset(4_ip)                               ! tell master slaves' list of witness nodes
     case( 46)
        call par_allgat(1_ip)                               ! AllGather integer
     case( 47)
        call par_allgat(2_ip)                               ! AllGatherv real
     case( 48)
        call par_sendat(7_ip)                               ! Send boundary condition structure
     case( 49)
        call par_allgat(3_ip)                               ! AllGatherv integer
     case( 50)
        call par_schurr()                                   ! Graphs for Schur complement type solvers
     case( 51)
        call par_finali(0_ip)                               ! Abort MPI - idem as 6 but now in a separate case when it is not OK
     case( 52)
        call par_memset(1_ip)                               ! Sets
     case(300)
        call par_mygather()                                 ! Gather slave arrays (Master -> Slave)
     case(301)
        call par_scatte()                                   ! Scatter slave arrays (Slaves -> Master)
     case(302)
        call par_scafil()                                   ! Scatter slave arrays with filter (Slaves -> Master)
     case(401)
        call par_sltake()                                   ! Exchange data between slaves
     case(402)
        call par_slesca()                                   ! Exchange data between slaves
     case(404)
        call par_slexib(1_ip)                               ! IB: identification of travesty nodes
     case(405)
        call par_slexib(2_ip)                               ! IB: interoplation for travesty node
     case(406)
        call par_slexib(3_ip)                               ! IB: hole node identification
     case(407)
        call par_slexca()                                   ! Exchange data between slaves (asynchronous)
     case(408)
        call par_slaves(1_ip)                               ! Exchange INT and REAL arrays between slaves
     case(409)
        call par_slaves(2_ip)                               ! Exchange fringe
     case(410)
        call par_slaves(3_ip)                               ! Exchange fringe
     case(411)
        call par_slaves(4_ip)                               ! Exchange between neighboring slaves
     case(421)
        call par_memory(18_ip)                              ! Mesh division
     case(422)
        call par_submsh()                                   ! Mesh division
     case(423)
        call par_checks()                                   ! Perform some checkings
     case(424)
        call par_slexmi()                                   ! Exchange data between slaves: take min
     case(425)
        call par_slexma()                                   ! Exchange data between slaves: take max
     case(426)
        call par_slesch()                                   ! Exchange data between slaves for Schur complement
     case(427)
        if( IMASTER ) then
           call par_cgschu(solve_sol(1)%ndofn,dummr,&       ! Schur complement solver
                dummr,dummr,dummr,dummr,dummr,dummr,dummr)
        else
           call par_cgschu(solve_sol(1)%ndofn,bbi,&
                bbb,xxi,xxb,aii,aib,abi,abb)
        end if
     case(428)
        call par_slexib(4_ip)                               ! IB: Hole identification in the fringe geometry
     case(429)
        call runend('PARALL: PAR_SLEXFR NO LONGER EXISTS')
        !call par_slexfr()                                  ! Take data from fringe geometry nodes
     case(430)
        call par_slonly()                                   ! Take data from only one subdomain
     case(431)
        call par_chkedg()                                   ! Check if an edge is shared by another subdomain
     case(432)
        call par_exgraf()                                   ! Exchange matrix on interface
     case(433)
        call par_extgra()                                   ! Extend matrix graph including all boundary nodes
     case(501)
        call par_deallo(one)                                ! Deallocate memory for parr1
     case(502)
        call par_deallo(two)                                ! Deallocate memory for parr2
     case(503)
        call par_deallo(three)                              ! Deallocate memory for parr3
     case(504)
        call par_deallo(four)                               ! Deallocate memory for pari1
     case(601)
        call par_alloca(one)                                ! Allocate memory for parr1
     case(602)
        call par_alloca(two)                                ! Allocate memory for parr2
     case(603)
        call par_alloca(three)                              ! Allocate memory for pari1
     case(604)
        call par_alloca(four)                               ! Allocate memory for parr1(PARII)
     case(605)
        call par_alloca(five)                               ! Allocate memory for pari1(NELEM)
     case(606)
        call par_lgface(one)                                ! Construct faces communication arrays
     case(607)
        call par_lgface(two)                                ! Exchange faces arrays
     case(608)
        call par_comzon(1_ip)                               ! Zone-wise Communication array: construct arrays
     case(609)
        call par_comzon(2_ip)                               ! Zone-wise Communication array: check if module is solved
     case(610)
        call par_comzon(3_ip)                               ! Zone-wise Communication array: recover global communication arrays
     case(611)
        call par_comzon(4_ip)                               ! Force global communication arrays and PAR_COMM_MYCODE
     case(612)
        call par_comzon(5_ip)                               ! Recover last communication arrays and MPI communicator
     case(613)
        call par_comzon(6_ip)                               ! Check if my zones solves these modules
     case(701)
        call par_gatr3p(one)                                !
     case(702)
        call par_gatr3p(two)                                !
     case(703)
        call par_gatr3p(three)                              !
     case(704)
        pard1 = nneig                                       ! Lagrangian particles: Bridge
     case(705)
        call par_lagran(1_ip)                               ! Lagrangian particles: size
     case(706)
        call par_lagran(2_ip)                               ! Lagrangian particles: dump
     case(707)
        call par_lagran(3_ip)                               ! Lagrangian particles: check repeated
     case(708)
        call par_lagran(4_ip)                               ! Gatherv for integers
     case(709)
        call par_lagran(5_ip)                               ! Lagrangian particles: size
     case(710)
        call par_lagran(6_ip)                               ! Lagrangian particles: dump
     case(711)
        call par_lagran(7_ip)                               ! Lagrangian particles: check repeated
     case(712)
        call par_lagran(8_ip)                               ! Gatherv for integers
     case(801)
        call par_cregrp(1_ip)                               !
     case(802)
        call par_cregrp(2_ip)                               !
     case(803)
        call par_cregrp(3_ip)                               !
     case(900)
        call par_slequa()                                   ! Equal solution at interface
     case(1001)
        call par_bridge(1_ip)
     case(1002)
        call par_bridge(2_ip)
     case(1003)
        call par_bridge(3_ip)
     case(1004)
        !call par_grodom(1_ip)
     case(1005)
        !call par_grodom(2_ip)
     case(3000)
        !call par_sledif()
     case default
        call runend("Parall: case out of range!")
     end select

  end if

  call cputim(time2)
  cpu_servi(1,5) =  cpu_servi(1,5) + time2-time0

end subroutine Parall

