!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    mod_performance.f90
!> @author  houzeaux
!> @date    2020-07-01
!> @brief   Performance
!> @details Memory and CPU performance
!-----------------------------------------------------------------------

module mod_performance

  use def_master
  use def_kermod
  use def_domain
  use def_solver
  use mod_communications, only : PAR_AVERAGE
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_GATHER
  use mod_parall,         only : PAR_CODE_SIZE
  use mod_messages,       only : messages_live
  use mod_outfor,         only : outfor
  use mod_maths,          only : maths_solve_overdetermined_system
  use mod_timings,        only : timings_output
  use mod_alya2talp,      only : alya2talp_parallel_efficiency
  use mod_alya2talp,      only : alya2talp_MonitoringRegionStop
  use mod_iofile,         only : iofile_open_unit
  use mod_outfor,         only : outfor
  use mod_memory,         only : kfl_varcount
  use mod_memory,         only : lun_varcount
  use mod_memory,         only : memory_output_variable_counter
  use mod_memory,         only : mem_maxim    
  use mod_perf,           only : init_perf
  use mod_perf,           only : write_perf_line
  use mod_perf,           only : write_perf_header
  use mod_strings,        only : integer_to_string
  use mod_std

  implicit none
  private

  real(rp) :: r_tomax                 ! Max memory
  real(rp) :: max_mem                 ! Should be almost=r_tomax... but taken directly from mod_memory
  real(rp) :: ave_mem                 ! Average maximum memory
  real(rp) :: r_mem_max_modul(mmodu)
  real(rp) :: r_mem_ave_modul(mmodu)
  
  real(rp) :: cpu_total
  real(rp) :: cpu_max_element(mmodu)  
  real(rp) :: cpu_ave_element(mmodu)  
  real(rp) :: lb_ave_element(mmodu)
  real(rp) :: cpu_ave_per_element(mmodu)
  real(rp) :: cpu_max_boundary(mmodu)  
  real(rp) :: cpu_ave_boundary(mmodu)  
  real(rp) :: lb_ave_boundary(mmodu)    
  real(rp) :: cpu_ave_per_boundary(mmodu)
  real(rp) :: cpu_max_node(mmodu)  
  real(rp) :: cpu_ave_node(mmodu)  
  real(rp) :: lb_ave_node(mmodu)    
  real(rp) :: cpu_ave_per_node(mmodu)
  real(rp) :: cpu_max_particle(mmodu)  
  real(rp) :: cpu_ave_particle(mmodu)  
  real(rp) :: lb_ave_particle(mmodu)    
  real(rp) :: cpu_max_min_element(mmodu)  
  real(rp) :: cpu_ave_min_element(mmodu)  
  real(rp) :: lb_ave_min_element(mmodu)
  real(rp) :: cpu_delta_per(mmodu)
  real(rp) :: cpu_max_solver(mmodu)
  real(rp) :: cpu_ave_solver(mmodu)
  real(rp) :: lb_ave_solver(mmodu)
  real(rp) :: cpu_max_post(mmodu)
  real(rp) :: cpu_ave_post(mmodu)
  real(rp) :: lb_ave_post(mmodu)        
  real(rp) :: lb_eff(3,0:mmodu)

  real(rp) :: time_comp_ave(0:mmodu)
  real(rp) :: time_comp_max(0:mmodu)
  real(rp) :: time_mpi_ave(0:mmodu)
  real(rp) :: time_mpi_max(0:mmodu)

  public :: performance_outcpu
  public :: performance_outmem
  public :: performance_output
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-01
  !> @brief   Performance
  !> @details Performance output in csv format
  !> 
  !-----------------------------------------------------------------------
  
  subroutine performance_output()
    
    integer(ip) :: imodu,ivari

    if( INOTSLAVE ) then
       call init_perf()
       call write_perf_header()
       call write_perf_line('kernel', 'all',                            'average time'              , cpu_total)
       call write_perf_line('kernel', 'all',                            'load balance'              , lb_eff(1,0))
       call write_perf_line('kernel', 'all',                            'communication efficiency'  , lb_eff(2,0))
       call write_perf_line('kernel', 'all',                            'parallel efficiency'       , lb_eff(3,0))
       call write_perf_line('kernel', 'all',                            'master memory'             , r_tomax)
       call write_perf_line('kernel', 'all',                            'average computation time'  , time_comp_ave(0))
       call write_perf_line('kernel', 'all',                            'average communication time', time_MPI_ave(0))
       call write_perf_line('kernel', 'all',                            'maximum computation time'  , time_comp_max(0))
       call write_perf_line('kernel', 'all',                            'maximum communication time', time_MPI_max(0))
       call write_perf_line('kernel', 'all',                            'maximum memory'            , max_mem)
       call write_perf_line('kernel', 'all',                            'average maximum memory'    , ave_mem)

       do imodu = 1,mmodu-1
          if( kfl_modul(imodu) /= 0 ) then
             !
             ! General
             !
             call write_perf_line(namod(imodu), 'all',                  'average time'              , cpu_modul(CPU_TOTAL_MODULE,imodu))
             call write_perf_line(namod(imodu), 'all',                  'load balance'              , lb_eff(1,imodu))
             call write_perf_line(namod(imodu), 'all',                  'communication efficiency'  , lb_eff(2,imodu))
             call write_perf_line(namod(imodu), 'all',                  'parallel efficiency'       , lb_eff(3,imodu))
             call write_perf_line(namod(imodu), 'all',                  'maximum memory'            , r_mem_max_modul(imodu))
             call write_perf_line(namod(imodu), 'all',                  'average memory'            , r_mem_ave_modul(imodu))
             call write_perf_line(namod(imodu), 'all',                  'average computation time'  , time_comp_ave(imodu))
             call write_perf_line(namod(imodu), 'all',                  'average communication time', time_MPI_ave(imodu))
             call write_perf_line(namod(imodu), 'all',                  'maximum computation time'  , time_comp_max(imodu))
             call write_perf_line(namod(imodu), 'all',                  'maximum communication time', time_MPI_max(imodu))
             !
             ! Assembly
             !
             if( cpu_modul(CPU_COUNT_ASSEMBLY,imodu) > 0.5_rp ) then
                call write_perf_line(namod(imodu), 'element assembly',  'average time', cpu_ave_element(imodu))
                call write_perf_line(namod(imodu), 'element assembly',  'maximum time', cpu_max_element(imodu))
                call write_perf_line(namod(imodu), 'element assembly',  'load balance', lb_ave_element(imodu))
                call write_perf_line(namod(imodu), 'element assembly',  'variability' , cpu_delta_per(imodu))
             end if
             if( cpu_modul(CPU_COUNT_BOUNDARY,imodu) > 0.5_rp ) then
                call write_perf_line(namod(imodu), 'boundary assembly', 'average time', cpu_ave_boundary(imodu))
                call write_perf_line(namod(imodu), 'boundary assembly', 'maximum time', cpu_max_boundary(imodu))
                call write_perf_line(namod(imodu), 'boundary assembly', 'load balance', lb_ave_boundary(imodu))
             end if
             if( cpu_modul(CPU_COUNT_NODE,imodu) > 0.5_rp ) then
                call write_perf_line(namod(imodu), 'node assembly',     'average time', cpu_ave_node(imodu))
                call write_perf_line(namod(imodu), 'node assembly',     'maximum time', cpu_max_node(imodu))
                call write_perf_line(namod(imodu), 'node assembly',     'load balance', lb_ave_node(imodu))
             end if
             if( cpu_modul(CPU_COUNT_PARTICLE,imodu) > 0.5_rp ) then
                call write_perf_line(namod(imodu), 'particle assembly', 'average time', cpu_ave_particle(imodu))
                call write_perf_line(namod(imodu), 'particle assembly', 'maximum time', cpu_max_particle(imodu))
                call write_perf_line(namod(imodu), 'particle assembly', 'load balance', lb_ave_particle(imodu))
             end if
             !
             ! Output
             !
             call write_perf_line(namod(imodu), 'output',               'average time', cpu_ave_post(imodu))
             call write_perf_line(namod(imodu), 'output',               'maximum time', cpu_max_post(imodu))
             !
             ! Solvers
             !
             solve_sol => momod(imodu) % solve
             if( associated(solve_sol) ) then
                do ivari = 1,size(solve_sol)
                   if( solve_sol(ivari) % kfl_algso/=-999 .and. solve_sol(ivari) % nsolv > 0 ) then
                      call write_perf_line(namod(imodu), 'solver '//integer_to_string(ivari), 'average time', lb_ave_solver(imodu))
                      call write_perf_line(namod(imodu), 'spmv '//integer_to_string(ivari)  , 'average time', solve_sol(1) % cpu_spmv(8))
                      call write_perf_line(namod(imodu), 'spmv '//integer_to_string(ivari)  , 'maximum time', solve_sol(1) % cpu_spmv(9))
                   end if
                end do
             end if
          end if
       end do
    end if

  end subroutine performance_output

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    15/07/2015
  !> @brief   Output general CPU time info
  !> @details Output general CPU time info about starting operation and
  !>          modules
  !>
  !-----------------------------------------------------------------------

  subroutine performance_outcpu()

    real(rp)             :: cpu_refer,cpu_minim,cpu_denom,error
    real(rp)             :: cpu_modut
    real(rp)             :: dummr
    integer(ip)          :: imodu,ipass,kelty,ielty,jelty,nn,mm
    real(rp)             :: PE,PE_rms,LB,LB_rms
    real(rp)             :: CE,CE_rms
    integer(ip), pointer :: lperm(:)
    integer(ip), pointer :: numel(:)
    integer(ip), pointer :: numel_gat(:,:)
    real(rp),    pointer :: time_gat(:)
    real(rp),    pointer :: aa(:,:)
    real(rp),    pointer :: bb(:)
    real(rp),    pointer :: xx(:)  

    !----------------------------------------------------------------------
    !
    ! Compute relative weights of elements during the assembly... this
    ! is just an approximation, as it uses the different timings coming from
    ! the partitions to compute a least square problem:
    !
    ! a11 * TET04 + a12 * PYR05 + a13 * PEN06 = t1
    ! a21 * TET04 + a22 * PYR05 + a23 * PEN06 = t2
    ! ...
    ! an1 * TET04 + an2 * PYR05 + an3 * PEN06 = tn
    !
    ! Note that it does not make sense in sequential
    !
    !----------------------------------------------------------------------

    if( IPARALL ) then

       nullify(lperm,numel,numel_gat,time_gat,aa,xx,bb)
       allocate( lperm(nelty) )
       kelty = 0
       nn    = PAR_CODE_SIZE-1
       do ielty = iesta_dom,iesto_dom
          if( lexis(ielty) /= 0 ) then
             kelty = kelty + 1
             lperm(kelty) = ielty
          end if
       end do
       if( kelty > 1 ) then
          allocate( numel(kelty)          )
          allocate( numel_gat(kelty,0:nn) )
          allocate( time_gat(0:nn)        )
          numel = 0
          if( INOTMASTER .and. nelem > 0 ) then
             do jelty = 1,kelty
                ielty = lperm(jelty)
                numel(jelty) = count(ltype(1:nelem)==ielty,KIND=ip)
             end do
          end if
          call PAR_GATHER(numel,numel_gat,'IN THE WORLD')
          mm = kelty
          if( INOTSLAVE ) then
             allocate( aa(nn,mm) )
             allocate( bb(nn)    )
             allocate( xx(mm)    )
          end if
          call outfor(96_ip)
          do imodu = 1,mmodu-1
             if( kfl_modul(imodu) /= 0 ) then
                dummr = cpu_modul(CPU_ASSEMBLY,imodu)
                call PAR_GATHER(dummr,time_gat,'IN THE WORLD')           
                if( INOTSLAVE ) then
                   bb(1:nn) = time_gat(1:nn)
                   do jelty = 1,kelty
                      aa(1:nn,jelty) = real(numel_gat(jelty,1:nn),rp)
                   end do
                   call maths_solve_overdetermined_system(nn,mm,aa,bb,xx,error)
                   if( error >= 0.0_rp ) then
                      xx            = xx/(xx(1)+zeror)
                      ioutp(1)      = imodu
                      ioutp(2)      = kelty
                      routp(1)      = error
                      routp(2:mm+1) = xx(1:mm)
                      call outfor(97_ip,INT_LIST=lperm)
                   end if
                end if
             end if
          end do
          deallocate( numel     )
          deallocate( numel_gat )
          deallocate( time_gat  )
          if( INOTSLAVE ) then
             deallocate( aa )
             deallocate( bb )
             deallocate( xx )
          end if
       end if
       deallocate( lperm     )
    end if

    !----------------------------------------------------------------------
    !
    ! Parallel efficiency
    !
    !----------------------------------------------------------------------

    if( IPARALL ) then

       call alya2talp_MonitoringRegionStop(GLOBAL_REGION=.true.)
       call alya2talp_parallel_efficiency(&
            PE,LB,CE,&
            time_comp_ave(0),time_mpi_ave(0),&
            time_comp_max(0),time_mpi_max(0),&
            GLOBAL_REGION=.true.) 
       routp(1)    = LB
       routp(2)    = CE
       routp(3)    = PE
       routp(4)    = time_comp_ave(0)
       routp(5)    = time_comp_max(0)
       routp(6)    = time_mpi_ave(0)
       routp(7)    = time_mpi_max(0)
       lb_eff(1,0) = LB
       lb_eff(2,0) = CE
       lb_eff(3,0) = PE
       call outfor(104_ip)
       do imodu = 1,mmodu-1
          if( kfl_modul(imodu) /= 0 ) then
             call alya2talp_parallel_efficiency(&
                  PE,LB,CE,&
                  time_comp_ave(imodu),time_mpi_ave(imodu),&
                  time_comp_max(imodu),time_mpi_max(imodu),&
                  MODULE_REGION=.true.,CURRENT_MODULE=imodu)
             ioutp(1)        = imodu
             routp(1)        = LB
             routp(2)        = CE
             routp(3)        = PE
             routp(4)        = time_comp_ave(imodu)
             routp(5)        = time_comp_max(imodu)
             routp(6)        = time_mpi_ave(imodu)
             routp(7)        = time_mpi_max(imodu)
             lb_eff(1,imodu) = LB
             lb_eff(2,imodu) = CE
             lb_eff(3,imodu) = PE
             call outfor(105_ip)
          end if
       end do

    end if

    !----------------------------------------------------------------------
    !
    ! Compute maximum times over slaves for assembly, solver and output
    !
    !----------------------------------------------------------------------

    do imodu = 1,mmodu
       if( kfl_modul(imodu) /= 0 ) then
          call timings_output(&
               imodu,&    
               cpu_max_element(imodu),&  
               cpu_ave_element(imodu),&  
               lb_ave_element(imodu),&      
               cpu_ave_per_element(imodu), &
               cpu_max_min_element(imodu),&  
               cpu_ave_min_element(imodu),&  
               lb_ave_min_element(imodu),&
               cpu_delta_per(imodu),&
               cpu_max_boundary(imodu),&  
               cpu_ave_boundary(imodu),&  
               lb_ave_boundary(imodu),&    
               cpu_ave_per_boundary(imodu), &
               cpu_max_node(imodu),&  
               cpu_ave_node(imodu),&  
               lb_ave_node(imodu),&    
               cpu_ave_per_node(imodu), &
               cpu_max_particle(imodu),&  
               cpu_ave_particle(imodu),&  
               lb_ave_particle(imodu),&    
               cpu_max_solver(imodu),&
               cpu_ave_solver(imodu),&
               lb_ave_solver(imodu),&
               cpu_max_post(imodu),&
               cpu_ave_post(imodu),&
               lb_ave_post(imodu))        
       end if
    end do
    !
    ! Starting operations
    !
    call PAR_MAX(9_ip,cpu_start,'IN MY CODE','INCLUDE MASTER')
    call cputim(cpu_refer)
    cpu_total =   cpu_refer - cpu_initi
    call PAR_MAX(cpu_total,'IN MY CODE','INCLUDE MASTER')

    !----------------------------------------------------------------------
    !
    ! Write CPU times
    !
    !----------------------------------------------------------------------
    !
    ! Initializations
    !
    cpu_modut = 0.0_rp
    ipass     = 0
    !
    ! Total CPU and CPU for starting operations
    !

    cpu_total = cpu_total + zeror
    routp( 1) =   cpu_total 

    routp( 2) =   cpu_start(CPU_READ_GEO)            + cpu_start(CPU_READ_SETS)           &
         &      + cpu_start(CPU_READ_BCS)            + cpu_start(CPU_READ_FIELDS)         &
         &      + cpu_start(CPU_MESH_PARTITION)      + cpu_start(CPU_MESH_MULTIPLICATION) &
         &      + cpu_start(CPU_CONSTRUCT_DOMAIN)    + cpu_start(CPU_ADDTIONAL_ARRAYS)
    routp( 3) = 100.0_rp * routp(2) / cpu_total

    routp( 4) = cpu_start(CPU_READ_GEO)
    routp( 5) = 100.0_rp * routp( 4) / cpu_total

    routp( 6) = cpu_start(CPU_READ_SETS)
    routp( 7) = 100.0_rp * routp( 6) / cpu_total

    routp( 8) = cpu_start(CPU_READ_BCS)
    routp( 9) = 100.0_rp * routp( 8) / cpu_total

    routp(20) = cpu_start(CPU_READ_FIELDS)
    routp(21) = 100.0_rp * routp(20) / cpu_total

    routp(10) = cpu_start(CPU_MESH_PARTITION)
    routp(11) = 100.0_rp * routp(10) / cpu_total

    routp(14) = cpu_start(CPU_MESH_MULTIPLICATION)
    routp(15) = 100.0_rp * routp(14) / cpu_total

    routp(16) = cpu_start(CPU_CONSTRUCT_DOMAIN)
    routp(17) = 100.0_rp * routp(16) / cpu_total

    routp(18) = cpu_start(CPU_ADDTIONAL_ARRAYS)
    routp(19) = 100.0_rp * routp(18) / cpu_total

    call outfor( 18_ip,lun_outpu,' ')
    !
    ! Module times
    !
    do imodu = 1,mmodu
       if( kfl_modul(imodu) /= 0 ) then
          if(  cpu_modul(CPU_COUNT_ASSEMBLY,imodu) > 0.5_rp .or. &
               cpu_modul(CPU_COUNT_BOUNDARY,imodu) > 0.5_rp .or. &
               cpu_modul(CPU_COUNT_NODE,imodu)     > 0.5_rp .or. &
               cpu_modul(CPU_COUNT_PARTICLE,imodu) > 0.5_rp ) then

             cpu_modut  = cpu_modul(CPU_TOTAL_MODULE,imodu)
             !&cpu_max_element(imodu) + cpu_max_boundary(imodu) &
             !     &     + cpu_max_node(imodu)    + cpu_max_particle(imodu) &
             !     &     + cpu_max_solver(imodu)  + cpu_max_post(imodu)
             cpu_minim  = 1.0e-6_rp
             cpu_denom  = max(cpu_modut,cpu_minim)
             coutp(1)   = namod(imodu)

             routp( 1)  = cpu_modut                               ! Total time
             routp( 2)  = 100.0_rp*routp(1)/cpu_total

             call outfor(19_ip)
             !
             ! Element
             !
             if( cpu_modul(CPU_COUNT_ASSEMBLY,imodu) > 0.5_rp ) then
                coutp(1)   = 'Matrix element assembly:'
                routp(1)   = cpu_ave_element(imodu)
                routp(2)   = 100.0_rp*routp(1)/cpu_denom
                routp(3)   = cpu_max_element(imodu)
                routp(4)   = 100.0_rp*routp(3)/cpu_denom
                routp(5)   = lb_ave_element(imodu)

                routp(6)   = cpu_ave_min_element(imodu)
                routp(7)   = cpu_max_min_element(imodu)
                routp(8)   = lb_ave_min_element(imodu)
                routp(9)   = cpu_delta_per(imodu)                   ! Max variability percentage              
                routp(10)  = cpu_ave_per_element(imodu)
                call outfor(99_ip)
             end if

             if( cpu_modul(CPU_COUNT_BOUNDARY,imodu) > 0.5_rp ) then
                ioutp(1)   = 1
                coutp(1)   = 'Matrix boundary assembly:'
                routp(1)   = cpu_ave_boundary(imodu)
                routp(2)   = 100.0_rp*routp(1)/cpu_denom
                routp(3)   = cpu_max_boundary(imodu)
                routp(4)   = 100.0_rp*routp(3)/cpu_denom
                routp(5)   = lb_ave_boundary(imodu)
                routp(10)  = cpu_ave_per_boundary(imodu)
                call outfor(100_ip)
             end if

             if( cpu_modul(CPU_COUNT_NODE,imodu) > 0.5_rp ) then
                ioutp(1)   = 1
                coutp(1)   = 'Node assembly:'
                routp(1)   = cpu_ave_node(imodu)
                routp(2)   = 100.0_rp*routp(1)/cpu_denom
                routp(3)   = cpu_max_node(imodu)
                routp(4)   = 100.0_rp*routp(3)/cpu_denom
                routp(5)   = lb_ave_node(imodu)
                routp(10)  = cpu_ave_per_node(imodu)
                call outfor(100_ip)
             end if

             if( cpu_modul(CPU_COUNT_PARTICLE,imodu) > 0.5_rp ) then
                ioutp(1)   = 0
                coutp(1)   = 'Particle:'
                routp(1)   = cpu_ave_particle(imodu)
                routp(2)   = 100.0_rp*routp(1)/cpu_denom
                routp(3)   = cpu_max_particle(imodu)
                routp(4)   = 100.0_rp*routp(3)/cpu_denom
                routp(5)   = lb_ave_particle(imodu)
                call outfor(100_ip)
             end if

             ioutp(1)   = 0
             coutp(1)   = 'Solver:'
             routp(1)   = cpu_ave_solver(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             routp(3)   = cpu_max_solver(imodu)
             routp(4)   = 100.0_rp*routp(3)/cpu_denom
             routp(5)   = lb_ave_solver(imodu)
             call outfor(100_ip)

             ioutp(1)   = 0
             coutp(1)   = 'Output'
             routp(1)   = cpu_ave_post(imodu)
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             routp(3)   = cpu_max_post(imodu)
             routp(4)   = 100.0_rp*routp(3)/cpu_denom
             routp(5)   = lb_ave_post(imodu)
             call outfor(100_ip)

             dummr = cpu_modul(CPU_BEGSTE,imodu)
             call PAR_MAX(dummr)
             ioutp(1)   = 2
             coutp(1)   = 'Begin time step'
             routp(1)   = dummr
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)

             dummr = cpu_modul(CPU_ENDSTE,imodu)
             call PAR_MAX(dummr)
             ioutp(1)   = 2
             coutp(1)   = 'End time step'
             routp(1)   = dummr
             routp(2)   = 100.0_rp*routp(1)/cpu_denom
             call outfor(100_ip)           
             !
             ! Write some reporting info
             !
             if( cpu_modul(CPU_COUNT_ASSEMBLY,imodu) > 0.5_rp ) then
                call messages_live('YOU HAVE A HIGH LOAD IMBALANCE IN YOUR MATRIX CONSTRUCTION IN MODULE '&
                     //trim(namod(imodu)),'REPORT')
             end if

          end if
       end if
    end do

  end subroutine performance_outcpu

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-12-30
  !> @brief   Output memory
  !> @details Output information on memory required
  !>
  !-----------------------------------------------------------------------

  subroutine performance_outmem()

    real(rp)     :: rgiga,rmega,rkilo,rbyte
    integer(8)   :: imodu
    character(6) :: lbyte
    integer(ip)  :: ipass,number_passes
    real(rp)     :: r_memor_dom,r_memor_sol,r_memor_els

    if( npart > 1 ) then
       number_passes = 2
    else
       number_passes = 1
    end if
    !
    ! First pass is for Master's max memory
    ! Second pass computes the max memory over the slaves
    !
    do ipass = 1,number_passes

       ioutp(50) = ipass     
       !
       ! Main memory: domain+master+solver
       !
       r_memor_dom = real(memor_dom(2),rp) 
       r_memor_sol = real(memma(2) + memdi(2) + memit(2),rp)     
       r_memor_els = relse(3)
       r_tomax     = r_memor_dom + r_memor_sol + r_memor_els
       max_mem     = mem_maxim    
       ave_mem     = mem_maxim    
       do imodu = 1,mmodu
          if( kfl_modul(imodu) /= 0 ) then
             r_mem_max_modul(imodu)     = real(mem_modul(2,imodu),rp)
             r_mem_ave_modul(imodu) = real(mem_modul(2,imodu),rp)
             r_tomax                = r_tomax + r_mem_max_modul(imodu)
          end if
       end do
       !
       ! Max values
       !
       if( ipass == 2 ) then 
          call PAR_MAX    (r_tomax     , 'IN MY CODE')
          call PAR_MAX    (r_memor_dom , 'IN MY CODE')
          call PAR_MAX    (r_memor_sol , 'IN MY CODE')
          call PAR_MAX    (r_memor_els , 'IN MY CODE')
          call PAR_MAX    (max_mem     , 'IN MY CODE')
          call PAR_AVERAGE(ave_mem     , 'IN MY CODE')
          do imodu = 1,mmodu
             if( kfl_modul(imodu) /= 0 ) then
                call PAR_MAX    (r_mem_max_modul(imodu)    ,'IN MY CODE')
                call PAR_AVERAGE(r_mem_ave_modul(imodu),'IN MY CODE')
             end if
          end do
       end if
       !
       ! Gbutes, Mbytes or bytes?
       !
       rgiga = 1024_rp*1024_rp*1024_rp
       rmega = 1024_rp*1024_rp
       rkilo = 1024_rp     
       if( r_tomax >= rgiga ) then
          rbyte = rgiga
          lbyte = 'Gbytes'
       else if( r_tomax >= rmega ) then 
          rbyte = rmega
          lbyte = 'Mbytes' 
       else if( r_tomax >= rkilo ) then 
          rbyte = rkilo
          lbyte = 'kbytes'          
       else  
          rbyte = 1.0_rp
          lbyte = ' bytes'     
       end if

       routp(1) = r_tomax     / rbyte
       coutp(1) = lbyte
       routp(2) = r_memor_dom / rbyte
       coutp(2) = lbyte
       routp(3) = r_memor_els / rbyte
       coutp(3) = lbyte

       if( INOTSLAVE ) call outfor(21_ip,lun_outpu,' ')
       !
       ! Memory depending on the module
       !
       do imodu = 1,mmodu
          if( kfl_modul(imodu) /= 0 ) then
             coutp(1) = trim(namod(imodu))
             routp(1) = r_mem_max_modul(imodu) / rbyte
             coutp(2) = lbyte
             if( INOTSLAVE ) call outfor(22_ip,lun_outpu,' ')
          end if
       end do
       !
       ! Solver and maximum memory
       !      
       routp(1) = r_memor_sol / rbyte
       coutp(1) = lbyte
       if( INOTSLAVE ) call outfor(24_ip,lun_outpu,' ')

    end do

    !----------------------------------------------------------------------
    !
    ! Variable memory counter
    !
    !----------------------------------------------------------------------

    if( kfl_varcount == 1 ) then
       call memory_output_variable_counter(lun_varcount,OUTPUT_FORMAT="('VARIABLE= ',a,'CALL= ',a,' MEMORY= ',e13.6)")
    end if

  end subroutine performance_outmem

end module mod_performance
!> @}
