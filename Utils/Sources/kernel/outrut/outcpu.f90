!-----------------------------------------------------------------------
!> @addtogroup CPU_Time
!> @{
!> @file    outcpu.f90
!> @author  Guillaume Houzeaux
!> @date    15/07/2015
!> @brief   Output general CPU time info
!> @details Output general CPU time info about starting operation and
!>          modules
!>
!> @}
!-----------------------------------------------------------------------

subroutine outcpu()

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
  use mod_std
  implicit none

  real(rp)             :: cpu_refer,cpu_minim,cpu_denom,error
  real(rp)             :: cpu_modut,cpu_servt,cpu_total
  real(rp)             :: dummr
  integer(ip)          :: imodu,ipass,kelty,ielty,jelty,nn,mm
  real(rp)             :: cpu_max_element(mmodu)  
  real(rp)             :: cpu_ave_element(mmodu)  
  real(rp)             :: lb_ave_element(mmodu)
  real(rp)             :: cpu_ave_per_element(mmodu)
  real(rp)             :: cpu_max_boundary(mmodu)  
  real(rp)             :: cpu_ave_boundary(mmodu)  
  real(rp)             :: lb_ave_boundary(mmodu)    
  real(rp)             :: cpu_ave_per_boundary(mmodu)
  real(rp)             :: cpu_max_node(mmodu)  
  real(rp)             :: cpu_ave_node(mmodu)  
  real(rp)             :: lb_ave_node(mmodu)    
  real(rp)             :: cpu_ave_per_node(mmodu)
  real(rp)             :: cpu_max_particle(mmodu)  
  real(rp)             :: cpu_ave_particle(mmodu)  
  real(rp)             :: lb_ave_particle(mmodu)    
  real(rp)             :: cpu_max_min_element(mmodu)  
  real(rp)             :: cpu_ave_min_element(mmodu)  
  real(rp)             :: lb_ave_min_element(mmodu)
  real(rp)             :: cpu_delta_per(mmodu)
  real(rp)             :: cpu_max_solver(mmodu)
  real(rp)             :: cpu_ave_solver(mmodu)
  real(rp)             :: lb_ave_solver(mmodu)
  real(rp)             :: cpu_max_post(mmodu)
  real(rp)             :: cpu_ave_post(mmodu)
  real(rp)             :: lb_ave_post(mmodu)        
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
                    xx            = xx/xx(1)
                    ioutp(1)      = imodu
                    ioutp(2)      = kelty
                    routp(1)      = error
                    routp(2:mm+1) = xx(1:mm)
                    call outfor(97_ip,ilist=lperm)
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

  if( INOTSLAVE ) then
     !
     ! Initializations
     !

     cpu_modut = 0.0_rp
     cpu_servt = 0.0_rp
     ipass     = 0
     !
     ! Total CPU and CPU for starting operations
     !

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

     call outfor(18_ip,lun_outpu,' ')
     !
     ! Module times
     !
     do imodu = 1,mmodu
        if( kfl_modul(imodu) /= 0 ) then
           if(  cpu_modul(CPU_COUNT_ASSEMBLY,imodu) > 0.5_rp .or. &
                cpu_modul(CPU_COUNT_BOUNDARY,imodu) > 0.5_rp .or. &
                cpu_modul(CPU_COUNT_NODE,imodu)     > 0.5_rp .or. &
                cpu_modul(CPU_COUNT_PARTICLE,imodu) > 0.5_rp ) then

              cpu_modut  = cpu_max_element(imodu) + cpu_max_boundary(imodu) &
                   &     + cpu_max_node(imodu)    + cpu_max_particle(imodu) &
                   + cpu_max_solver(imodu)  + cpu_max_post(imodu)
              cpu_minim  = 1.0e-6_rp
              cpu_denom  = max(cpu_modut,cpu_minim)
              coutp(1)   = namod(imodu)

              routp( 1)   = cpu_modut                              ! Total time
              routp( 2)   = 100.0_rp*routp(1)/cpu_total
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

  end if

end subroutine outcpu
