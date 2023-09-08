!------------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_boundary_loop.f90
!> @author  Guillaume Houzeaux
!> @date    10/10/2006
!> @brief   Define arrays for hybrid parallelization
!> @details Define arrays for hybrid parallelization: includes OpenMP, 
!>          OmpSs, vectorization and CUDA.
!>          For element loops, elements are ordered according to
!>          pnode, pgaus, pmate
!>
!>          1. Without race condition:
!>          --------------------------
!>
!>             1.1 Classical element loop
!>             1.2 Otherwise:
!>                 do isubd = 1,num_subd_norace_par
!>                   do ipack = 1,num_pack_norace_par(isubd)
!>                     list_elements_norace_par
!>
!>          2. With race condition
!>          ----------------------
!>             Choose a hybrid method with flag 
!>             par_hybrid = PAR_OPENMP_COLORING
!>             par_hybrid = PAR_OPENMP_NO_COLORING
!>             par_hybrid = PAR_OMPSS
!>
!>             2.1 PAR_OPENMP_NO_COLORING:
!>                 define NO_COLORING macro and do not forget to put
!>                 ATOMIC pragma or define CRITICAL sections
!>                 using #ifdef NO_COLORING. Otherwise, DO NOT define
!>                 this macro.
!>             2.2 PAR_OMPSS:
!>                 Define ALYA_OMPSS macro
!>                 DO NOT define NO_COLORING macro!
!>             2.3 PAR_OPENMP_COLORING:
!>                 Do not do anything...
!>                 
!>          If no vectorization is used, add the following inner loop:
!>          do kelem = 1,VECTOR_SIZE
!>             iboun = list_elements_par(isubd) % packs(ipack) % l(kelem)
!>             if( iboun > 0 ) call element_assembly(iboun)
!>
!>
!> @} 
!------------------------------------------------------------------------
subroutine par_boundary_loop()

  use def_kintyp,         only : ip,rp,lg,i1p
#ifndef VECTOR_SIZE
  use def_master,         only : VECTOR_SIZE
#endif

  use def_master,         only : INOTMASTER,kfl_paral,intost
  use def_domain,         only : nboun,lnnob,lmate,ngaus,lelbo
  use def_domain,         only : mgaub,mnodb,nmate,ltypb,lnnod
  use def_domain,         only : ompss_boundaries
  use mod_parall,         only : par_memor
  use mod_parall,         only : num_subd_nboun_par
  use mod_parall,         only : num_pack_nboun_par
  use mod_parall,         only : list_boundaries_par
  use mod_parall,         only : num_subd_norace_nboun_par
  use mod_parall,         only : num_pack_norace_nboun_par
  use mod_parall,         only : list_boundaries_norace_par
  use mod_parall,         only : par_omp_num_threads
  use mod_parall,         only : par_omp_nboun_num_colors
  use mod_parall,         only : par_omp_nboun_ia_colors 
  use mod_parall,         only : par_omp_nboun_ja_colors
  use mod_parall,         only : par_hybrid
  use mod_parall,         only : PAR_OPENMP_COLORING 
  use mod_parall,         only : PAR_OPENMP_NO_COLORING
  use mod_parall,         only : PAR_OMPSS
  use mod_parall,         only : PAR_HYBRID_OFF
  use mod_maths,          only : maths_mapping_1d_to_3d_x
  use mod_maths,          only : maths_mapping_1d_to_3d_y
  use mod_maths,          only : maths_mapping_1d_to_3d_z
  use mod_maths,          only : maths_mapping_3d_to_1d
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_size
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SUM
  use mod_messages,       only : messages_live
  implicit none

  integer(ip)             :: isubd,ipack,itype
  integer(ip)             :: ivect,kelem,iboun
  integer(ip)             :: pnode,pgaub,pmate
  integer(ip)             :: imate,inode,igaus
  integer(ip)             :: jtype,jpack,ielem
  integer(ip)             :: VECTOR_SIZE_LOC
  integer(ip)             :: hybrid_method,pnodb

  integer(ip)             :: num_types
  integer(ip)             :: num_holes
  integer(ip), pointer    :: num_element_types(:)
  logical(lg), pointer    :: consider_element(:) 
  integer(ip), pointer    :: element_types(:)
  logical(lg)             :: adjust_vector_size
  logical(lg)             :: if_vector_size
  
  if_vector_size = .true.
  !
  ! Deallocate memory if necessary
  !
  call memory_deallo(par_memor,'NUM_PACK_NBOUN_PAR','par_element_loop',num_pack_nboun_par)
  if( associated(list_boundaries_par) ) then
     do isubd = 1,size(list_boundaries_par) 
        call memory_deallo(par_memor,'LIST_BOUNDARIES_PAR(ISUBD) % PACKS','par_element_loop',list_boundaries_par(isubd) % packs)
     end do
     deallocate(list_boundaries_par)
  end if
  call memory_deallo(par_memor,'NUM_PACK_NORACE_NBOUN_PAR','par_element_loop',num_pack_norace_nboun_par)
  if( associated(list_boundaries_norace_par) ) then
     do isubd = 1,size(list_boundaries_norace_par) 
        call memory_deallo(par_memor,'LIST_BOUNDARIES_NORACE_PAR(ISUBD) % PACKS','par_element_loop',list_boundaries_norace_par(isubd) % packs)
     end do
     deallocate(list_boundaries_norace_par)
  end if
  !
  ! Errors
  !
  if( par_hybrid == PAR_HYBRID_OFF .and. par_omp_num_threads /= 0 ) then
     call runend('PARALL: CHOOSE A HYBRID PARALLELIZATION STRATEGY')
  end if
  if( par_hybrid /= PAR_HYBRID_OFF .and. par_omp_num_threads == 0 ) then
     call runend('PARALL: SET ENVIRONMENT VARIABLE OMP_NUM_THREAD')
  end if
  !
  ! Vector size
  !
  adjust_vector_size = .false.
  if( if_vector_size ) then
     if( VECTOR_SIZE == -1 ) then
        adjust_vector_size = .true.
     else
        VECTOR_SIZE_LOC = max(1_ip,int(VECTOR_SIZE,ip))
     end if
  else
     VECTOR_SIZE_LOC = 1
  end if

  call messages_live('HYBRID PARALLELIZATION & VECTORIZATION (BOUNDARIES)','START SECTION')
  if( if_vector_size ) then
     if( VECTOR_SIZE > 1 ) then
        call messages_live('VECTORIZATION WITH VECTOR SIZE= '//trim(intost(int(VECTOR_SIZE_LOC,ip))))
     else if( VECTOR_SIZE == -1 ) then
        call messages_live('VECTORIZATION AUTOMATICALLY ADJUSTED')
     else     
        call messages_live('NO VECTORIZATION')
     end if
  end if

  num_types = 0
  num_holes = 0
  jtype     = 0
  
  if( INOTMASTER ) then
     !
     ! Initialization
     !
     num_subd_nboun_par = 0
     num_subd_norace_nboun_par = 0
     nullify(num_pack_nboun_par)
     nullify(list_boundaries_par)
     nullify(num_pack_norace_nboun_par)
     nullify(list_boundaries_norace_par)
     !
     ! Local variables
     !
     nullify(num_element_types)
     nullify(consider_element)
     nullify(element_types)
     !
     ! Number and list of element types: NUM_ELEMENT_TYPES, ELEMENT_TYPES
     !
     if( VECTOR_SIZE_LOC == 1 .and. .not. adjust_vector_size ) then
        num_types = 1
        num_holes = 0
     else
        num_types = nmate * mnodb * mgaub
        num_holes = 0
     end if

     call memory_alloca(par_memor,'NUM_ELEMENT_TYPES','par_element_loop',num_element_types,num_types)
     call memory_alloca(par_memor,'ELEMENT_TYPES'    ,'par_element_loop',element_types,nboun)
     
     if( VECTOR_SIZE_LOC == 1 ) then
        num_element_types(1) = 0
        do iboun = 1,nboun
           if( ltypb(iboun) > 0 ) then
              element_types(iboun) = 1
              num_element_types(1) = num_element_types(1) + 1
           else
              num_holes = num_holes + 1
           end if
        end do
        jtype = 1
     else
        !
        ! This could fail because the same boudnary type could belong to different element types
        !
        do iboun = 1,nboun
           if( ltypb(iboun) > 0 ) then
              ielem                    = lelbo(iboun)
              pmate                    = lmate(ielem)
              pnode                    = lnnod(ielem)
              pnodb                    = lnnob(iboun)
              pgaub                    = ngaus(ltypb(iboun))
              itype                    = maths_mapping_3d_to_1d(nmate,mnodb,mgaub,pmate,pnodb,pgaub)
              element_types(iboun)     = itype
              num_element_types(itype) = num_element_types(itype) + 1
           else
              num_holes = num_holes + 1
           end if
        end do
        jtype = 0
        do itype = 1,num_types
           jtype = jtype + min(1_ip,num_element_types(itype))
        end do
     end if

  end if

  call PAR_MAX(jtype)
  call PAR_SUM(num_holes)
  call messages_live('MAX NUMBER OF BOUNDARY TYPES DETECTED= '//trim(intost(jtype)))
  if( num_holes > 0 ) then
     call messages_live('NUMBER OF HOLE BOUNDARIES= '//trim(intost(num_holes)))     
  end if

  !----------------------------------------------------------------------
  !
  ! Loops without race condition: classical OpenMP
  !
  ! NUM_SUBD_NORACE_NBOUN_PAR = 1
  ! NUM_PACK_NORACE_NBOUN_PAR
  ! LIST_BOUNDARIES_NORACE_PAR
  !
  !----------------------------------------------------------------------
  
  if( par_hybrid /= PAR_HYBRID_OFF ) then
     call messages_live('NUMBER OF THREADS= '//trim(intost(par_omp_num_threads)))
     call messages_live('LOOP WITHOUT RACE CONDITION: OPENMP WITHOUT COLORING')
  else
     call messages_live('OPENMP/OMPSS DESACTIVATED IN ALL LOOPS')     
  end if
  
  if( INOTMASTER ) then

     num_subd_norace_nboun_par = 1
     allocate( list_boundaries_norace_par(num_subd_norace_nboun_par) )
     do isubd = 1,num_subd_norace_nboun_par
        nullify(list_boundaries_norace_par(isubd) % packs)
     end do
     call memory_alloca(par_memor,'NUM_PACK_NORACE_NBOUN_PAR','par_element_loop',num_pack_norace_nboun_par,num_subd_norace_nboun_par)

     isubd = 1
     num_element_types = 0
     do iboun = 1,nboun
        if( ltypb(iboun) > 0 ) then
           itype = element_types(iboun)
           num_element_types(itype) = num_element_types(itype) + 1
        end if
     end do

     if( adjust_vector_size ) then
        num_pack_norace_nboun_par(isubd) = 0 
        do itype = 1,num_types
           if( num_element_types(itype) > 0 ) then
              num_pack_norace_nboun_par(isubd) = num_pack_norace_nboun_par(isubd) + 1
           end if
        end do
     else
        num_pack_norace_nboun_par(isubd) = 0 
        do itype = 1,num_types
           if( num_element_types(itype) > 0 ) then
              num_pack_norace_nboun_par(isubd) = num_pack_norace_nboun_par(isubd) + ceiling(real(num_element_types(itype),rp) / real(VECTOR_SIZE_LOC,rp), ip )
           end if
        end do
     end if     
     call memory_alloca(par_memor,'LIST_BOUNDARIES_NORACE_PAR','par_element_loop',list_boundaries_norace_par(isubd) % packs,num_pack_norace_nboun_par(isubd))
     
     if( adjust_vector_size ) then
        ipack = 0
        do itype = 1,num_types
           if( num_element_types(itype) > 0 ) then
              ipack = ipack + 1
              call memory_alloca(par_memor,'LIST_BOUNDARIES_NORACE_PAR','par_element_loop',list_boundaries_norace_par(isubd) % packs(ipack) % l,num_element_types(itype))
           end if
        end do
     else
        do ipack = 1,num_pack_norace_nboun_par(isubd)
           call memory_alloca(par_memor,'LIST_BOUNDARIES_NORACE_PAR','par_element_loop',list_boundaries_norace_par(isubd) % packs(ipack) % l,VECTOR_SIZE_LOC) 
        end do
     end if
          
     ivect = 0
     ipack = 0 
     jtype = 0
     do itype = 1,num_types
        if( num_element_types(itype) > 0 ) then
           do iboun = 1,nboun
              if( ltypb(iboun) > 0 .and. element_types(iboun) == itype ) then
                 ivect = ivect + 1
                 if( adjust_vector_size ) then
                    if( ivect > num_element_types(itype) .or. itype /= jtype ) then                    
                       ivect = 1
                       ipack = ipack + 1
                    end if
                 else
                    if( ivect > VECTOR_SIZE_LOC .or. itype /= jtype ) then                    
                       ivect = 1
                       ipack = ipack + 1
                    end if
                 end if
                 if( ipack > size(list_boundaries_norace_par(isubd) % packs) ) &
                      call runend('WRONG PACK NUMBER: CALL GUILLAUME TO FIX IT')
                 list_boundaries_norace_par(isubd) % packs(ipack) % l(ivect) = iboun
                 jtype = itype
              end if
           end do
        end if
     end do
     !
     ! Number of subdomains NUM_SUBD_NBOUN_PAR
     !
     if(      par_hybrid == PAR_HYBRID_OFF ) then
        num_subd_nboun_par = 1
     else if( par_hybrid == PAR_OPENMP_NO_COLORING ) then
        num_subd_nboun_par = 1
     else if( par_hybrid == PAR_OPENMP_COLORING ) then
        num_subd_nboun_par = par_omp_nboun_num_colors
     else if( par_hybrid == PAR_OMPSS ) then
        if( associated(ompss_boundaries) ) num_subd_nboun_par = size(ompss_boundaries)           
     end if

     allocate( list_boundaries_par(num_subd_nboun_par) )
     do isubd = 1,num_subd_nboun_par
        nullify(list_boundaries_par(isubd) % packs)
     end do
     call memory_alloca(par_memor,'NUM_PACK_NBOUN_PAR','par_element_loop',num_pack_nboun_par,num_subd_nboun_par)

  end if

  if( par_hybrid == PAR_HYBRID_OFF .or. par_hybrid == PAR_OPENMP_NO_COLORING ) then

     !-------------------------------------------------------------------
     !
     ! No hybrid parallelization 
     !
     !-------------------------------------------------------------------

     if( par_hybrid /= PAR_HYBRID_OFF ) then
        call messages_live('LOOP WITH RACE CONDITION: OPENMP WITHOUT COLORING')
     end if

     if( INOTMASTER ) then

        isubd = 1
        num_element_types = 0
        do iboun = 1,nboun
           if( ltypb(iboun) > 0 ) then
              itype = element_types(iboun)
              num_element_types(itype) = num_element_types(itype) + 1
           end if
        end do
        
        if( adjust_vector_size ) then
           num_pack_nboun_par(isubd) = 0 
           do itype = 1,num_types
              if( num_element_types(itype) > 0 ) then
                 num_pack_nboun_par(isubd) = num_pack_nboun_par(isubd) + 1
              end if
           end do
        else
           num_pack_nboun_par(isubd) = 0 
           do itype = 1,num_types
              if( num_element_types(itype) > 0 ) then
                 num_pack_nboun_par(isubd) = num_pack_nboun_par(isubd) + ceiling(real(num_element_types(itype),rp) / real(VECTOR_SIZE_LOC,rp), ip )
              end if
           end do
        end if
        
        call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR','par_element_loop',list_boundaries_par(isubd) % packs,num_pack_nboun_par(isubd))
        if( adjust_vector_size ) then
           ipack = 0
           do itype = 1,num_types
              if( num_element_types(itype) > 0 ) then
                 ipack = ipack + 1
                 call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR','par_element_loop',list_boundaries_par(isubd) % packs(ipack) % l,num_element_types(itype))
              end if
           end do
        else
           do ipack = 1,num_pack_nboun_par(isubd)
              call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR','par_element_loop',list_boundaries_par(isubd) % packs(ipack) % l,VECTOR_SIZE_LOC) 
           end do
        end if
        
        ivect = 0
        ipack = 0 
        jtype = 0
        do itype = 1,num_types
           if( num_element_types(itype) > 0 ) then
              do iboun = 1,nboun
                 if( ltypb(iboun) > 0 .and. element_types(iboun) == itype ) then
                    ivect = ivect + 1
                    if( adjust_vector_size ) then
                       if( ivect > num_element_types(itype) .or. itype /= jtype ) then                    
                          ivect = 1
                          ipack = ipack + 1
                       end if
                    else
                       if( ivect > VECTOR_SIZE_LOC .or. itype /= jtype ) then                    
                          ivect = 1
                          ipack = ipack + 1
                       end if
                    end if
                    if( ipack > size(list_boundaries_par(isubd) % packs) ) call runend('WRONG PACK NUMBER')
                    list_boundaries_par(isubd) % packs(ipack) % l(ivect) = iboun
                    jtype = itype
                 end if
              end do
           end if
        end do

     end if

  else if( par_hybrid == PAR_OPENMP_COLORING  ) then

     !-------------------------------------------------------------------
     !
     ! OpenMP with coloring
     !
     !-------------------------------------------------------------------

     call messages_live('LOOP WITH RACE CONDITION: OPENMP WITH COLORING')

     if( INOTMASTER ) then

        do isubd = 1,num_subd_nboun_par
           num_element_types = 0
           do kelem = par_omp_nboun_ia_colors(isubd),par_omp_nboun_ia_colors(isubd+1)-1
              iboun = par_omp_nboun_ja_colors(kelem)
              if( ltypb(iboun) > 0 ) then
                 itype = element_types(iboun)
                 num_element_types(itype) = num_element_types(itype) + 1
              end if
           end do

           if( adjust_vector_size ) then
              num_pack_nboun_par(isubd) = num_types
           else
              num_pack_nboun_par(isubd) = 0
              do itype = 1,num_types
                 if( num_element_types(itype) > 0 ) then
                    num_pack_nboun_par(isubd) = num_pack_nboun_par(isubd) + ceiling(real(num_element_types(itype),rp) / real(VECTOR_SIZE_LOC,rp), ip )
                 end if
              end do
           end if
           
           call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR','par_element_loop',list_boundaries_par(isubd) % packs,num_pack_nboun_par(isubd))
           if( adjust_vector_size ) then
              do ipack = 1,num_pack_nboun_par(isubd)
                 itype = ipack
                 call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR','par_element_loop',list_boundaries_par(isubd) % packs(ipack) % l,num_element_types(itype))
              end do
           else
              do ipack = 1,num_pack_nboun_par(isubd)
                 call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR','par_element_loop',list_boundaries_par(isubd) % packs(ipack) % l,VECTOR_SIZE_LOC) 
              end do
           end if
           ivect = 0
           ipack = 0 
           jtype = 0

           do itype = 1,num_types
              if( num_element_types(itype) > 0 ) then
                 do kelem = par_omp_nboun_ia_colors(isubd),par_omp_nboun_ia_colors(isubd+1)-1
                    iboun = par_omp_nboun_ja_colors(kelem)
                    if( ltypb(iboun) > 0 .and. element_types(iboun) == itype ) then
                       ivect = ivect + 1
                       if( adjust_vector_size ) then
                          if( ivect > num_element_types(itype) .or. itype /= jtype ) then                    
                             ivect = 1
                             ipack = ipack + 1
                          end if
                       else
                          if( ivect > VECTOR_SIZE_LOC .or. itype /= jtype ) then                    
                             ivect = 1
                             ipack = ipack + 1
                          end if
                       end if
                       if( ipack > size(list_boundaries_par(isubd) % packs) ) call runend('WRONG PACK NUMBER')
                       list_boundaries_par(isubd) % packs(ipack) % l(ivect) = iboun
                       jtype = itype
                    end if
                 end do
              end if
           end do
        end do
     end if

  else if( par_hybrid == PAR_OMPSS ) then

     !-------------------------------------------------------------------
     !
     ! OmpSs
     !
     !-------------------------------------------------------------------

     call messages_live('LOOP WITH RACE CONDITION: OMPSS WITH MULDIDEPENDENCIES')

     if( INOTMASTER .and. nboun > 0 ) then
        do isubd = 1,num_subd_nboun_par
           num_element_types = 0
           do kelem = 1,size(ompss_boundaries(isubd) % elements)
              iboun = ompss_boundaries(isubd) % elements(kelem)
              if( ltypb(iboun) > 0 ) then
                 itype = element_types(iboun)
                 num_element_types(itype) = num_element_types(itype) + 1
              end if
           end do
           num_pack_nboun_par(isubd) = 0 
           do itype = 1,num_types
              if( num_element_types(itype) > 0 ) then
                 num_pack_nboun_par(isubd) = num_pack_nboun_par(isubd) + ceiling(real(num_element_types(itype),rp) / real(VECTOR_SIZE_LOC,rp), ip )
              end if
           end do
           call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR','par_element_loop',list_boundaries_par(isubd) % packs,num_pack_nboun_par(isubd))
           do ipack = 1,num_pack_nboun_par(isubd)
              call memory_alloca(par_memor,'LIST_BOUNDARIES_PAR','par_element_loop',list_boundaries_par(isubd) % packs(ipack) % l,VECTOR_SIZE_LOC) 
           end do

           ivect = 0
           ipack = 0 
           jtype = 0
           do itype = 1,num_types
              if( num_element_types(itype) > 0 ) then
                 do kelem = 1,size(ompss_boundaries(isubd) % elements)
                    iboun = ompss_boundaries(isubd) % elements(kelem)
                    if( ltypb(iboun) > 0 .and. element_types(iboun) == itype ) then
                       ivect = ivect + 1
                       if( ivect > VECTOR_SIZE_LOC .or. itype /= jtype ) then                    
                          ivect = 1
                          ipack = ipack + 1
                       end if
                       if( ipack > size(list_boundaries_par(isubd) % packs) ) call runend('WRONG PACK NUMBER')
                       list_boundaries_par(isubd) % packs(ipack) % l(ivect) = iboun
                       jtype = itype
                    end if
                 end do
              end if
           end do
        end do

     end if

  end if
  !
  ! Check
  !
  if( INOTMASTER .and. nboun > 0 ) then

     if( 1 == 1 ) then
        call memory_alloca(par_memor,'CONSIDER_ELEMENT','par_element_loop',consider_element  ,nboun)
        do isubd = 1,num_subd_nboun_par
           do ipack = 1,num_pack_nboun_par(isubd)
              do kelem = 1,memory_size(list_boundaries_par(isubd) % packs(ipack) % l)
                 iboun = list_boundaries_par(isubd) % packs(ipack) % l(kelem)
                 if( iboun > 0 ) then
                    if( ltypb(iboun) > 0 ) then
                       if( consider_element(iboun)) call runend('ELEMENT ALREADY TOUCHED')
                       consider_element(iboun) = .true.
                    end if
                 end if
              end do
           end do
        end do
        kelem = 0
        inode = 0
        do iboun = 1,nboun
           if( .not. consider_element(iboun) .and. ltypb(iboun) > 0 ) kelem=kelem+1
           if( consider_element(iboun) ) inode = inode + 1
        end do
        if( inode /= nboun ) then
           print*,'error=',kfl_paral,nboun,inode
           call runend('WE ARE IN TROUBLE')
        end if
        if(kelem>0) call runend('ERROR: '//trim(intost(kelem))//' ELEMENTS NOT TOUCHED')
        do isubd = 1,num_subd_nboun_par
           do ipack = 1,num_pack_nboun_par(isubd)
              iboun = list_boundaries_par(isubd) % packs(ipack) % l(1)
              itype = element_types(iboun)
              do kelem = 2,memory_size(list_boundaries_par(isubd) % packs(ipack) % l)
                 iboun = list_boundaries_par(isubd) % packs(ipack) % l(kelem)
                 if( iboun > 0 ) then
                    if( element_types(iboun) /= itype ) then
                       print*,'elem 1= ',list_boundaries_par(isubd) % packs(ipack) % l(1),itype
                       print*,'elem 2= ',iboun,element_types(iboun)
                       call runend('ERROR: ELEMENT OF DIFFERENT TYPES IN THE SAME LIST')
                    end if
                 end if
              end do
           end do
        end do
        call memory_deallo(par_memor,'CONSIDER_ELEMENT','par_element_loop',consider_element)
     end if
     !
     ! Deallocate memory
     !
     call memory_deallo(par_memor,'ELEMENT_TYPES'    ,'par_element_loop',element_types)
     call memory_deallo(par_memor,'NUM_ELEMENT_TYPES','par_element_loop',num_element_types)

  end if

  call messages_live('HYBRID PARALLELIZATION & VECTORIZATION (BOUNDARIES)','END SECTION')

end subroutine par_boundary_loop
