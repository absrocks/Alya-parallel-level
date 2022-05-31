!------------------------------------------------------------------------
!> @addtogroup ChemicMatrixAssembly
!> @{
!> @file    chm_elmope_all.f90
!> @author  Guillaume Houzeaux
!> @brief   Scalar transport element assembly and other element 
!>          calculations
!> @details Elemental operations, to include OpenMP, OmpSs, 
!>          vectorization and CUDA
!>
!> @} 
!------------------------------------------------------------------------
subroutine chm_elmope_all(order)
  
  use def_kintyp,                 only : ip,rp
  use def_domain,                 only : ompss_domain   ! Required by OMPSS
  use def_domain,                 only : ompss_domains
  use def_domain,                 only : lnnod
  use def_domain,                 only : lgaus
  use mod_parall,                 only : par_omp_nelem_chunk
  use mod_parall,                 only : num_subd_par
  use mod_parall,                 only : num_pack_par
  use mod_parall,                 only : list_elements_norace_par
  use mod_parall,                 only : num_subd_norace_par
  use mod_parall,                 only : num_pack_norace_par
  use mod_parall,                 only : list_elements_par
  use mod_parall,                 only : typ_list_elements_par
  use mod_communications,         only : PAR_BARRIER
  use mod_chm_element_operations, only : chm_element_operations
  use mod_chm_finiteRate,         only : chm_element_operations_finiteRate
  use mod_chm_finiteRate,         only : chm_calc_div_enthalpy_transport_finiteRate
  use mod_chm_finiteRate,         only : chm_calc_hk_grad_Yk_others
  use def_chemic,                 only : kfl_model_chm, imixf_rk
  use mod_chm_operations_CMC,     only : chm_element_operations_CMC

  implicit none

  integer(ip), intent(in)              :: order                     !< What to do
  integer(ip)                          :: isubd,ipack,ielem
  integer(ip)                          :: pnode,pgaus,jsubd
  integer(ip)                          :: num_neigh
  integer(ip)                          :: num_subd    
  integer(ip),                 pointer :: num_pack(:)
  type(typ_list_elements_par), pointer :: list_elements(:)

  if( order == 1 .or. order == 6 ) then
     !
     ! Element assembly (ORDER=1) and Laplacian assembly (ORDER=6)
     !
     num_subd      =  num_subd_par
     num_pack      => num_pack_par
     list_elements => list_elements_par

  else if( order == 4 ) then
     !
     ! Subgrid scale update
     !
     num_subd      =  num_subd_norace_par
     num_pack      => num_pack_norace_par
     list_elements => list_elements_norace_par 

  else

     call runend('CHM_ELMOPE_ALL: NOT CODED')

  end if

  !
  ! Compute gradient of species Yk at gauss points projected to the nodes
  !
  if (kfl_model_chm == 3_ip ) & 
      call chm_calc_hk_grad_Yk_others

  if( order == 4 ) then

     !-------------------------------------------------------------------
     !
     ! Subgrid scale: no race condition
     !
     !-------------------------------------------------------------------

     do isubd = 1,num_subd  

        !--------------------------------------------------------------------------     
        !$OMP PARALLEL DO                                                         &
        !$OMP SCHEDULE     ( DYNAMIC , par_omp_nelem_chunk )                      & 
        !$OMP SHARED       ( isubd, par_omp_nelem_chunk )                         &
        !--------------------------------------------------------------------------     
        !--------------------------------------------------------------------------     
        !$OMP DEFAULT      ( NONE )                                               &
        !$OMP PRIVATE      ( ipack,pnode,pgaus,ielem )                            &
        !$OMP SHARED       ( list_elements,num_pack,lnnod,lgaus,kfl_model_chm )   &
        !--------------------------------------------------------------------------     
        !$OMP SHARED       ( order, imixf_rk                                      )                                                                
        !--------------------------------------------------------------------------          

        do ipack = 1,num_pack(isubd)   !pack = vector

           ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
           pnode = lnnod(ielem)                                ! Number of nodes
           pgaus = lgaus(ielem)                                ! Number of Gauss points

           if( kfl_model_chm == 1 ) then
              !
              ! Flamelet model 
              !
              call chm_element_operations(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)

           else if ( kfl_model_chm == 3 ) then
              !
              ! Finite rate kinetics
              !
              call chm_element_operations_finiteRate(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)

           else if ( kfl_model_chm == 4 ) then
              !
              ! CMC model
              !
              call chm_element_operations_CMC(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l, imixf_rk)

           end if

        end do

        !$OMP END PARALLEL DO
     end do

  else

     !-------------------------------------------------------------------
     !
     ! Element assembly and projections: race condition
     !
     !-------------------------------------------------------------------

     do isubd = 1,num_subd
        
#ifdef ALYA_OMPSS
        num_neigh = size(ompss_domains(isubd) % neighbours)

        !-----------------------------------------------------------------------------    
        !$OMP TASK         COMMUTATIVE(                                              &
        !$OMP              [ompss_domains(ompss_domains(isubd) % neighbours(jsubd)), &
        !$OMP              jsubd = 1,num_neigh] ) PRIORITY(num_neigh)                &
        !$OMP FIRSTPRIVATE ( num_neigh,jsubd,isubd )                                 &
        !$OMP SHARED       ( ompss_domains )                                         &
        !-----------------------------------------------------------------------------         
#else
        !-----------------------------------------------------------------------------    
        !$OMP PARALLEL DO                                                            &
        !$OMP SCHEDULE     ( DYNAMIC , par_omp_nelem_chunk )                         &
        !$OMP SHARED       ( isubd,par_omp_nelem_chunk )                             &
        !-----------------------------------------------------------------------------    
#endif
        !-----------------------------------------------------------------------------    
        !$OMP DEFAULT      ( SHARED )                                                &
        !$OMP PRIVATE      ( ipack,pnode,pgaus,ielem )                               &
        !$OMP SHARED       ( list_elements,num_pack,lnnod,lgaus,kfl_model_chm )      &
        !-----------------------------------------------------------------------------    
        !$OMP SHARED       ( order, imixf_rk                                         )               
        !-----------------------------------------------------------------------------         

        do ipack = 1,num_pack(isubd)

           ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
           pnode = lnnod(ielem)                                ! Number of nodes
           pgaus = lgaus(ielem)                                ! Number of Gauss points

           if( kfl_model_chm == 1 ) then
              !
              ! Flamelet model 
              !
              call chm_element_operations(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)

           else if ( kfl_model_chm == 3 ) then
              !
              ! Finite rate kinetics
              !
              call chm_element_operations_finiteRate(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l)

           else if ( kfl_model_chm == 4 ) then
              !
              ! CMC model
              !
              call chm_element_operations_CMC(&
                               order,pnode,pgaus,list_elements(isubd) % packs(ipack) % l, imixf_rk)

           end if

        end do

#ifdef ALYA_OMPSS
        !$OMP END TASK
#else
        !$OMP END PARALLEL DO
#endif
     end do
#ifdef ALYA_OMPSS
     !$OMP  TASKWAIT
#endif
  end if

  !
  ! Compute divergence of enthalpy transport
  !
  if (kfl_model_chm == 3_ip ) & 
      call chm_calc_div_enthalpy_transport_finiteRate


end subroutine chm_elmope_all
