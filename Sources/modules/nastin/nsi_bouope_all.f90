!-----------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    nsi_bouope_all.f90
!> @author  Guillaume Houzeaux
!> @brief   Matrix assembly: boundary contribution
!> @details Boundary operations
!> @}
!-----------------------------------------------------------------------
subroutine nsi_bouope_all(itask)
  
  use def_parame
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use def_nastin 
  use mod_ker_proper
  use mod_windk,                        only : mod_windk_interact
  use mod_parall,                       only : par_omp_nboun_chunk
  use mod_parall,                       only : num_subd_nboun_par
  use mod_parall,                       only : num_pack_nboun_par
  use mod_parall,                       only : list_boundaries_par
  use mod_parall,                       only : typ_list_elements_par
  use mod_nsi_boundary_operations,      only : nsi_boundary_operations
  use mod_nsi_boundary_operations_fast, only : nsi_boundary_operations_fast
  use mod_nsi_boundary_operations_fast, only : nsi_boundary_operations_fast5
  use mod_nsi_boundary_operations_fast, only : nsi_boundary_operations_fast_bck1
  use mod_nsi_boundary_operations_fast, only : nsi_boundary_operations_fast_bck2
  use mod_memory
  implicit none 
  integer(ip), intent(in)              :: itask
  integer(ip)                          :: isubd,ipack,iboun
  integer(ip)                          :: pnodb,pgaub,jsubd
  integer(ip)                          :: pblty,ivect
  integer(ip)                          :: num_neigh
  integer(ip)                          :: num_subd
  integer(ip),                 pointer :: num_pack(:)
  type(typ_list_elements_par), pointer :: list_boundaries(:)
  integer(ip)                          :: ifunc
  real(rp)                             :: p_return
  !
  ! Loop indices
  !
  num_subd        =  num_subd_nboun_par
  num_pack        => num_pack_nboun_par
  list_boundaries => list_boundaries_par

  if( IMASTER ) return
  !
  ! Initialize nodal traction and boundary mass
  !
  if( itask == 1_ip ) then
     if( associated(notra_nsi) ) notra_nsi = 0.0_rp
     if( associated(massb_nsi) ) massb_nsi = 0.0_rp
  end if
  !
  ! Loop over boundaries
  !
  do isubd = 1,num_subd

#ifndef OPENACCHHH
#ifdef ALYA_OMPSS
     num_neigh = size(ompss_boundaries(isubd) % neighbours,KIND=ip)

     !-----------------------------------------------------------------------------------
     !$OMP TASK         COMMUTATIVE(                                                    &
     !$OMP              [ompss_boundaries(ompss_boundaries(isubd) % neighbours(jsubd)), &
     !$OMP              jsubd = 1,num_neigh] ) PRIORITY(num_neigh)                      &
     !$OMP FIRSTPRIVATE ( num_neigh,jsubd,isubd )                                       &
     !$OMP SHARED       ( ompss_boundaries )                                            &
     !-----------------------------------------------------------------------------------
#else 
     !-----------------------------------------------------------------------------------
     !$OMP PARALLEL DO                                                                  &
     !$OMP SCHEDULE     ( DYNAMIC , par_omp_nboun_chunk )                               &
     !$OMP SHARED       ( isubd,par_omp_nboun_chunk )                                   &
     !-----------------------------------------------------------------------------------
#endif
     !-----------------------------------------------------------------------------------
     !$OMP DEFAULT      ( SHARED )                                                      &
     !$OMP PRIVATE      ( ipack,pnodb,pgaub,iboun,pblty,ivect )                         &
     !$OMP SHARED       ( list_boundaries,num_pack,ltypb,lnnob,ngaus,lelbo)             &
     !-----------------------------------------------------------------------------------
     !$OMP SHARED       ( itask                                                         )
     !-----------------------------------------------------------------------------------
#endif

     do ipack = 1,num_pack(isubd)

        iboun = list_boundaries(isubd) % packs(ipack) % l(1)                     ! Select first element
        pnodb = lnnob(iboun)                                                     ! Number of boundary nodes
        pblty = ltypb(iboun)                                                     ! Type of boundary
        pgaub = ngaus(pblty)                                                     ! Number of Gauss points

        if ( kfl_asbou_nsi == 5 ) then
           call nsi_boundary_operations_fast5(&
                itask,size(list_boundaries(isubd) % packs(ipack) % l,KIND=ip),pnodb,&
                pgaub,list_boundaries(isubd) % packs(ipack) % l)           
        else
           do ivect = 1,size(list_boundaries(isubd) % packs(ipack) % l,KIND=ip)
              iboun = list_boundaries(isubd) % packs(ipack) % l(ivect)           ! Boundary
              if( iboun > 0 ) &
                   call nsi_boundary_operations(&
                   itask,pnodb,pgaub,list_boundaries(isubd) % packs(ipack) % l(ivect:))
           end do
        end if
        
     end do

#ifndef OPENACCHHH
#ifdef ALYA_OMPSS
     !$OMP END TASK
#else
     !$OMP END PARALLEL DO
#endif
#endif     
  end do
#ifndef OPENACCHHH  
#ifdef ALYA_OMPSS
  !$OMP  TASKWAIT
#endif
#endif
  
end subroutine nsi_bouope_all
