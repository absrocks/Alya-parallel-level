!------------------------------------------------------------------------
!> @addtogroup NastinMatrixAssembly
!> @{
!> @file    nsi_elmope_all.f90
!> @author  Guillaume Houzeaux
!> @brief   Navier-Stokes system element assembly and other element
!>          calculations
!> @details Elemental operations, to include OpenMP, OmpSs,
!>          vectorization and CUDA
!>
!> @}
!------------------------------------------------------------------------
subroutine nsi_elmope_all(itask)
  use def_kintyp,                 only : ip,rp
  use def_master,                 only : IEMPTY
  use def_domain,                 only : ompss_domain   ! Required by OMPSS
  use def_domain,                 only : ompss_domains
  use def_domain,                 only : lnnod
  use def_domain,                 only : lgaus
  use def_nastin,                 only : resis_nsi
  use def_nastin,                 only : itsta_nsi
  use def_nastin,                 only : resgs_nsi
  use def_nastin,                 only : tamin_nsi
  use def_nastin,                 only : rmsgs_nsi
  use def_nastin,                 only : tamax_nsi
  use def_nastin,                 only : dtmax_nsi
  use def_nastin,                 only : cputi_assembly_nsi
  use def_nastin,                 only : kfl_assem_nsi
  use mod_parall,                 only : par_omp_nelem_chunk
  use mod_parall,                 only : num_subd_par
  use mod_parall,                 only : num_pack_par
  use mod_parall,                 only : list_elements_norace_par
  use mod_parall,                 only : num_subd_norace_par
  use mod_parall,                 only : num_pack_norace_par
  use mod_parall,                 only : list_elements_par
  use mod_parall,                 only : typ_list_elements_par
  use def_master,                 only : ittim,kfl_paral

  use mod_communications,              only : PAR_BARRIER
  use mod_nsi_element_operations,      only : nsi_element_operations
  use mod_nsi_element_operations_fast, only : nsi_element_operations_fast
  use mod_nsi_element_operations_test, only : nsi_element_operations_fast8
  use def_master,                      only : veloc,rhsid
  use def_nastin,                      only : dt_rho_nsi,mass_rho_nsi
  use def_parall,                      only : kfl_streams_per_gpu
  use def_kermod,                      only : turmu_ker 
  use def_master,                      only : amatr,rhsid

#ifdef OPENACCHHH
    use openacc 
#endif

  implicit none

  integer(ip), intent(in)              :: itask                     !< What to do
  integer(ip)                          :: isubd,ipack,ielem
  integer(ip)                          :: pnode,pgaus,jsubd
  integer(ip)                          :: num_neigh
  integer(ip)                          :: num_subd,i
  integer(ip),                 pointer :: num_pack(:)
  type(typ_list_elements_par), pointer :: list_elements(:)
  real(rp)                             :: time_detail(10)
  integer(ip)                          :: sizetot
  integer(ip)                          :: streamid
  integer(ip)                          :: ngpus

  if( IEMPTY ) return
  dtmax_nsi =-1.0_rp
  
  if( itask == 1 .or. itask == 6 ) then
     !
     ! Element assembly (ITASK=1) and Laplacian assembly (ITASK=6)
     !
     num_subd      =  num_subd_par
     num_pack      => num_pack_par
     list_elements => list_elements_par

  else if( itask == 4 ) then
     !
     ! Subgrid scale update
     !
     num_subd      =  num_subd_norace_par
     num_pack      => num_pack_norace_par
     list_elements => list_elements_norace_par

  else

     call runend('NSI_ELMOPE_ALL: NOT CODED')

  end if
  
  time_detail = 0.0_rp

  if( itask == 4 ) then

     !-------------------------------------------------------------------
     !
     ! Subgrid scale: no race condition
     !
     !-------------------------------------------------------------------

     do isubd = 1,num_subd    !colors

        !--------------------------------------------------------------------------
        !$OMP PARALLEL DO                                                         &
        !$OMP SCHEDULE     ( DYNAMIC , par_omp_nelem_chunk )                      &
        !$OMP SHARED       ( isubd, par_omp_nelem_chunk )                         &
        !--------------------------------------------------------------------------
        !--------------------------------------------------------------------------
        !$OMP DEFAULT      ( NONE )                                               &
        !$OMP PRIVATE      ( ipack,pnode,pgaus,ielem )                            &
        !$OMP SHARED       ( list_elements,num_pack,lnnod,lgaus )                 &
        !--------------------------------------------------------------------------
        !$OMP SHARED       ( itask )                                              &
        !$OMP REDUCTION    ( +:resis_nsi,itsta_nsi,resgs_nsi,time_detail )        &
        !$OMP REDUCTION    ( MIN:tamin_nsi )                                      &
        !$OMP REDUCTION    ( MAX:rmsgs_nsi,tamax_nsi,dtmax_nsi                    )
        !--------------------------------------------------------------------------

        do ipack = 1,num_pack(isubd)   !pack = vector

           ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
           pnode = lnnod(ielem)                                ! Number of nodes
           pgaus = lgaus(ielem)                                ! Number of Gauss points

           call nsi_element_operations(&
                itask,pnode,pgaus,list_elements(isubd) % packs(ipack) % l,&
                resis_nsi,itsta_nsi,resgs_nsi,tamin_nsi,rmsgs_nsi,tamax_nsi,&
                dtmax_nsi,time_detail)

        end do

        !$OMP END PARALLEL DO
     end do

  else !itask /= 4

     !-------------------------------------------------------------------
     !
     ! Element assembly and projections: race condition
     !
     !-------------------------------------------------------------------

     !$acc enter data create(dt_rho_nsi,mass_rho_nsi,rhsid,veloc)

     !$acc update device(dt_rho_nsi,mass_rho_nsi,rhsid,veloc)

#ifndef OPENACCHHH
#endif

     do isubd = 1,num_subd

#ifndef OPENACCHHH
#ifdef ALYA_OMPSS
        num_neigh = size(ompss_domains(isubd) % neighbours,KIND=ip)

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
        !$OMP SHARED       ( list_elements,num_pack,lnnod,lgaus,kfl_assem_nsi )      &
        !-----------------------------------------------------------------------------
        !$OMP SHARED       ( itask )                                                 &
        !$OMP REDUCTION    ( +:resis_nsi,itsta_nsi,resgs_nsi,time_detail )           &
        !$OMP REDUCTION    ( MIN:tamin_nsi )                                         &
        !$OMP REDUCTION    ( MAX:rmsgs_nsi,tamax_nsi,dtmax_nsi                       )
        !-----------------------------------------------------------------------------
#endif
        
        do ipack = 1,num_pack(isubd)

           ielem = list_elements(isubd) % packs(ipack) % l(1)  ! Select first element
           pnode = lnnod(ielem)                                ! Number of nodes
           pgaus = lgaus(ielem)                                ! Number of Gauss points
 

           if ( kfl_assem_nsi == 5 ) then
               if(kfl_streams_per_gpu /= 1) then
              
#ifdef OPENACCHHH

                  ngpus = acc_get_num_devices(acc_device_nvidia)
                  streamid=(mod(ipack,kfl_streams_per_gpu) + 1)
                  call nsi_element_operations_fast8(&
                   size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
                   pnode,pgaus,list_elements(isubd) % packs(ipack) % l,time_detail,streamid)
#else
                    call nsi_element_operations_fast(&
                   size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
                   pnode,pgaus,list_elements(isubd) % packs(ipack) %l,time_detail)
#endif
              else
                   call nsi_element_operations_fast(&
                   size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
                   pnode,pgaus,list_elements(isubd) % packs(ipack) % l,time_detail)
              end if
             !#ifndef OPENACCHHH
!#ifdef NINJA
!              call nsi_element_operations_cuda(&
!                   size(list_elements(isubd) % packs(ipack) % l,KIND=ip),&
!                   pnode,pgaus,list_elements(isubd) % packs(ipack) % l,time_detail)
!#endif
!#else
!              call nsi_element_operations_fast5(&
!                   size(list_elements(isubd) % packs(ipack) % l,kind=ip),&
!                   pnode,pgaus,list_elements(isubd) % packs(ipack) % l,time_detail)
!!#endif
!                 end if
           else

              call nsi_element_operations(&
                   itask,pnode,pgaus,list_elements(isubd) % packs(ipack) % l,&
                   resis_nsi,itsta_nsi,resgs_nsi,tamin_nsi,rmsgs_nsi,tamax_nsi,&
                   dtmax_nsi,time_detail)

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
#else

     !$acc wait 
 
#endif
     !
     ! Accumulate CPU time
     !
     cputi_assembly_nsi = cputi_assembly_nsi + time_detail

  end if

end subroutine nsi_elmope_all
