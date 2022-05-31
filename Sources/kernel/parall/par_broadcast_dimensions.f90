!-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_broadcast_dimensions.f90
!> @author  houzeaux
!> @date    2019-10-18
!> @brief   Broadcast dimensions
!> @details Broadcast meta data to slave. There are all dimensions
!>          required to dimension domain arrays read in reageo, reabcs,
!>          reaset and reafie
!> @}
!-----------------------------------------------------------------------

subroutine par_broadcast_dimensions()

  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use mod_memory,          only : memory_alloca
  use mod_memory,          only : memory_deallo
  use mod_communications,  only : PAR_EXCHANGE
  use mod_communications,  only : PAR_BROADCAST
  use mod_parall,          only : par_memor
  use def_mpio,            only : mpio_flag_geometry_export, PAR_MPIO_ON
  use mod_domain,          only : domain_memory_allocate
  implicit none

  integer(ip) :: ii,jj,ifiel
  integer(ip) :: nsteps

  nullify(parin)
  nullify(parre)

  !----------------------------------------------------------------------
  !
  ! Fixed dimensions
  !
  !----------------------------------------------------------------------

  do parii = 1,2
     npari = 0
     nparr = 0
     !
     ! Read in readim
     !
#ifndef NDIMEPAR 
     call PAR_EXCHANGE(ndime,parin,npari,parii)
#endif
     !
     ! Slaves do not have geometry yet! do not exchange mesh data
     !
     call PAR_EXCHANGE(nhang,parin,npari,parii)
     call PAR_EXCHANGE(nzone,parin,npari,parii)
     call PAR_EXCHANGE(nfiel,parin,npari,parii)
     call PAR_EXCHANGE(mcodb,parin,npari,parii)
     call PAR_EXCHANGE(mnode,parin,npari,parii)
     call PAR_EXCHANGE(kfl_autbo,parin,npari,parii)
     call PAR_EXCHANGE(nexis,parin,npari,parii)
     call PAR_EXCHANGE(utype,parin,npari,parii)
     do jj = 1,size(kfl_field,2,KIND=ip)
        do ii = 1,size(kfl_field,1,KIND=ip)
           call PAR_EXCHANGE(kfl_field(ii,jj),parin,npari,parii)
        end do
     end do
     do ii = 1,nelty
        call PAR_EXCHANGE(lexis(ii),parin,npari,parii)
     end do
     !
     ! Read in reastr
     !
     do ii = 1,nelty
        call PAR_EXCHANGE(lquad(ii),parin,npari,parii)
     end do
     do ii = 1,nelty
        call PAR_EXCHANGE(ngaus(ii),parin,npari,parii)
     end do
     call PAR_EXCHANGE(ngrou_dom,parin,npari,parii)
     call PAR_EXCHANGE(ngrou_dom_target,parin,npari,parii)
     call PAR_EXCHANGE(kfl_ngrou,parin,npari,parii)
     call PAR_EXCHANGE(ngrou_boxes_coarse,parin,npari,parii)
     call PAR_EXCHANGE(ngrou_boxes_fine,parin,npari,parii)

     do ii = 1,3
        call PAR_EXCHANGE(xscal(ii),parre,nparr,parii)
     end do
     do ii = 1,3
        call PAR_EXCHANGE(trans(ii),parre,nparr,parii)
     end do

     call PAR_EXCHANGE(kfl_extra,parin,npari,parii)
     call PAR_EXCHANGE(kfl_geome,parin,npari,parii)
     call PAR_EXCHANGE(kfl_convx,parin,npari,parii)
     call PAR_EXCHANGE(kfl_frees,parin,npari,parii)
     call PAR_EXCHANGE(awind,parre,nparr,parii)
     call PAR_EXCHANGE(tolan,parre,nparr,parii)
     call PAR_EXCHANGE(geoan,parre,nparr,parii) 
     do ii = 1,size(npbcs,KIND=ip)
        call PAR_EXCHANGE(npbcs(ii),parin,npari,parii)
     end do
     do jj = 1,size(lsbcs,2,KIND=ip)
        do ii = 1,size(lsbcs,1,KIND=ip)
           call PAR_EXCHANGE(lsbcs(ii,jj),parin,npari,parii)
        end do
     end do

     call PAR_EXCHANGE(kfl_chege,parin,npari,parii)
     call PAR_EXCHANGE(kfl_naxis,parin,npari,parii)
     call PAR_EXCHANGE(kfl_spher,parin,npari,parii)
     call PAR_EXCHANGE(kfl_bouel,parin,npari,parii)
     call PAR_EXCHANGE(kfl_divid,parin,npari,parii)
     call PAR_EXCHANGE(curvatureDataField,parin,npari,parii)
     call PAR_EXCHANGE(curvatureField,parin,npari,parii)

     do ii = 1,size(materials_icode,KIND=ip)
        call PAR_EXCHANGE(materials_nlaye(ii),parin,npari,parii)
        call PAR_EXCHANGE(materials_icode(ii),parin,npari,parii)
        call PAR_EXCHANGE(materials_imate(ii),parin,npari,parii)
     end do
     !
     ! Allocate and broadcast
     !
     if( parii == 1 ) then
        call memory_alloca(par_memor,'PARIN','ker_parall',parin,npari,'DO_NOT_INITIALIZE')
        call memory_alloca(par_memor,'PARRE','ker_parall',parre,nparr,'DO_NOT_INITIALIZE')
     end if
     if( ( ISLAVE .and. parii == 1 ) .or. ( IMASTER .and. parii == 2 ) ) then
        call PAR_BROADCAST(parin)
        call PAR_BROADCAST(parre)
     end if
  end do

  call memory_deallo(par_memor,'PARRE','ker_parall',parre)
  call memory_deallo(par_memor,'PARIN','ker_parall',parin)

  !----------------------------------------------------------------------
  !
  ! Variable dimensions (Special treatment of fields)
  !
  !----------------------------------------------------------------------

  if( nfiel > 0 ) then

     if( ISLAVE ) then
        do ifiel = 1,nfiel
           if( kfl_field(4,ifiel) > 0 ) then
              call domain_memory_allocate('TIME_FIELD % A',NUMBER1=ifiel)
           end if
        end do
     end if

     do parii = 1,2
        npari = 0
        nparr = 0
        do ifiel = 1,nfiel
           if( kfl_field(4,ifiel) > 0 ) then

              if ( (kfl_field(6,ifiel) /= 1) .OR. (mpio_flag_geometry_export == PAR_MPIO_ON) ) then 
                 nsteps = kfl_field(4,ifiel)
              else
                 nsteps = nsteps_fiel_ondemand
              end if

              do ii = 1,nsteps
                 call PAR_EXCHANGE(time_field(ifiel)%a(ii),parre,nparr,parii)
              end do
           end if
        end do
        if( parii == 1 ) then
           call memory_alloca(par_memor,'PARIN','ker_parall',parin,npari,'DO_NOT_INITIALIZE')
           call memory_alloca(par_memor,'PARRE','ker_parall',parre,nparr,'DO_NOT_INITIALIZE')
        end if
        if( ( ISLAVE .and. parii == 1 ) .or. ( IMASTER .and. parii == 2 ) ) then
           call PAR_BROADCAST(parin)
           call PAR_BROADCAST(parre)
        end if
     end do

     call memory_deallo(par_memor,'PARRE','ker_parall',parre)
     call memory_deallo(par_memor,'PARIN','ker_parall',parin)
  end if

end subroutine par_broadcast_dimensions
