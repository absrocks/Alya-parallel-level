subroutine nsi_restar(itask)
  !------------------------------------------------------------------------
  !****f* Nastin/nsi_restar
  ! NAME 
  !    nsi_restar
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... Reads the initial values from the restart file
  !            2 ... Writes restart file
  ! USES
  ! USED BY
  !    nsi_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_nastin
  use mod_postpr
  use mod_memchk
  use def_kermod, only : kfl_wlaav_ker, velav_ker, kfl_waexl_ker, kfl_twola_ker, kfl_noslw_ker, avta1_nsw_ker,&
       kount_nsw_ele_ker,fact_nsw_ker,avupo_ker,kfl_fixbo_nsw_ker
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: icomp,kfl_gores,idime,igaus
  integer(ip)             :: ielem,pgaus,pelty,pblty
!  character(5)            :: wopos(3)
  integer(ip)             :: ipoin,icount,iboun,igaub,pgaub
  real(rp), pointer       :: avvre(:,:)
  real(rp), pointer       :: grnor(:), ubpre(:)
  real(rp), pointer       :: bpess_aux(:,:)
  logical(lg)             :: lauxi
  !
  ! Check if restrt file should be read or written
  !
  call respre(itask,kfl_gores)
  if( kfl_gores == 0 ) return

  if( itask == READ_RESTART_FILE ) then
     icomp = nprev_nsi
  else
     icomp = 1
  end if

  !----------------------------------------------------------------------
  !
  ! Velocity and pressure 
  !
  !----------------------------------------------------------------------
  gevec => veloc(:,:,icomp)
  call postpr(gevec,postp(1)%wopos(1:3,1),ittim,cutim)
  gesca => press(:,icomp)
  call postpr(gesca,postp(1)%wopos(1:3,2),ittim,cutim)

  if (kfl_tisch_nsi==2) then    ! BDF
     if (ncomp_nsi > 3) then    ! BDF2 and bigger
        gevec => veloc(:,:,4)
        call postpr(gevec,postp(1)%wopos(1:3,41),ittim,cutim)

        gesca => press(:,4)
        call postpr(gesca,postp(1)%wopos(1:3,42),ittim,cutim)
     end if

     if (ncomp_nsi > 4) then  ! BDF3 and bigger ! beware for BDF4 and bigger more components should be written
        gevec => veloc(:,:,5)
        call postpr(gevec,postp(1)%wopos(1:3,43),ittim,cutim)

        gesca => press(:,5)
        call postpr(gesca,postp(1)%wopos(1:3,44),ittim,cutim)
     end if
  end if

  if (kfl_tisch_nsi==4) then       ! Multistep FS
     if(kfl_fscon_nsi == 0) then   ! NON-CONSISTENT - Capuano approximation this is where we need the n-1 pressure
        gesca => press(:,4)
        call postpr(gesca,postp(1)%wopos(1:3,42),ittim,cutim)
     end if

  end if
  

  !----------------------------------------------------------------------
  !
  ! Hydrostatic values
  !
  !----------------------------------------------------------------------

  if( kfl_hydro_nsi /= 0 ) then
     bpess_aux => bpess_nsi(:,:,1)
     call postpr(bpess_aux,postp(1) % wopos(1:3,31),ittim,cutim)
  end if
  if( kfl_hydro_gravity_nsi /= 0 ) then
     call runend('NSI_RESTAR: NOT CODED')
     ! should postprocess NSI_HYDRO_DENSITY or not
     !wopos(1) = 'FLEVE'
     !wopos(2) = 'SCALA'
     !wopos(3) = 'NPOIN'
     !call postpr(,wopos,ittim,cutim)    
  end if

  !----------------------------------------------------------------------
  !
  ! Time-averaged velocity at boundary gauss points for computing wall-law
  !
  !----------------------------------------------------------------------

  if( kfl_wlaav_ker /= 0 ) then 
      allocate(avvre(ndime,nboun*mgaub))
      avvre = 0.0_rp
      if (itask == WRITE_RESTART_FILE .and. INOTMASTER ) then
         do iboun=1,nboun
              lauxi = .false.        ! to allocate kfl_fixbo_nsw_ker only if kfl_noslw_ker /= 0
              if( kfl_noslw_ker /= 0 ) then
                 if( kfl_fixbo_nsw_ker(iboun) /= 0_ip ) lauxi = .true. ! no slip wall law uses the variational force not the one from gradients
              end if

            
            if(  kfl_fixbo_nsi(iboun) ==  3 .or. &    ! Wall law
               & kfl_fixbo_nsi(iboun) == 13 .or. &    ! Wall law + open pressure
               & kfl_fixbo_nsi(iboun) == 18 .or. &    ! u.n in weak form
               & lauxi)  then                         ! no slip walllaw

               pblty = ltypb(iboun)
               pgaub = ngaus(pblty)

               do igaub=1,pgaub
                  icount = igaub + pgaub * (iboun-1)
                  avvre(1:ndime,icount)=velav_ker(1:ndime,igaub,iboun)
               end do
            end if
         end do
      end if
      call postpr(avvre,postp(1)%wopos(1:3,72),ittim,cutim, pleng_opt=nboun*mgaub )
      if ( kfl_waexl_ker == 0_ip .or. kfl_noslw_ker == 1_ip ) then
         call postpr(avupo_ker,postp(1)%wopos(1:3,73),ittim,cutim)
      end if
      if (itask == READ_RESTART_FILE .and. INOTMASTER ) then
         do iboun=1,nboun
            lauxi = .false.        ! to allocate kfl_fixbo_nsw_ker only if kfl_noslw_ker /= 0
            if( kfl_noslw_ker /= 0 ) then
               if( kfl_fixbo_nsw_ker(iboun) /= 0_ip ) lauxi = .true. ! no slip wall law uses the variational force not the one from gradients
            end if
            if(  kfl_fixbo_nsi(iboun) ==  3 .or. &    ! Wall law
               & kfl_fixbo_nsi(iboun) == 13 .or. &    ! Wall law + open pressure
               & kfl_fixbo_nsi(iboun) == 18 .or. &    ! u.n in weak form
               & lauxi)  then                         ! no slip walllaw

               pblty = ltypb(iboun)
               pgaub = ngaus(pblty)

               do igaub=1,pgaub
                  icount = igaub + pgaub * (iboun-1)
                  do idime=1,ndime
                     if (abs(avvre(idime,icount)) > 1e-15_rp ) then  ! Needs improving
                        velav_ker(idime,igaub,iboun)=avvre(idime,icount)
                        kfl_wlare_nsi = 1_ip ! Restart file exists, no need to initialize
                     end if
                  end do
               end do
            end if
         end do
      end if
      deallocate(avvre)
   end if


  !----------------------------------------------------------------------
  !
  ! Time-averaged values used for no slip wall law viscous part of tangential stress
  !
  !----------------------------------------------------------------------

  if( kfl_wlaav_ker /= 0 .and. kfl_noslw_ker /= 0_ip ) then 
     call postpr(avta1_nsw_ker,postp(1)%wopos(1:3,100),ittim,cutim, pleng_opt=kount_nsw_ele_ker )  ! running average viscous(lam+tur) part of tangential stress
     call postpr(avvaf_nsi,postp(1)%wopos(1:3,102),ittim,cutim)     ! running average  VARIATIONAL FORCE 
     call postpr(avntr_nsi,postp(1)%wopos(1:3,104),ittim,cutim)     ! running average  VARIATIONAL Tangential traction 
     call postpr(avgtr_nsi,postp(1)%wopos(1:3,105),ittim,cutim)     ! running average  gradient based Tangential traction 
     call postpr(fact_nsw_ker,postp(1)%wopos(1:3,106),ittim,cutim)  ! factor for viscosity aclculation 
  end if
 
  !----------------------------------------------------------------------
  !
  ! Projections for OSS stabilization
  !
  !----------------------------------------------------------------------

  if( kfl_stabi_nsi > 0 ) then ! oss method
     !
     ! Velocity
     !
     call postpr(vepro_nsi,postp(1)%wopos(1:3,23),ittim,cutim)
     !
     ! Velocity divergence
     !
     call postpr(prpro_nsi,postp(1)%wopos(1:3,24),ittim,cutim)

  end if

  if( kfl_stabi_nsi == 2 ) then
     !
     ! Pressure gradients
     !
     call postpr(grpro_nsi,postp(1)%wopos(1:3,27),ittim,cutim)

  end if

  !----------------------------------------------------------------------
  !
  ! SGS for tracking
  !
  !----------------------------------------------------------------------

  if( kfl_sgsco_nsi /= 0 .or. kfl_sgsti_nsi /= 0 ) then
     !
     ! Velocity SGS
     ! 
     call postpr(vesgs,postp(1)%wopos(1:3,5),ittim,cutim,1_ip)

     if(itask == READ_RESTART_FILE .and. INOTMASTER ) then
        do ielem = 1,nelem
           pelty = abs(ltype(ielem))
           pgaus = ngaus(pelty)
           do igaus = 1,pgaus
              do idime = 1,ndime
                 vesgs(ielem)%a(idime,igaus,2) = vesgs(ielem)%a(idime,igaus,1)
              end do
           end do
        end do
     end if

  end if
 
  !----------------------------------------------------------------------
  !
  ! Mass flow control restart
  !
  !---------------------------------------------------------------------- 

  if ( kfl_mfrco_nsi > 0 ) then ! Workaround way - Could be improved
     allocate(grnor(npoin))
     allocate(ubpre(npoin))
     if(itask == WRITE_RESTART_FILE .and. INOTMASTER ) then
        do ipoin = 1,npoin
           grnor(ipoin) = grnor_nsi
           ubpre(ipoin) = ubpre_nsi
        end do
     end if
     call postpr(grnor,postp(1)%wopos(1:3,81),ittim,cutim)
     call postpr(ubpre,postp(1)%wopos(1:3,82),ittim,cutim)
     if (itask == READ_RESTART_FILE .and. INOTMASTER ) then
        grnor_nsi = grnor(1)
        ubpre_nsi = ubpre(1)
     end if
     deallocate(grnor)
     deallocate(ubpre)
  end if
 
  !----------------------------------------------------------------------
  !
  ! Restart for the RANS/LES two layer model
  !
  !----------------------------------------------------------------------
 
  if ( kfl_twola_ker > 0 ) then
     if ( kfl_twola_ker == 1_ip ) then
        call postpr(btrac_nsi,postp(1)%wopos(1:3,96),ittim,cutim)
        call postpr(tluav_nsi,postp(1)%wopos(1:3,98),ittim,cutim)
     else if ( kfl_twola_ker == 2_ip ) then
        call postpr(tracr_nsi,postp(1)%wopos(1:3,97),ittim,cutim)
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! assign constant forward velocity and pressure values for adjoint
  !
  !----------------------------------------------------------------------  
!   if( itask == READ_RESTART_FILE .and. kfl_adj_prob == 1 ) then
!     do ipoin = 1, npoin
!       press_forw(ipoin,1) = press(ipoin,nprev_nsi)
!       press(ipoin,nprev_nsi) = 0.0_rp
!       do idime = 1,ndime
! 	  veloc_forw(idime,ipoin,1) = veloc(idime,ipoin,nprev_nsi)
! 	  veloc(idime,ipoin,nprev_nsi) = 0.0_rp
!       end do
!     end do
!   endif
  
  !
  ! Finish
  !
  call respre(3_ip,kfl_gores)
  
end subroutine nsi_restar
 
