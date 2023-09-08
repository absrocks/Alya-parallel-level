subroutine chm_updunk(itask)
  !-----------------------------------------------------------------------
  !****f* partis/chm_updunk
  ! NAME 
  !    chm_updunk
  ! DESCRIPTION
  !    This routine performs several types of updates 
  ! USED BY
  !    chm_begste (itask=1)
  !    chm_begite (itask=2)
  !    chm_endite (itask=3, inner loop) 
  !    chm_endite (itask=4, outer loop) 
  !    chm_endste (itask=5)
  !    chm_restar (itask=6)
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_ADR,    only : ADR_end_time_step
  use mod_ADR,    only : ADR_begin_inner_iteration
  use mod_ADR,    only : ADR_end_inner_iteration
  use mod_ADR,    only : ADR_begin_time_step
  use mod_ADR,    only : ADR_after_restart

  
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iclas,itime,kpoin,ipoin,icomp,ielem,igaus
  integer(ip)             :: pgaus
  real(rp)                :: total,rela1
  integer(ip)             :: monolithic_dim

  if( INOTMASTER ) then

     select case (itask)

     case(1_ip)
        !
        ! Assign initial guess for outer iterations
        !
        icomp = min(3_ip,ncomp_chm)
        do iclas=1,nspec_chm
           !!DMM-ADRdo ipoin=1,npoin
           !!DMM-ADR   conce(ipoin,iclas,2) = conce(ipoin,iclas,icomp)
           !!DMM-ADRend do
           call ADR_begin_time_step(ADR_chm(iclas),conce(:,iclas,:))
        end do

     case(2_ip)
        !
        ! Assign initial guess for inner iterations
        !
        do iclas = 1,nclas_chm
           !!DMM-ADR do ipoin = 1,npoin
           !!DMM-ADR    conce(ipoin,iclas,1) = conce(ipoin,iclas,2)
           !!DMM-ADR end do
           call ADR_begin_inner_iteration(ADR_chm(iclas),conce(:,iclas,:))
        end do

        if( kfl_coupl_chm /= 0 ) then
           kpoin = 0
           do ipoin = 1,npoin
              do iclas = 1,nclas_chm
                 kpoin = kpoin + 1
                 !!DMM-ADR conce(ipoin,iclas,1) = conce(ipoin,iclas,2)
                 unkno(kpoin)         = conce(ipoin,iclas,2)
              end do
           end do
           do iclas=nclas_chm+1,nspec_chm
              !!DMM-ADR do ipoin=1,npoin
              !!DMM-ADR    conce(ipoin,iclas,1) = conce(ipoin,iclas,2)
              !!DMM-ADR end do
              call ADR_begin_inner_iteration(ADR_chm(iclas),conce(:,iclas,:))
           end do
        end if

     case(3_ip)
        !
        ! Assign update (PERHAPS OPTIMIZE FOR MONOLITHIC)
        !
        kpoin = 0
        do ipoin = 1,npoin
           do iclas = iclai_chm,iclaf_chm
              kpoin = kpoin+1
              conce(ipoin,iclas,1) = unkno(kpoin)
           end do
        end do

     case(4_ip)
        !
        ! Assign 
        !        
        do iclas = 1,nspec_chm
           !!DMM-ADR do ipoin = 1,npoin
           !!DMM-ADR    conce(ipoin,iclas,2) = conce(ipoin,iclas,1)
           !!DMM-ADR end do
           call ADR_end_inner_iteration(ADR_chm(iclas),conce(:,iclas,:))
        end do

     case(5_ip)
        !
        ! Obtain c^n and assign c^{n-1} <- c^n
        !        
        !!DMM-ADR if( kfl_tisch_chm == 1 .and. kfl_tiacc_chm == 2  ) then
        !!DMM-ADR    !
        !!DMM-ADR    ! Crank-Nicolson method 
        !!DMM-ADR    !        
        !!DMM-ADR    do iclas = 1,nspec_chm
        !!DMM-ADR       do ipoin = 1,npoin
        !!DMM-ADR          conce(ipoin,iclas,1) = 2.0_rp*conce(ipoin,iclas,1)-conce(ipoin,iclas,3)
        !!DMM-ADR       end do
        !!DMM-ADR    end do

        !!DMM-ADR else if( kfl_tisch_chm == 2 ) then
        !!DMM-ADR    !
        !!DMM-ADR    ! BDF scheme
        !!DMM-ADR    !
        !!DMM-ADR    do itime = 2+kfl_tiacc_chm,4,-1
        !!DMM-ADR       do iclas = 1,nspec_chm
        !!DMM-ADR          do ipoin = 1,npoin
        !!DMM-ADR             conce(ipoin,iclas,itime) = conce(ipoin,iclas,itime-1)
        !!DMM-ADR          end do
        !!DMM-ADR       end do
        !!DMM-ADR    end do
        !!DMM-ADR end if

        !!DMM-ADR if( kfl_dttyp_chm == 2) then
        !!DMM-ADR    !
        !!DMM-ADR    ! Adaptive time step: save previous concentration
        !!DMM-ADR    !
        !!DMM-ADR    if( kfl_tisch_chm == 1 ) then
        !!DMM-ADR       do iclas = 1,nspec_chm
        !!DMM-ADR          do ipoin = 1,npoin
        !!DMM-ADR             conce(ipoin,iclas,5) = conce(ipoin,iclas,4)
        !!DMM-ADR             conce(ipoin,iclas,4) = conce(ipoin,iclas,3)
        !!DMM-ADR          end do
        !!DMM-ADR       end do
        !!DMM-ADR    else
        !!DMM-ADR       call runend('CHM_UPDUNK: ADAPTIVE DT AND BDF NOT CODED')
        !!DMM-ADR    end if
        !!DMM-ADR end if

        do iclas = 1,nspec_chm
           !!DMM-ADR do ipoin = 1,npoin
           !!DMM-ADR    conce(ipoin,iclas,3) = conce(ipoin,iclas,1)
           !!DMM-ADR end do
           call ADR_end_time_step(ADR_chm(iclas),conce(:,iclas,:))
        end do
 
       !!DMM-ADR if( kfl_sgsti_chm == 1 ) then
       !!DMM-ADR     !
       !!DMM-ADR     ! Time tracking of the subscales
       !!DMM-ADR     !        
       !!DMM-ADR     if( kfl_tiacc_chm == 2 ) then 
       !!DMM-ADR        do iclas = 1,nclas_chm
       !!DMM-ADR           do ielem = 1,nelem
       !!DMM-ADR              pgaus = ngaus(ltype(ielem))
       !!DMM-ADR              do igaus = 1,pgaus
       !!DMM-ADR                 cosgs(ielem,iclas)%a(igaus,1) = &
       !!DMM-ADR                      2.0_rp * cosgs(ielem,iclas)%a(igaus,1) - cosgs(ielem,iclas)%a(igaus,2)               
       !!DMM-ADR              end do
       !!DMM-ADR           end do
       !!DMM-ADR        end do
       !!DMM-ADR     end if

       !!DMM-ADR     do iclas = 1,nclas_chm
       !!DMM-ADR        do ielem = 1,nelem
       !!DMM-ADR           pgaus = ngaus(ltype(ielem))
       !!DMM-ADR           do igaus = 1,pgaus
       !!DMM-ADR              cosgs(ielem,iclas)%a(igaus,2) = cosgs(ielem,iclas)%a(igaus,1)
       !!DMM-ADR           end do
       !!DMM-ADR        end do
       !!DMM-ADR     end do
       !!DMM-ADR  end if

     case(6_ip) 
        !
        ! Assign c(n,:,*) <-- c(n-1,:,*), when reading from restart file
        ! 
        icomp = min(3_ip,ncomp_chm) 

        do iclas = 1,nspec_chm
           !!DMM-ADR do ipoin = 1,npoin
           !!DMM-ADR    conce(ipoin,iclas,1) = conce(ipoin,iclas,icomp)
           !!DMM-ADR end do
           call ADR_after_restart(ADR_chm(iclas),conce(:,iclas,:))  !!DMM Restart should be checked
        end do

     case(7_ip)
        !
        ! Assign update
        !
        if( kfl_coupl_chm == 2 ) then
           kpoin = 0
           do ipoin = 1,npoin
              do iclas = iclai_chm,iclaf_chm             
                 kpoin = kpoin+1
                 unkno(kpoin) = conce(ipoin,iclas,1)
              end do
           end do
        else
           do iclas = iclai_chm,iclaf_chm
              do ipoin = 1,npoin
                 unkno(ipoin) = conce(ipoin,iclas,1)
              end do
           end do
        end if


     case(8_ip)
        !
        ! Prevent undershoots
        !
        if( kfl_coupl_chm == 2 ) then
           monolithic_dim = nspec_chm
        else
           monolithic_dim = 1
        end if
        if (kfl_negat_chm==1) then
           do ipoin = 1,npoin*monolithic_dim
              if( unkno(ipoin) < 0.0_rp ) then
                 kfl_under_chm = kfl_under_chm + 1
                 !unkno(ipoin) = conce(ipoin,iclas_chm,2)
                 unkno(ipoin) = 0.0_rp
              end if
           end do
        endif
       if (kfl_posit_chm==1) then
           do ipoin = 1,npoin*monolithic_dim
              if( unkno(ipoin) > 1.0_rp ) then
                 !kfl_under_chm = kfl_under_chm + 1
                 !unkno(ipoin) = conce(ipoin,iclas_chm,2)
                 unkno(ipoin) = 1.0_rp
              end if
           end do
        endif

     case(9_ip)
        !
        ! Normalize one species in case of errors
        !
        do ipoin =1,npoin 
           total = 0.0_rp
           do iclas = 1,nspec_chm
              if (conce(ipoin,iclas,1) < 0.0_rp) then 
                 kfl_under_chm = kfl_under_chm+1
                 conce(ipoin,iclas,1) = 0.0_rp
              endif
              total = total + conce(ipoin,iclas,1)
           enddo
           if (total .ne. 1.0_rp) then
              kfl_overs_chm = kfl_overs_chm + 1
              total = total - conce(ipoin,kfl_norma_chm,1)  ! We substract the corrected species
              if (total .ge. 0.0_rp .and. total .le. 1.0) then
                 conce(ipoin,kfl_norma_chm,1) = 1.0_rp - total
              else
                 print *,' CHEMIC WARNING: Unfixable overshoot'
              endif
           endif
        end do

     case(10_ip)
        !
        ! Relax update
        !
        if (relax_chm < 1.0_rp) then 
           rela1 = 1.0_rp - relax_chm
           if( kfl_coupl_chm == 2 ) then
              kpoin = 0
              do ipoin = 1,npoin
                 do iclas = iclai_chm,iclaf_chm
                    kpoin = kpoin +1
!!                    unkno(ipoin+(iclas-1)*npoin) = relax_chm * unkno(ipoin+(iclas-1)*npoin) + rela1 * conce(ipoin,iclas,1)
!!                    unkno((ipoin-1)*iclaf_chm+iclas) = relax_chm * unkno((ipoin-1)*iclaf_chm+iclas) + rela1 * conce(ipoin,iclas,1)
                    unkno(kpoin) = relax_chm * unkno(kpoin) + rela1 * conce(ipoin,iclas,1)
                 enddo
              enddo
           else
              do iclas = iclai_chm,iclaf_chm
                 do ipoin = 1,npoin
                    unkno(ipoin) = relax_chm * unkno(ipoin) + rela1 * conce(ipoin,iclas,1)
                 enddo
              enddo
           end if
           
        end if

     case(11_ip)

        do iclas = 1,nspec_chm
           do ipoin = 1,npoin
              conce(ipoin,iclas,3) = conce(ipoin,iclas,1)
           end do
        end do

     end select

  end if

end subroutine chm_updunk

