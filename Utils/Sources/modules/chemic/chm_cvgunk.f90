subroutine chm_cvgunk(itask)
  !-----------------------------------------------------------------------
  !****f* partis/chm_cvgunk
  ! NAME 
  !    chm_cvgunk
  ! DESCRIPTION
  !    Convergence criterion
  ! USES
  ! USED BY
  !    chm_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_postpr, only : postpr ! temporary
  use mod_outfor, only : outfor

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass=0,jpass=0
  integer(ip)             :: iclas,jclas,kpoin,ipoin
  real(rp),    save       :: cpuit_chm=0.0_rp
  real(rp)                :: time1,comin,comax
  real(rp)                :: rtpts,dummr

  select case(itask)

  case(1_ip)
     !
     ! Check convergence of the inner iterations:
     ! || c(n,i,j) - c(n,i,j-1)|| / ||c(n,i,j)||
     !     
     if (kfl_coupl_chm /= 2) then !! Not monolithic
        kpoin = 1
        do iclas = iclai_chm,iclaf_chm
           if( kfl_normc_chm /= 3 ) then
              call residu(&
                   kfl_normc_chm,one,one,unkno(kpoin),conce(1,iclas,1),&
                   one,one,one,1.0_rp,ripts_chm(iclas))     
              if( INOTMASTER ) kpoin = kpoin + npoin
           end if
           if( ripts_chm(iclas) < cotol_chm ) then
              kfl_gocla_chm = kfl_spite_chm ! End intra-species iterations
              kfl_goit2_chm = kfl_goit2_chm + 1
           endif
           rtpts_chm = rtpts_chm + ripts_chm(iclas)
        end do
     else !!Monolithic
        do iclas = iclai_chm,iclaf_chm
           if( kfl_normc_chm /= 3 ) then
              call residu(&
                   kfl_normc_chm,nspec_chm,one,unkno(1),conce(1,iclas,1),&
                   iclas,one,one,1.0_rp,ripts_chm(iclas))     
           end if
           if( ripts_chm(iclas) < cotol_chm ) then
              kfl_gocla_chm = kfl_spite_chm ! End intra-species iterations
              kfl_goit2_chm = kfl_goit2_chm + 1
           endif
           rtpts_chm = rtpts_chm + ripts_chm(iclas)
        end do
     endif
        !write(90+iclas,'(10(1x,e12.6))') ripts_chm(iclas)
        !if( iclas == 1 ) then
        !   ii=ii+1
        !   wopos(1) = 'RESID'
        !   wopos(2) = 'SCALA'
        !   wopos(3) = 'NPOIN'
        !   allocate(gesca(npoin))
        !   do ipoin = 1,npoin
        !      gesca(ipoin) = abs(unkno(ipoin)-conce(ipoin,iclas,1))
        !   end do
        !   call postpr(gesca,wopos,ii,real(ii,rp))   
        !   deallocate(gesca)
        !end if
     !
     ! Stop or go on
     !
     if( kfl_goit2_chm == nclas_chm ) kfl_goite_chm = 0
     if ( itinn(modul) >= miinn_chm ) kfl_goite_chm = 0
     !
     ! Compute min and max 
     !
     kpoin    =  1
     do iclas = iclai_chm,iclaf_chm
        if( IMASTER ) then
           call minmax(one,npoin,zero,dummr,comin,comax)
        else
           call minmax(one,npoin,zero,unkno(kpoin),comin,comax)
        end if
        comin_chm = min(comin_chm,comin)
        comax_chm = max(comax_chm,comax)
        if( INOTMASTER ) kpoin = kpoin + npoin
     end do
     !
     ! Write convergence
     !
     if( INOTSLAVE ) then
        if (kfl_model_chm==4 .or. kfl_model_chm==5 ) then ! Combustion saves species convergence by default   
           call cputim(time1)
           if( ipass == 0 ) then
              time1 = time1-cpu_initi
           else
              time1 = time1-cpuit_chm
           end if
           if( ipass == 0 .and. kfl_rstar /= 2 ) write(momod(modul)%lun_conve,110)
           write(momod(modul)%lun_conve,111) &
                ittim,itcou,itinn(modul),cutim,rtpts_chm,comin_chm,&
                comax_chm,time1,(ripts_chm(iclas),iclas = 1,nclas_chm)
           call cputim(cpuit_chm)
           flush(momod(modul)%lun_conve) 
           ipass = ipass + 1
        else 
           if( ( kfl_coupl_chm == 0 .and. iclas_chm == nclas_chm ) .or. kfl_coupl_chm == 1 ) then
              call cputim(time1)
              if( ipass == 0 ) then
                 time1 = time1-cpu_initi
              else
                 time1 = time1-cpuit_chm
              end if
              if( ipass == 0 .and. kfl_rstar /= 2 ) write(momod(modul)%lun_conve,100)
              write(momod(modul)%lun_conve,101) &
                   ittim,itcou,itinn(modul),cutim,rtpts_chm,comin_chm,&
                   comax_chm,time1
              call cputim(cpuit_chm)
              flush(momod(modul)%lun_conve) 
              ipass = ipass + 1
           end if
        endif
     end if
     !
     ! Write species convergence
     !
     if( kfl_model_chm == 3 ) then
        if( INOTSLAVE ) then
           if( ( kfl_coupl_chm == 0 .and. iclas_chm == 1 ) .or. kfl_coupl_chm == 1 ) then
              if( jpass == 0 .and. kfl_rstar /= 2 ) write(lun_spcvg_chm,102)
              write(lun_spcvg_chm,103,advance='no') &
                   ittim,itcou,itinn(modul),cutim
              jpass = 1
           end if
        end if
        write(lun_spcvg_chm,104,advance='no') ripts_chm(iclai_chm)
     end if

  case(2_ip)
     !
     ! Check convergence of the outer iterations:
     ! || c(n,i,*) - c(n,i-1,*)|| / ||c(n,i,*)||
     !
     resid_chm = 0.0_rp
     do iclas = 1,nclas_chm
        if( IMASTER ) then 
           jclas = 1
        else
           jclas = iclas
        end if
        call residu(kfl_normc_chm,one,one,conce(1,jclas,1),conce(1,jclas,2),one,one,one,1.0_rp,ripts_chm(iclas))
        resid_chm = resid_chm + ripts_chm(iclas)
     end do

  case(3_ip)
     !
     ! Check residual of the time iterations:
     ! || c(n,*,*) - c(n-1,*,*)|| / ||c(n,*,*)||
     !
     rtpts = 0.0_rp
     do iclas = 1,nclas_chm
        if( IMASTER ) then 
           jclas = 1
        else
           jclas = iclas
        end if
        call residu(kfl_normc_chm,one,one,conce(1,jclas,1),conce(1,jclas,3),one,one,one,1.0_rp,ripts_chm(iclas))
        rtpts = rtpts + ripts_chm(iclas)
     end do
     rtpts = rtpts / real(nclas_chm,rp)
     if( rtpts <= sstol_chm ) then
        momod(modul) % kfl_stead = 1
        call outfor(28_ip,momod(modul)%lun_outpu,' ')
     end if

  end select
  !
  ! Formats
  !
100 format('# --| ALYA Convergence '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   ' ,/,&
       &   '# --| 4. Current time      5. Concentration      6. Min. value        ' ,/,& 
       &   '# --| 7. Max. value        8. Elapsed CPU Time  ' ,/,&
       &   '# ','          1','          2','          3',&
       &        '             4','             5','             6','             7',&
       &        '             8','             9')
102 format('# --| ALYA Convergence '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   '  ,/,&
       &   '# --| 4. Current time      5. Concentration 1...                      '  ,/,&
       &   '# ','          1','          2','          3',&
       &        '             4','             5...')
101 format(4x,i9,2x,i9,2x,i9,100(2x,e12.6))
103 format(4x,i9,2x,i9,2x,i9,2x,e12.6)
104 format(2x,e12.6)

110 format('# --| ALYA Convergence, Combustion model '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   ' ,/,&
       &   '# --| 4. Current time      5. Concentration      6. Min. value        ' ,/,& 
       &   '# --| 7. Max. value        8. Elapsed CPU Time   9. Residual Spec 1,2,... ' ,/,&
       &   '# ','          1','          2','          3',&
       &        '             4','             5','             6','             7',&
       &        '             8','             9...')
111 format(4x,i9,2x,i9,2x,i9,100(2x,e12.6))

end subroutine chm_cvgunk

