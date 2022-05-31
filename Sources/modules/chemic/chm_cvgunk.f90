subroutine chm_cvgunk(itask)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_cvgunk
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
  use mod_array_operations, only : array_operations_residual_norm

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass=0,jpass=0
  integer(ip)             :: iclas,jclas,kpoin,ipoin
  real(rp),    save       :: cpuit_chm=0.0_rp
  real(rp)                :: time1,comin,comax
  real(rp)                :: rtpts,dummr

  select case(itask)

  case(ITASK_ENDINN)

     !
     ! Check convergence of the inner iterations:
     ! || c(n,i,j) - c(n,i,j-1)|| / ||c(n,i,j)||
     !     
     do iclas = iclai_chm,iclaf_chm

        if( kfl_normc_chm /= 3 ) then
           ripts_chm(iclas) = array_operations_residual_norm(kfl_normc_chm,nclas_chm,1_ip,unkno,conce,iclas-1,npoin*(iclas-1),1_ip) 
           !!!!!!!!!!!!!!!!!!!!!! CHECK IF THIS IS CORRECT FOR CMC MODEL
           !!!if (kfl_model_chm == 4 .and. kfl_solve_enth_CMC_chm /= 0) then
           !!!   if (iclas > nclas_chm) then
           !!!      call residu(&
           !!!           kfl_normc_chm,nvar_CMC_chm,one,unkno(1),therm(1,1),&
           !!!           iclas,one,one,1.0_rp,ripts_chm(iclas))
           !!!   else
           !!!      call residu(&
           !!!           kfl_normc_chm,nvar_CMC_chm,one,unkno(1),conce(1,iclas,1),&
           !!!           iclas,one,one,1.0_rp,ripts_chm(iclas))
           !!!   end if

           !!!else
           !!!   call residu(&
           !!!        kfl_normc_chm,nclas_chm,one,unkno(1),conce(1,iclas,1),&
           !!!        iclas,one,one,1.0_rp,ripts_chm(iclas))
           !!!end if
        end if

        if( ripts_chm(iclas) < cotol_chm ) then
           kfl_goit2_chm = kfl_goit2_chm + 1
        endif

        rtpts_chm = rtpts_chm + ripts_chm(iclas)

     end do

     !
     ! Stop or go on
     !
     if( kfl_goit2_chm == nclas_chm ) kfl_goite_chm = 0

     !!DMM if( kfl_spray_chm == 2 ) miinn_chm = 0 

     !! DMM if ( itinn(modul) >= miinn_chm ) kfl_goite_chm = 0
     kfl_goite_chm = 0
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

        if ( kfl_model_chm==1 ) then ! Combustion saves species convergence by default   
           !================!
           ! FLAMELET MODEL !
           !================!
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

        elseif ( kfl_model_chm==3 ) then
           !
           !=============================!
           ! FINITE RATE CHEMISTRY MODEL !
           !=============================!
           ! 
           call cputim(time1)
           if( ipass == 0 ) then
              time1 = time1-cpu_initi
           else
              time1 = time1-cpuit_chm
           end if
           if( ipass == 0 .and. kfl_rstar /= 2 ) write(momod(modul)%lun_conve,112)
           write(momod(modul)%lun_conve,113) &
                ittim,itcou,itinn(modul),cutim,rtpts_chm,comin_chm,&
                comax_chm,time1,hrr_int_chm,(ripts_chm(iclas),iclas = 1,nclas_chm)
           call cputim(cpuit_chm)
           flush(momod(modul)%lun_conve)
           ipass = ipass + 1

        endif

     end if
     
  case(ITASK_ENDITE)
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
        ripts_chm(iclas) = array_operations_residual_norm(kfl_normc_chm,1_ip,1_ip,conce,conce,0_ip,1_ip*npoin*nclas_chm,1_ip)       
        !call residu(kfl_normc_chm,one,one,conce(1,jclas,1),conce(1,jclas,2),one,one,one,1.0_rp,ripts_chm(iclas))
        resid_chm = resid_chm + ripts_chm(iclas)
     end do

  case(ITASK_ENDSTE)
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
        ripts_chm(iclas) = array_operations_residual_norm(kfl_normc_chm,1_ip,1_ip,conce,conce,0_ip,2_ip*npoin*nclas_chm,1_ip)       
        !call residu(kfl_normc_chm,one,one,conce(1,jclas,1),conce(1,jclas,3),one,one,one,1.0_rp,ripts_chm(iclas))
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

112 format('# --| ALYA Convergence, Combustion model '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   ' ,/,&
       &   '# --| 4. Current time      5. Concentration      6. Min. value        ' ,/,&
       &   '# --| 7. Max. value        8. Elapsed CPU Time   9. Heat release (W)  ' ,/,&
       &   '# --|10. Residual Spec 1,2,... '                                        ,/,&
       &   '# ','          1','          2','          3',&
       &        '             4','             5','             6','             7',&
       &        '             8','             9','            10...')
113 format(4x,i9,2x,i9,2x,i9,100(2x,e12.6))

end subroutine chm_cvgunk

