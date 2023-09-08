subroutine wav_cvgunk(itask)
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_cvgunk
  ! NAME 
  !    wav_cvgunk
  ! DESCRIPTION
  !    This routine compute the residual
  ! USES
  !    residu
  ! USED BY
  !    wav_endite
  !    wav_endste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_wavequ
  use mod_outfor, only : outfor
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass=0
  real(rp)                :: riwav,time1

  select case(itask)

  case(1)
     !
     ! Check convergence of the inner iterations:
     ! || u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||
     !
     call residu(kfl_normc_wav,one,one,unkno,wavam,one,one,one,1.0_rp,riwav)     
     if((riwav<cotol_wav).or.(itinn(modul)>=miinn_wav)) kfl_goite_wav = 0
     !
     ! Compute min and max of the wave amplitude
     !
     call minmax(one,npoin,zero,unkno,wamin_wav,wamax_wav)

     if(kfl_paral<=0) then
        call cputim(time1)
        if(ipass==0.and.kfl_rstar/=2) then
           ipass=1
           write(momod(modul)%lun_conve,100)
        else
           time1=time1-cpuit_wav           
        end if
        write(momod(modul)%lun_conve,101) ittim,itcou,itinn(modul),cutim,riwav,&
             wamin_wav,wamax_wav,time1
        call cputim(cpuit_wav)
        flush(momod(modul)%lun_conve)
     end if

  case(2)
     !
     ! Check convergence of the outer iterations:
     ! || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||
     !
     call residu(kfl_normc_wav,one,one,wavam(1,1),wavam(1,2),one,one,one,1.0_rp,resid_wav)

  case(3)
     !
     ! Check residual of the time iterations:
     ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||
     !
     call residu(kfl_normc_wav,one,one,wavam(1,1),wavam(1,3),one,one,one,1.0_rp,riwav)
     if(riwav<=sstol_wav) then
        kfl_stead_wav = 1
        call outfor(28_ip,momod(modul)%lun_outpu,' ')
     end if
  end select
  !
  ! Formats
  !
100 format('# --| ALYA Convergence  ' ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   '  ,/,&
       &   '# --| 4. Current time      5. Wave amplitude     6. Min. wave ampl.   ' ,/,& 
       &   '# --| 7. Max. wave ampl.   8. Elapsed CPU time   ',//,&
       &   '$ ','          1','          2','          3',&
       &        '             4','             5','             6','             7',&
       &        '             8')
101 format(4x,i9,2x,i9,2x,i9,11(2x,e12.6))

end subroutine wav_cvgunk


