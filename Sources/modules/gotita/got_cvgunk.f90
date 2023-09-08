subroutine got_cvgunk(itask)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_cvgunk
  ! NAME 
  !    got_cvgunk
  ! DESCRIPTION
  !    This routine performs several convergence checks for GOTITA
  ! USES
  !    got_endite (itask=1,2)
  !    got_endste (itask=3)
  ! USED BY
  !    Gotita
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_gotita
  use mod_communications
  use mod_outfor, only : outfor
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass=0
  integer(ip)             :: icomp,ii,mcomp
  real(rp)                :: numer,denom,rialp,time1,vemin,vemax
  real(rp)                :: almin,almax
  real(rp),    save       :: rivel,cpuit_got=0.0_rp
  real(rp),    target     :: taumm(1)
  !
  ! Initializations
  !
  numer = 0.0_rp
  denom = 0.0_rp

  if(itask==1.or.itask==4) then
     !
     ! Check convergence of the inner iterations:
     ! || u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||
     !
     if(itask==1) then
        if(kfl_algor_got==4) then
           call residu(kfl_normc_got,ndime,ndime,unkno,vdrop,one,one,ndime,relax_got,rivel)
           call residu(kfl_normc_got,one,one,unkno(ndime*npoin+1),cdrop,one,one,one,1.0_rp,rialp)
        else
           if(kfl_probl_got==1) then
              call residu(kfl_normc_got,ndofn_got,ndime,unkno,vdrop,one,one,ndime,relax_got,rivel)
              call residu(kfl_normc_got,ndofn_got,one,  unkno,cdrop,ndime+1,one,one,1.0_rp,rialp)
           else if(kfl_probl_got==2) then
              call residu(kfl_normc_got,ndime,ndime,unkno,vdrop,one,one,ndime,relax_got,rivel)
           else if(kfl_probl_got==3) then
              call residu(kfl_normc_got,one,one,unkno,cdrop,one,one,one,relax_got,rialp)
           end if
        end if
        if(rivel<cotol_got.or.itinn(modul)>=miinn_got) kfl_goite_got = 0 
     else
        if(ivari_got==1) then
           call residu(kfl_normc_got,ndime,ndime,unkno,vdrop,one,one,ndime,relax_got,rivel)
           if(rivel<cotol_got.or.itinn(modul)>=miinn_got) kfl_goite_got = 0 
        else
           call residu(kfl_normc_got,one,one,unkno(ndime*npoin+1),cdrop,one,one,one,relax_got,rialp)
           if(rialp<cotol_got.or.itinn(modul)>=miinn_got) kfl_goite_got = 0 
        end if     
     end if
     !
     ! Compute min and max of the velocity and alpha
     !
     if(kfl_algor_got==1) then
        if(kfl_probl_got==1) then
           call minmax(ndofn_got(3),npoin,        ndime,unkno,vemin,vemax)
           call minmax(ndofn_got(3),npoin,-ndofn_got(3),unkno,almin,almax)
        else if(kfl_probl_got==2) then
           call minmax(ndime,npoin,ndime,unkno,vemin,vemax)
        else if(kfl_probl_got==3) then
           call minmax(one,npoin,one,unkno,almin,almax)
        end if
     else
        call minmax(ndime,npoin,     ndime,unkno,vemin,vemax)
        call minmax( 1_ip,npoin,      1_ip,unkno,almin,almax)        
     end if
     !
     ! Minimum and maximum tau
     !
     if(kfl_paral>=0) then
        call PAR_MIN(tamin_got)
        call PAR_MAX(tamax_got)
     end if

     if(kfl_paral<=0) then
        call cputim(time1)
        if(ipass==0.and.kfl_rstar/=2) write(momod(modul)%lun_conve,100)
        if(ipass==1) then
           time1=time1-cpuit_got
        else
           time1=time1-cpu_initi
        end if

        write(momod(modul)%lun_conve,101)&
             ittim,itcou,itinn(modul),cutim,rivel,rialp,&
             vemin,vemax,almin,almax,tamin_got,tamax_got,&
             dimin_got,dimax_got,time1
        call cputim(cpuit_got)
        flush(momod(modul)%lun_conve)
        !
        ! Subgrid scale statistics for convergence
        !
        if(kfl_sgsco_got/=0) then
           if(ipass==0.and.kfl_rstar/=2) then
              write(lun_stasg_got,200) tosgs_got,misgs_got
              write(lun_cvgsg_got,300)               
           end if
           ii=0
           do icomp=1,misgs_got
              ii=ii+itsta_got(icomp)
           end do
           if(ii==0) ii=1
           write(lun_stasg_got,201) &
                (100.0_rp*real(itsta_got(icomp))/real(ii),icomp=1,misgs_got)
           do mcomp=misgs_got,1,-1
              if(itsta_got(mcomp)/=0) exit
           end do
           do icomp=1,mcomp
              if(resis_got(2,icomp)<1.0e-10) resis_got(2,icomp)=1.0_rp
              resis_got(1,icomp)=100.0_rp*sqrt(resis_got(1,icomp)/resis_got(2,icomp))
              write(lun_cvgsg_got,301) &
                   ittim,itcou,itinn(modul),icomp,cutim,resis_got(1,icomp)
           end do
           flush(lun_stasg_got)
           flush(lun_cvgsg_got)
        end if
        ipass=1       
     end if

  else if(itask==2) then
     !
     ! Check convergence of the outer iterations:
     ! || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||
     !
     if(kfl_probl_got==1) then
        call residu(kfl_normc_got,ndime,ndime,vdrop(1,1,1),vdrop(1,1,2),one,one,ndime,1.0_rp,resiv_got)
        call residu(kfl_normc_got,one,  one,  cdrop(1,1),  cdrop(1,2),  one,one,one,  1.0_rp,resic_got)
     else if(kfl_probl_got==2) then
        call residu(kfl_normc_got,ndime,ndime,vdrop(1,1,1),vdrop(1,1,2),one,one,ndime,1.0_rp,resiv_got)
     else if(kfl_probl_got==3) then
        call residu(kfl_normc_got,one,  one,  cdrop(1,1),  cdrop(1,2),  one,one,one,  1.0_rp,resic_got)
     end if

  else if(itask==3) then
     !
     ! Check residual of the time iterations:
     ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||
     !
     if(kfl_probl_got/=3) then
        call residu(kfl_normc_got,ndime,ndime,vdrop,vdrop(1,1,3),one,one,ndime,1.0_rp,rivel)
     else
        call residu(kfl_normc_got,one,one,cdrop,cdrop(1,3),one,one,one,1.0_rp,rivel)
     end if
     if(rivel<=sstol_got.and.ittim>=2) then
        kfl_stead_got = 1
        call outfor(28_ip,momod(modul)%lun_outpu,' ') 
     end if

  end if
  !
  ! Formats
  !
100 format('# --| ALYA convergence  ' ,/,&
         & '# --| Columns displayed:' ,/,&
         & '# --|  1. Time Step        2. Global Iteration  3. Inner Iteration   ',/,&
         & '# --|  4. Current time     5. Drop. velocity    6. Water vol. frac   ',/,& 
         & '# --|  7. Min. velocity    8. Max. velocity     9. Min. alpha        ',/,&
         & '# --| 10. Max. alpha      11. Min. tau         12. Max. tau          ',/,&
         & '# --| 13. Min. diffusion  14. Max diffusion    15. Elapsed CPU time  ',/,&
         & '# ','          1','          2','          3',&
         &      '             4','             5','             6','             7',&
         &      '             8','             9','            10','            11',&
         &      '            12','            13') 
101 format(4x,i9,2x,i9,2x,i9,17(2x,e12.6))
200 format(&
       & '# Percentage of the number of iterations to achieve the L2 tolerance= ',e12.6,&
       & ' with a maxiumum number of iterations= ',i3,/,&
       &      '#            1','             2','             3','             4',&
       &      '             5','             6','             7','             8',&
       &      '             9','            10','            11','            12',&
       &      '            13','            14','            15','           >16')
201 format(100(2x,e12.6))
300 format('# --| ALYA convergence '  ,/,&
         & '# --| Columns displayed:' ,/,&
         & '# --| 1. Time Step         2. Global Iteration  3. Inner Iteration   ',/,&
         & '# --| 4. Inner SGS         5. Current time      6. Velocity SGS      ',//,&
         & '# ','          1','          2','          3','          4',&
         &      '             5','             6') 
301 format(2x,4(2x,i9),2(2x,e12.6))

end subroutine got_cvgunk


