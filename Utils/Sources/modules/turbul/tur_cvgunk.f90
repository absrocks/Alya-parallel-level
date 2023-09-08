subroutine tur_cvgunk(itask)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_cvgunk
  ! NAME 
  !    tur_cvgunk
  ! DESCRIPTION
  !    This routine compute the residual
  ! USES
  !    residu
  ! USED BY
  !    tur_endite
  !    tur_endste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use def_kermod,     only : kfl_logva,kfl_adj_prob,kfl_calc_sens
  use mod_ker_regularization, only : regul_k, regul_e, kfl_regularization
  use mod_outfor, only : outfor
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iturb
  integer(ip), save       :: ipass=0
  real(rp),    save       :: ritur(4),ritl2(4), turmi(2), turma(2), nutmi, nutma
  real(rp)                :: time1
  select case(itask)

  case(1_ip)
     !
     ! Check convergence of the inner iterations, always in the L2 norm:
     ! || f(n,i,j) - f(n,i,j-1)|| / ||f(n,i,j)||
     ! 
     if(kfl_algor_tur==1) then
        call residu(&
             kfl_normc_tur,one,nturb_tur,unkno,untur(:,:,1),one,iunkn_tur,&
             one,relax_tur,ritl2(iunkn_tur))     
     else
        do iunkn_tur=1,nturb_tur
           call residu(&
                kfl_normc_tur,nturb_tur,nturb_tur,unkno,untur(:,:,2),iunkn_tur,iunkn_tur,&
                one,relax_tur,ritl2(iunkn_tur))                
        end do
        iunkn_tur=nturb_tur
     end if 

!     if(iunkn_tur==nturb_tur) then

        if(kfl_normc_tur==3) then              ! Algebraic residual
           do iturb=1,nturb_tur
              ritur(iturb)=solve(iturb)%resin
           end do
        else                                   ! L2 residual
           do iturb=1,nturb_tur
              ritur(iturb)=ritl2(iturb)
           end do
        end if
        ! obtains mins and maxs values for turbul
        call tur_minmax(turmi, turma, nutmi, nutma)
        if (kfl_logva==1) then
           turmi = exp(turmi)
           turma = exp(turma)
        else if (kfl_regularization) then
           turmi(1) = regul_k(turmi(1))
           turmi(2) = regul_e(turmi(2))
           turma(1) = regul_k(turma(1))
           turma(2) = regul_e(turma(2))
        end if
        call cputim(time1)
        if((ritur(1)<cotol_tur.and.ritur(2)<cotol_tur).or.&
             (itinn(modul)>=miinn_tur)) kfl_goite_tur = 0
        if(kfl_paral<=0) then
           if(ipass==0.and.kfl_rstar/=2) then
              ipass=1
              if(nturb_tur==1) then
                 write(momod(modul)%lun_conve,100)
                 ritur(2)=0.0_rp
              else if(nturb_tur==2) then
                 write(momod(modul)%lun_conve,110)
              else 
                 write(momod(modul)%lun_conve,120)
              end if
           else
!              time1=time1-cpuit_tur
           end if
           if(nturb_tur<=2) then
              write(momod(modul)%lun_conve,101) ittim,itinn(modul),itera_tur,cutim,&
                   &                   ritur(1),ritur(2),time1, turmi(1), turma(1),&
                   &                   turmi(2), turma(2), nutmi, nutma
           else
              write(momod(modul)%lun_conve,101) ittim,itcou,itinn(modul),cutim,&
                   &                   ritur(1),ritur(2),ritur(3),ritur(4),time1
           end if
           call cputim(cpuit_tur)
           flush(momod(modul)%lun_conve)
        end if
 !    end if

  case(2_ip)
     !
     ! Check convergence of the outer iterations in the norm selected by the user:
     ! || f(n,i,*) - f(n,i-1,*)|| / ||f(n,i,*)||
     !     
     resid_tur=0.0_rp
     do iturb=1,nturb_tur
        call residu(&
             kfl_normc_tur,nturb_tur,nturb_tur,untur(:,:,1),untur(:,:,2),&
             iturb,iturb,one,1.0_rp,ritur(iturb))
        resid_tur=resid_tur+ritur(iturb)*ritur(iturb)
     end do
     resid_tur=sqrt(resid_tur)


  case(3_ip)
     !
     ! Check residual of the time evolution, always in the L2 norm:
     ! || f(n,*,*) - f(n-1,*,*)|| / ||f(n,*,*)||
     !     
     ritur=0.0_rp
     do iturb=1,nturb_tur
        call residu(&
             kfl_normc_tur,nturb_tur,nturb_tur,untur(:,:,1),untur(:,:,3),&
             iturb,iturb,one,1.0_rp,ritur(iturb))
     end do

     if(maxval(ritur(1:nturb_tur))<=sstol_tur .and. kfl_adj_prob == 0) then
        kfl_stead_tur = 1
        call outfor(28_ip,momod(modul)%lun_outpu,' ')
     elseif ((maxval(ritur(1:nturb_tur))<=sstol_tur .and. kfl_adj_prob == 1 .and. ittim >= 100) .or. kfl_calc_sens == 1) then
        kfl_stead_tur = 1
        call outfor(28_ip,momod(modul)%lun_outpu,' ')
     end if
     
  end select
  !
  ! Formats
  !
100 format('# --| ALYA Convergence '       ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   '  ,/,&
       &   '# --| 4. Current time      5. Turbulence 1       6. Nothing           ',//,&
       &   '# ','          1','          2','          3',&
       &        '             4','             5','             6')
110 format('# --| Convergence '       ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   ' ,/,&
       &   '# --| 4. Current time      5. Turbulence 1       6. Turbulence 2      ' ,/,&
       &   '# --| 7. Cputime           8. Mintur1            9. Maxtur1           ' ,/,&
       &   '# --|10. Mintur2          11. Maxtur2           12. Min_viscot        ' ,/,&
       &   '# --|13. Max_viscot                                                   ' ,//,&
       &   '# ','          1','          2','          3',&
       &        '             4','             5','             6',&
       &        '             7','             8','             9',&
       &        '            10','            11','            12',&
       &        '            13')
120 format('# --| Convergence '       ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   '  ,/,&
       &   '# --| 4. Current time      5. Turbulence 1       6. Turbulence 2      ' ,//,&
       &   '# --| 7. Turbulence 3      8. Turbulence 4       ',//,&
       &   '# ','          1','          2','          3',&
       &        '             4','             5','             6','             7',&
       &        '             8')
101 format(4x,i9,2x,i9,2x,i9,11(2x,e12.6))

end subroutine tur_cvgunk

