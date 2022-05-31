subroutine nsa_cvgunk(itask)
!-----------------------------------------------------------------------
!****f* Nastal/nsa_cvgunk
! NAME
!    nsa_cvgunk
! DESCRIPTION
!    This routine performs several convergence checks
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_kermod
  use      def_domain
  use      def_nastal
  use      def_solver, only : iters,resi1

  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass=0
  integer(ip)             :: idofn,ipoin
  real(rp)                :: rinsa(4),cpnew,cpdif,vauxi
  real(rp),save           :: cpold=0.0
  real(rp),save           :: cpcum=0.0

  select case(itask)

  case(1)
     !
     ! Check convergence of the inner iterations:
     ! || u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||
     !
     rinsa=0.0_rp
     call residu(kfl_normc_nsa,ndofn_nsa,ndime,unkno,&
          umome,one,one,ndime,1.0_rp,rinsa(1))
     call residu(kfl_normc_nsa,ndofn_nsa,one,unkno,&
          densi,ndime+1,one,one,1.0_rp,rinsa(2))
     call residu(kfl_normc_nsa,ndofn_nsa,one,unkno,&
          energ,ndime+2,one,one,1.0_rp,rinsa(3))     

!     rinsa = rinsa * dtinv_nsa  ! convergence will be checked against DU / Dt
     rinsa(4) = rinsa(1) + rinsa(2) + rinsa(3)

     if(itinn(modul) == 1) then
        do idofn= 1,4
           resin_first_nsa(idofn) = 1.0_rp
           if (rinsa(idofn) > zensa) then ! to avoid division by zero
              resin_first_nsa(idofn) = rinsa(idofn)
           end if
        end do
     end if

     if (kfl_refre_nsa(itask) == 0) then
        do idofn=1,4
           if (rinsa(idofn) > zensa) then ! to avoid division by zero
              resin_nsa(idofn) = rinsa(idofn)
           else
              resin_nsa(idofn) = zensa
           end if
        end do
        kfl_refre_nsa(itask) = 1
     end if
!!$     rinsa(1)= rinsa(1)/resin_nsa(1)
!!$     rinsa(2)= rinsa(2)/resin_nsa(2)
!!$     rinsa(3)= rinsa(3)/resin_nsa(3)
!!$     rinsa(4)= rinsa(4)/resin_nsa(4)

     !
     ! check stopping condition for subiterations
     !
     if (kfl_adres_nsa == 0) then
        if(rinsa(4)<cotol_nsa.or.itinn(modul)>=miinn_nsa) kfl_goite_nsa = 0
     else
        if(itinn(modul)>=miinn_nsa) then
           kfl_goite_nsa = 0
        else
           if (rinsa(4) < cotol_nsa) kfl_goite_nsa = 0                     ! firstly, absolute check
           if (rinsa(4)/resin_first_nsa(4) < corat_nsa) kfl_goite_nsa = 0  ! then, relative check
        end if
     end if
     call minmax(one,npoin,zero,densi,vamxm_nsa(1,1),vamxm_nsa(2,1))
     call minmax(one,npoin,zero,press,vamxm_nsa(1,2),vamxm_nsa(2,2))
     call minmax(one,npoin,zero,tempe,vamxm_nsa(1,3),vamxm_nsa(2,3))
     call minmax(one,npoin,zero,vmach_nsa,vamxm_nsa(1,4),vamxm_nsa(2,4))
     
     !
     ! Sum up slaves contribution for extensive quantities
     !
     call pararr('SUM',0_ip,1_ip,envol_nsa)
     call pararr('SUM',0_ip,1_ip,devol_nsa)
     
     if( INOTSLAVE .and. miinn_nsa > 1) then   ! only done when subiterations are being done
        if (writeliveinfo()) then
           if(ipass==0.and.kfl_rstar/=2) then
              ipass=1
              write(momod(modul)%lun_conve,200)
              write(lun_maxmi_nsa,300)
              write(lun_force_nsa,400)
           end if
        end if
        !        if (ittim == 1) cpold = 0.0_rp
        call cputim(cpnew)
        cpdif = cpnew - cpold
        cpold = cpnew
        cpcum = cpcum + cpdif

        last_iters_nsa = solve_sol(1)%iters
!
!       Write convergence files
!
        if (writeliveinfo()) then
           write(momod(modul)%lun_conve,1000) ittim,itcou,itinn(modul),last_iters_nsa,cutim,&
                rinsa(4),rinsa(1),rinsa(2),rinsa(3),1.0_rp/dtinv,1.0_rp/dtinv_nsa,envol_nsa,&
                devol_nsa,cpdif,cpcum,resi1,safet_nsa,uinlet_nsa
           flush(momod(modul)%lun_conve)
           write(lun_maxmi_nsa,1000) ittim,itcou,itinn(modul),iters,&
                vamxm_nsa(1,1),vamxm_nsa(2,1),vamxm_nsa(1,2),vamxm_nsa(2,2),&
                vamxm_nsa(1,3),vamxm_nsa(2,3),vamxm_nsa(1,4),vamxm_nsa(2,4),&
                1.0_rp/dtinv_nsa,dtcri_nsa(1),dtcri_nsa(2),dtcri_nsa(3),&
                dtmax_nsa(1),dtmax_nsa(2),dtmax_nsa(3),resi1
           flush(lun_maxmi_nsa)
        end if

     end if

!
!       Write aerodynamic body coefficient forces file
!
     call nsa_aebody(3,1,1,1,1)


  case(2)
     !
     ! Check convergence of the outer iterations:
     ! || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||
     !
     call residu(kfl_normc_nsa,ndime,ndime,umome,umome(1,1,ITER_AUX),&
          one,one,ndime,1.0_rp,resid_nsa(1))
     call residu(kfl_normc_nsa,one,one,densi,densi(1,ITER_AUX),&
          one,one,one,1.0_rp,resid_nsa(2))
     call residu(kfl_normc_nsa,one,one,energ,energ(1,ITER_AUX),&
          one,one,one,1.0_rp,resid_nsa(3))
     resid_nsa(4) = resid_nsa(1) + resid_nsa(2) + resid_nsa(3)

  case(3)
     !
     ! Check residual of the time evolution:
     ! || u(n,*,*) - u(n-1,*,*)|| / ||u(n,*,*)||
     !
     rinsa=0.0_rp
     call residu(kfl_normc_nsa,ndime,ndime,umome(1,1,ITER_K),umome(1,1,TIME_N),&
          one,one,ndime,1.0_rp,rinsa(1))
     call residu(kfl_normc_nsa,  one,  one,densi(1,  ITER_K),densi(1,  TIME_N),&
          one,one,one,1.0_rp,rinsa(2))
     call residu(kfl_normc_nsa,  one,  one,energ(1,  ITER_K),energ(1,  TIME_N),&
          one,one,one,1.0_rp,rinsa(3))
     rinsa(4) = rinsa(1) + rinsa(2) + rinsa(3)

     call nsa_corens(ITER_K,TIME_N)

     if (kfl_refre_nsa(itask) == 0) then
        do idofn=1,4
           if (rinsa(idofn) > zensa) then ! to avoid division by zero
              resou_nsa(idofn) = rinsa(idofn)
           else
              resou_nsa(idofn) = zensa
           end if
        end do
        kfl_refre_nsa(itask) = 1
     end if
!!$     rinsa(1)= rinsa(1)/resou_nsa(1)
!!$     rinsa(2)= rinsa(2)/resou_nsa(2)
!!$     rinsa(3)= rinsa(3)/resou_nsa(3)
!!$     rinsa(4)= rinsa(4)/resou_nsa(4)

     if( INOTSLAVE ) then
        if (writeliveinfo()) then
           if(ipass==0.and.kfl_rstar/=2) then
              ipass=1
              write(momod(modul)%lun_conve,200)
              write(lun_maxmi_nsa,300)
              write(lun_force_nsa,400)
           end if
        end if
!        if (ittim == 1) cpold = 0.0_rp
        call cputim(cpnew)
        cpdif = cpnew - cpold
        cpold = cpnew
        cpcum = cpcum + cpdif

        last_iters_nsa = solve_sol(1)%iters
!
!       Write convergence files
!

        if (writeliveinfo()) then
           write(momod(modul)%lun_conve,1000) ittim,itcou,itinn(modul),last_iters_nsa,cutim,&
                rinsa(4),rinsa(1),rinsa(2),rinsa(3),1.0_rp/dtinv,1.0_rp/dtinv_nsa,envol_nsa,&
                devol_nsa,cpdif,cpcum,resi1,safet_nsa,uinlet_nsa
           flush(momod(modul)%lun_conve)
           write(lun_maxmi_nsa,1000) ittim,itcou,itinn(modul),iters,&
                vamxm_nsa(1,1),vamxm_nsa(2,1),vamxm_nsa(1,2),vamxm_nsa(2,2),&
                vamxm_nsa(1,3),vamxm_nsa(2,3),vamxm_nsa(1,4),vamxm_nsa(2,4),&
                1.0_rp/dtinv_nsa,dtcri_nsa(1),dtcri_nsa(2),dtcri_nsa(3),&
                dtmax_nsa(1),dtmax_nsa(2),dtmax_nsa(3),resi1
           flush(lun_maxmi_nsa)
        end if

     end if


  end select
!
! Formats
!
102 format('# >>>  PROBLEM BECOMES STATIONARY AT TIME STEP ',i5)

200 format('# ',' --|  Nastal Convergence File '       ,/,&
         & '# ',' --|  Columns displayed:' ,/,&
         & '# ',' --|  1. Time Step     2. Global Iteration     3. Sub Iteration   4. Solver iterations'  ,/,&
         & '# ',' --|  5. Current time  ',/,&
         & '# ',' --|  6. Total Res.    7. Mom. Res.            8. Cont. Res.        9. Ener. Res.   ' ,/,&
         & '# ',' --| 10. Delta t      11. Pseudo Delta t      12. Vol Energy       13. Vol. Density ' ,/,&   
         & '# ',' --| 14. CPU time     15. Cumulative CPU time 16. Solver resi      17. Safety factor   18. Uinlet ') 

300 format('# ',' --|  MinMax File '       ,/,&
         & '# ',' --|  Columns displayed:' ,/,&
         & '# ',' --|  1. Time Step   2. Global Iteration   3. Inner Iteration'  ,/,&
         & '# ',' --|  4. Dens Min    5. Dens Max   6. Press Min   7. Press Max' ,/,&
         & '# ',' --|  8. Temp Min    9. Temp Max  10. Mach  Min  11. Mach  Max' ,/,&
         & '# ',' --| 12. Delta t    13. Dt Min Momen  14. Dt Min Densi  15. Dt Min Ener' ,/,&
         & '# ',' --| 16. Dt Max Momen  17. Dt Max Densi  18. Dt Max Ener')

400 format('# ',' --|  Aerodynamic body force coefficients File '       ,/,&
         & '# ',' --|  Columns displayed:' ,/,&
         & '# ',' --|  1. Time Step   2. Global Iteration   3. Inner Iteration   4. Current time'  ,/,&
         & '# ',' --|  5-7.   Press CL/D (total area)  X,Y,Z'  ,/,&
         & '# ',' --|  8-10.  Visco CL/D (total area)  X,Y,Z'  ,/,&
         & '# ',' --|  11-13. Press CX (frontal area)  X,Y,Z'  ,/,&
         & '# ',' --|  14-16. Visco CX (frontal area)  X,Y,Z'  ,/,&
         & '# ',' --|  17-19. Press force  X,Y,Z'  ,/,&
         & '# ',' --|  20-22. Visco force  X,Y,Z'  ,/,&
         & '# ',' --|  ') 

1000 format(4x,i9,2x,i9,2x,i9,2x,i9,25(2x,e12.6))


end subroutine nsa_cvgunk
