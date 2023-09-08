subroutine rad_cvgunk(itask)

  !-----------------------------------------------------------------------
  !
  ! This routine performs several convergence checks for the radiation
  !
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_radiat
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass=0
  real(rp),    save       :: cpuit_rad=0.0_rp
  real(rp)                :: ritem,time1,ritl2
  integer                 :: ipoin

#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  select case(itask)

  case(1)
     !
     ! Check convergence of the inner iterations:
     ! || G(n,i,j) - G(n,i,j-1)|| / ||G(n,i,j)||
     !
     call residu(kfl_normc_rad,one,one,unkno,radav_rad(:,1),one,one,one,1.0_rp,ritl2)      
     if(kfl_normc_rad==3) then  
        ritem=solve(1)%resin       ! Algebraic residual 
     else
        ritem=ritl2                ! L2 residual
     end if
!!$     if (INOTMASTER) then
!!$        ritem=0.0_rp
!!$        do ipoin=1,npoin
!!$           ritem=ritem+(unkno(ipoin)-radav_rad(ipoin,1))**2
!!$        enddo
!!$        ritem = sqrt(ritem)
!!$        print *,' RES:',ritem
!!$     else
!!$        print *,' yo master :',ritem
!!$     endif
     if((ritem<cotol_rad).or.(itinn(modul)>=miinn_rad)) kfl_goite_rad = 0   !! We have achieved convergence or exceeded max # iterations
     !
     ! Compute min and max of the radiation average intensity
     !
     call minmax(one,npoin,zero,unkno,ramin_rad,ramax_rad)
     !
     ! Compute SGS residual
     !
     if(kfl_sgsno_rad==1) then
        if(kfl_paral>=0) then
           nparr    =  2
           parre    => resgs_rad
           call Parall(9_ip)
        end if
        if(resgs_rad(2)>0.0_rp) resgs_rad(1)=sqrt(resgs_rad(1)/resgs_rad(2))
     end if
     !
     ! Write convergence
     !
     if(kfl_paral<=0) then
        call cputim(time1)
        if(ipass==0) then
           time1=time1-cpu_initi
        else
           time1=time1-cpuit_rad
        end if
        if(ipass==0.and.kfl_rstar/=2) write(momod(modul)%lun_conve,100)
        write(momod(modul)%lun_conve,101) ittim,itcou,itinn(modul),cutim,ritem,resgs_rad(1),&
             ramin_rad,ramax_rad,time1
        call cputim(cpuit_rad)
        flush(momod(modul)%lun_conve)
        ipass=1

     end if

  case(2)
     !
     ! Check convergence of the outer iterations:
     ! || G(n,i,*) - G(n,i-1,*)|| / ||G(n,i,*)||
     !
     call residu(kfl_normc_rad,one,one,radav_rad(1,1),radav_rad(1,2),one,one,one,1.0_rp,resid_rad)

  case(3)
     !
     ! Check residual of the time iterations:
     ! || T(n,*,*) - T(n-1,*,*)|| / ||T(n,*,*)||
     !
     call residu(kfl_normc_rad,one,one,radav_rad(1,1),radav_rad(1,3),one,one,one,1.0_rp,ritem)
  end select
  !
  ! Formats
  !
100 format('# --| ALYA Convergence '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   '  ,/,&
       &   '# --| 4. Current time      5. Radiation Av. Int.        6. Radiation SGS    ' ,/,& 
       &   '# --| 7. Min. rad. intensity  8. Max. rad. intensity   9. Elapsed CPU Time  ' ,//,&
       &   '$ ','          1','          2','          3',&
       &        '             4','             5','             6','             7',&
       &        '             8','             9')
101 format(4x,i9,2x,i9,2x,i9,11(2x,e12.6))

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine rad_cvgunk

