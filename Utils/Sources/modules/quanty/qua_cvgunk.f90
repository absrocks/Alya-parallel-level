subroutine qua_cvgunk(itask)

  !-----------------------------------------------------------------------
  !
  ! This routine performs several convergence checks for the solution
  !
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_quanty
  use mod_outfor, only : outfor
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip), save       :: ipass=0
  integer(ip)             :: ipoin
  real(rp),    save       :: cpuit_qua=0.0_rp
  real(rp)                :: ritem,time1,error,denom,numer
  real(rp)                :: dummr(2)
#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  select case(itask)

  case(1)
     !
     ! Check convergence of the inner iterations:
     ! || rho(2) - rho(1)|| / ||rho(1)||

     if(kfl_dftgs_qua==1 .or. kfl_alele_qua ==1) then

        numer = 0.0_rp
        denom = 0.0_rp
        error = 0.0_rp
        if( INOTMASTER ) then
           do ipoin = 1,npoin
              numer = numer + (rhoon(ipoin,2)-rhoon(ipoin,1))*(rhoon(ipoin,2)-rhoon(ipoin,1))
              denom = denom +  rhoon(ipoin,1)*rhoon(ipoin,1)
           enddo
        end if
        dummr(1) = numer
        dummr(2) = denom
        call pararr('SUM',0_ip,2_ip,dummr)
        error = sqrt(dummr(1)/dummr(2))

        write(6,*) '###### ##### DFT convergence = = ',error      
        if( error < cotol_qua .or. itinn(modul) >= miinn_qua ) kfl_goite_qua = 0

        !
        ! Compute min and max of the density
        !
        call minmax(one,npoin,zero,rhoon(:,2),rhomin_qua,rhomax_qua)

     else

        rhomin_qua    = 0.0_rp
        rhomax_qua    = 0.0_rp
        error         = 0.0_rp
        kfl_goite_qua = 0

     endif

     !
     ! Write convergence
     !
     if( INOTSLAVE ) then
        call cputim(time1)
        if(ipass==0) then
           time1=time1-cpu_initi
        else
           time1=time1-cpuit_qua
        end if
        if(ipass==0.and.kfl_rstar/=2) write(lun_conve_qua,100)

        write(lun_conve_qua,101) ittim,itcou,itinn(modul),cutim,error,&
             rhomin_qua,rhomax_qua,time1

        call cputim(cpuit_qua)
        flush(lun_conve_qua)
        ipass=1

     end if

  case(2)
     !
     ! Check convergence of the outer iterations: Todavia no hay tiempo asi que no hace nada

  case(3)
     !
     ! Check residual of the time iterations:
     ! || T(n,*,*) - T(n-1,*,*)|| / ||T(n,*,*)||
     !
     !call residu(kfl_normc_qua,one,one,rhoon(1,1),rhoon(1,3),one,one,one,1.0_rp,ritem)

     if(ritem<=sstol_qua) then
        kfl_stead_qua = 1
        call outfor(28_ip,lun_outpu_qua,' ')
     end if
  end select
  !
  ! Formats
  !
100 format('# --| ALYA Convergence '  ,/,&
       &   '# --| Columns displayed:' ,/,&
       &   '# --| 1. Time step         2. Global Iteration   3. Inner Iteration   '  ,/,&
       &   '# --| 4. Current time      5. Solution           6. Solution SGS    ' ,/,& 
       &   '# --| 7. Min. rho          8. Max. rho           9. Elapsed CPU Time  ' ,//,&
       &   '$ ','          1','          2','          3',&
       &        '             4','             5','             6','             7',&
       &        '             8','             9')
101 format(4x,i9,2x,i9,2x,i9,11(2x,e12.6))

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine qua_cvgunk

