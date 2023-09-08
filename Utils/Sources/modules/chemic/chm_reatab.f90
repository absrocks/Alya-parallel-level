subroutine chm_reatab()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_reatab
  ! NAME 
  !    chm_reatab
  ! DESCRIPTION
  !    Read table properties for CFI combustion model
  ! USES
  ! USED BY
  !    chm_iniunk: initialize flow field
  !    chm_endite: update table properties to start doiter in temper with updated coefficients 
  !    chm_endste: properties available for all modules at the end of time-step
  !                (there is no update of properties during chemic iteration)   
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      : ip,rp
  use def_master, only      : wmean,conce,inotmaster,sphek,visck,condk,tempe,densi, &
                              sphec,therm,kfl_paral,table_cfi,encfi,lescl,massk
  use def_domain, only      : npoin,ltypb,nnode,lnodb,nboun
  use def_chemic, only      : tempe_chm,kfl_radia_chm,rspec_chm,kfl_wallc_chm,yscale_chm, &
                              cvar_chm
  implicit none
  integer(ip)               :: ipoin,iclas,ivalu,trang
  integer(ip)               :: iboun,pblty,inodb,pnodb
  real(rp)                  :: retva(table_cfi(1)%nccfi)      ! Values read in from CFI table: 
                                                              ! (1) S_c (2) Wmean (3) Lambda (4) Mu (5-10) cp_low (11-16) cp_high (17) S_c*c
  real(rp)                  :: cploc                          ! Cp calculated from table using local flow field temperature
  real(rp)                  :: teloc                          ! Local flow field temperature
  real(rp)                  :: tab_conce(table_cfi(1)%ndcfi)  ! Concentration for each variable at each node

  if (INOTMASTER) then

    yscale_chm = 0.0_rp
    cvar_chm   = 0.0_rp

    do ipoin = 1,npoin
      !
      ! Initialization
      !
      do iclas = 1,table_cfi(1)%nclas
        massk(ipoin,iclas) = 0.0_rp
        tab_conce(iclas)   = min(1.0_rp,max(0.0_rp,(conce(ipoin,iclas,1))))
      end do
      !
      ! compute normalized enthalpy in case of non-adiabatic combustion
      !
      if(table_cfi(1)%ndcfi == 3 .or. table_cfi(1)%ndcfi == 5) then 
        tab_conce(table_cfi(1)%ndcfi) = ((table_cfi(1)%imima(1) - therm(ipoin,1)) / &
                                        (table_cfi(1)%imima(1) - table_cfi(1)%imima(2)))
        encfi(ipoin) = tab_conce(table_cfi(1)%ndcfi)
      end if
      !
      ! get values from table 
      !    
      if (table_cfi(1)%nclas == 2) then
        call chm_getval(1_ip,tab_conce,retva,yscale_chm(ipoin),cvar_chm(ipoin))
      else
        if (tab_conce(3) > table_cfi(1)%fmima(1) .and.  tab_conce(3) < table_cfi(1)%fmima(2)) then
          call chm_getval(1_ip,tab_conce,retva,yscale_chm(ipoin),cvar_chm(ipoin))
        else
          call chm_getval(3_ip,tab_conce,retva,yscale_chm(ipoin),cvar_chm(ipoin))
        end if
      end if
      !
      ! clip of the temperature for the use in the polynomial function
      ! and set temperature range
      !   
      if (tempe_chm(ipoin) < 200.0_rp) then
        teloc = 200.0_rp
      else if (tempe_chm(ipoin) > 3000.0_rp) then
        teloc = 3000.0_rp
      else
        teloc = tempe_chm(ipoin)
      end if
      if (teloc < 1000.0_rp) then
        trang = 1_ip
      else
        trang = 2_ip
      end if
      !
      ! update properties 
      !
      massk(ipoin,1) = retva(1)
      wmean(ipoin,1) = retva(2)

      do iclas = 1,table_cfi(1)%nclas
        condk(ipoin,iclas) = retva(3)
        visck(ipoin,iclas) = retva(4)
      end do
      do ivalu = 1,6 
        sphec(ipoin,ivalu,1) = retva(4+ivalu)
        sphec(ipoin,ivalu,2) = retva(10+ivalu)
      end do
      lescl(ipoin) = retva(17)  ! filter{w_c * c}

      if (kfl_radia_chm > 0) then
        rspec_chm(1,ipoin) = retva(18)
        rspec_chm(2,ipoin) = retva(19)
      end if
      
    end do
    !
    ! Imposing reaction source term to zero at walls 
    !
    if (kfl_wallc_chm == 1_ip ) then
     
      boundaries: do iboun = 1,nboun
         pblty = ltypb(iboun)
         pnodb = nnode(pblty)
         do inodb = 1,pnodb
            ipoin = lnodb(inodb,iboun)
            massk(ipoin,1_ip) = 0.0_rp
         end do
      end do boundaries

    endif

  endif
  
end subroutine chm_reatab
