subroutine tem_reatab()
  !-----------------------------------------------------------------------
  ! NAME 
  !    tem_reatab
  ! DESCRIPTION
  !    Read table properties for CFI combustion model
  ! USES
  ! USED BY
  !   
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      : ip,rp
  use def_master, only      : wmean,conce,inotmaster,sphek,visck,condk,tempe,densi, &
                              sphec,therm,kfl_paral,table_cfi,lescl,encfi,massk
  use def_domain, only      : npoin
  implicit none
  integer(ip)               :: ipoin,iclas,inodb,ivalu
  real(rp)                  :: retva(table_cfi(1)%nccfi)           ! Values read in from CFI table: 
                                                          ! (1) S_c (2) Wmean (3) Lambda (4) Mu (5-10) cp_low (11-16) cp_high (17) S_c*c
  real(rp)                  :: tab_conce(table_cfi(1)%ndcfi)       ! Concentration for each variable at each node

  if (INOTMASTER) then

!  write(*,*) 'read table'

    do ipoin = 1,npoin
      !
      ! Initialization
      !
      do iclas = 1,table_cfi(1)%nclas
        massk(ipoin,iclas) = 0.0_rp
        tab_conce(iclas)       = conce(ipoin,iclas,1)
      end do
      !
      ! compute normalized enthalpy in case of non-adiabatic combustion
      !
      if(table_cfi(1)%ndcfi == 3 .or. table_cfi(1)%ndcfi == 5) then 
        tab_conce(table_cfi(1)%ndcfi) = 1.0_rp - ((table_cfi(1)%imima(1) - therm(ipoin,1)) / &
                                         (table_cfi(1)%imima(1) - table_cfi(1)%imima(2)))
        encfi(ipoin) = tab_conce(table_cfi(1)%ndcfi)
      end if
      !
      ! get values from table 
      !    
      if (table_cfi(1)%nclas == 2) then
        call tem_getval(1_ip,tab_conce,retva)
      else
        if (tab_conce(3) > table_cfi(1)%fmima(1) .and.  tab_conce(3) < table_cfi(1)%fmima(2)) then
          call tem_getval(1_ip,tab_conce,retva)
        else
!          call tem_getval(2_ip,tab_conce,retva)
          call tem_getval(3_ip,tab_conce,retva)
        end if
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

    end do

  endif

end subroutine tem_reatab
