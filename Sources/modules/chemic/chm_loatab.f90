subroutine chm_loatab()
    !-----------------------------------------------------------------------
    !****f* chemic/chm_loatab
    ! NAME 
    !    chm_loatab
    ! DESCRIPTION
    !    Load thermochemical database for CFI combustion model at start
    ! USES
    ! USED BY
    !    chm_reaphy
    !   
    !***
    !-----------------------------------------------------------------------
    use def_inpout, only      : param,words
    use def_kintyp, only      : ip,rp
    use def_master, only      : table_cfi,speci
    use mod_ecoute, only      : ecoute
    implicit none
    integer(ip)               :: ii,npara,ipara,icoef,jcoef
    real(rp)                  :: rpara

    if (words(1) == 'NVALU') then
       do ii=1,table_cfi(1)%nclas
          table_cfi(1)%nvcfi(ii) = int(param(ii))
       end do
    end if
    if (table_cfi(1)%ndcfi == 3_ip .or. table_cfi(1)%ndcfi == 5_ip) table_cfi(1)%nvcfi(table_cfi(1)%ndcfi) = int(param(5))

    if ((table_cfi(1)%nclas == 2_ip .and. (table_cfi(1)%nvcfi(1) == 0_ip .or. table_cfi(1)%nvcfi(2) == 0_ip)) .or. &
       (table_cfi(1)%nclas == 4_ip .and. (table_cfi(1)%nvcfi(1) == 0_ip .or. table_cfi(1)%nvcfi(2) == 0_ip .or. &
       table_cfi(1)%nvcfi(3) == 0_ip .or. table_cfi(1)%nvcfi(4) == 0_ip)) .or. &
       ((table_cfi(1)%ndcfi == 3_ip .or. table_cfi(1)%ndcfi == 5_ip) .and. table_cfi(1)%nvcfi(table_cfi(1)%ndcfi) == 0_ip) ) &
       call runend('CHEMIC REAPHY: Wrong table setup for use with the CFI combustion model')
    !
    ! Get full size of database: % nrcfi = c * c_var * f * f_var [database entries]
    !
    table_cfi(1)%nrcfi = 1_ip
    do ii = 1,table_cfi(1)%ndcfi
      table_cfi(1)%nrcfi = table_cfi(1)%nrcfi * table_cfi(1)%nvcfi(ii)
    end do
    !
    ! Read reaction progress variable, variance, mixture fraction and variance spacing: c, c_var, f & f_var
    ! 
    call chm_memphy(5_ip)
    
    call ecoute('chm_reaphy')
    do ii = 1,table_cfi(1)%nclas
      if (words(1) == speci(ii)%name) then
        call ecoute('chm_reaphy')
        rpara = real(table_cfi(1)%nvcfi(ii)) / 15.0_rp
        npara = ceiling(rpara)
        ipara = 1_ip
        do icoef = 1,npara
          do jcoef = 1,15
            table_cfi(1)%ivcfi(ii,ipara) = param(jcoef)
            if (ipara == table_cfi(1)%nvcfi(ii)) exit
            ipara = ipara + 1_ip
          end do
          call ecoute('chm_reaphy')
        end do
      else 
         call runend('CHEMIC REAPHY: Wrong table setup for use with the CFI combustion model')
      end if
    end do
    !
    ! Read normalized enthalpy spacing: i 
    !
    if (table_cfi(1)%ndcfi == 3_ip .or. table_cfi(1)%ndcfi == 5_ip) then
      if (words(1) == 'IMEAN') then
        call ecoute('chm_reaphy')
        rpara = real(table_cfi(1)%nvcfi(table_cfi(1)%ndcfi)) / 15.0_rp
        npara = ceiling(rpara)
        ipara = 1_ip
        do icoef = 1,npara
          do jcoef = 1,15
            table_cfi(1)%ivcfi(table_cfi(1)%ndcfi,ipara) = param(jcoef)
            if (ipara == table_cfi(1)%nvcfi(table_cfi(1)%ndcfi)) exit
            ipara = ipara + 1_ip
          end do
        call ecoute('chm_reaphy')
        end do
      else
         call runend('CHEMIC REAPHY: Wrong table setup for use with the CFI combustion model')
      end if
    end if
    !
    ! Read additional properties for non-premixed flames 
    !
    if (table_cfi(1)%nclas == 4_ip) then
      !
      ! Take min and max mixture fraction: f_min, f_max and extract minimum value from all 
      !
      table_cfi(1)%fmima(1) = table_cfi(1)%ivcfi(3,1)
      table_cfi(1)%fmima(2) = table_cfi(1)%ivcfi(3,table_cfi(1)%nvcfi(3))
      do ipara = 1,table_cfi(1)%nvcfi(3) 
        table_cfi(1)%ivcfi(3,ipara) = table_cfi(1)%ivcfi(3,ipara) - table_cfi(1)%fmima(1)
      end do
      !
      ! Read inlet mixture properties 
      !
      call ecoute('chm_reaphy')
      do icoef = 1,table_cfi(1)%nfcfi
        table_cfi(1)%inval(1,icoef) = param(icoef)
      end do
      call ecoute('chm_reaphy')
      call ecoute('chm_reaphy')
      do icoef = 1,table_cfi(1)%nfcfi
        table_cfi(1)%inval(2,icoef) = param(icoef)
      end do
      call ecoute('chm_reaphy')
    end if
    !
    ! Read max and min enthalpy levels: h_max & h_min
    !
    if (table_cfi(1)%ndcfi == 3_ip .or. table_cfi(1)%ndcfi == 5_ip) then
      call ecoute('chm_reaphy')
      table_cfi(1)%imima(1) = param(1)
      table_cfi(1)%imima(2) = param(2)
      call ecoute('chm_reaphy')
    end if
    !
    ! Read full database 
    !
    call ecoute('chm_reaphy')
    do ii=1,table_cfi(1)%nrcfi
       do ipara=1,10
         table_cfi(1)%table(ii,ipara) = param(ipara+5)
       end do
       call ecoute('chm_reaphy')
       do ipara=1,table_cfi(1)%nccfi-10
         table_cfi(1)%table(ii,ipara+10) = param(ipara)
       end do
       call ecoute('chm_reaphy')
    end do
    
end subroutine chm_loatab
