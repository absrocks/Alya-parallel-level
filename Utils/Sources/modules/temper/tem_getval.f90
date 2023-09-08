subroutine tem_getval(itask,inval_global,reval) 
  !-----------------------------------------------------------------------
  !**** 
  ! NAME 
  !    tem_getval
  ! DESCRIPTION
  !    Returns properties from table for a given combination of reaction 
  !    progress variable, mixture fraction and their variances (stored
  !    in conce)
  ! USES
  ! USED BY
  !   
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      : ip,rp
  use def_master, only      : wmean,conce,table_cfi
  implicit none
  integer(ip)               :: iclas,iedge,icube,jclas,ncube,itask
  integer(ip)               :: ivalu(table_cfi(1)%ndcfi),vcube(table_cfi(1)%ndcfi,16),numpo(table_cfi(1)%ndcfi+1)
  integer(ip)               :: order(table_cfi(1)%ndcfi),itemp,idumi,dimdx(table_cfi(1)%ndcfi+1,table_cfi(1)%ndcfi)
  integer(ip)               :: nline(16),nmult,nvalu(table_cfi(1)%ndcfi+1)
  real(rp)                  :: local_mf
  real(rp)                  :: dxval(table_cfi(1)%ndcfi),acval(table_cfi(1)%ndcfi),rtemp
  real(rp)                  :: fuval(table_cfi(1)%ndcfi+1,table_cfi(1)%nccfi)
  real(rp)                  :: inval(table_cfi(1)%ndcfi)
  real(rp),     intent(in)  :: inval_global(table_cfi(1)%ndcfi)
  real(rp),     intent(out) :: reval(table_cfi(1)%nccfi)

  select case (itask)

  !
  ! premixed or mean of mixture fraction is inside the flamability limit
  !
  case(1_ip)

    ncube = int(2.0_rp**real(table_cfi(1)%ndcfi,rp),ip)
  
    do iclas = 1, table_cfi(1)%ndcfi
      nvalu(iclas) =  table_cfi(1)%nvcfi(iclas)
      inval(iclas) = inval_global(iclas)
    end do
    nvalu(table_cfi(1)%ndcfi+1) = 1_ip

    if (table_cfi(1)%nclas == 4_ip) inval(3) = inval(3) - table_cfi(1)%fmima(1)
    do iclas = 1,table_cfi(1)%ndcfi
      !
      ! clip values to within allowable range 
      !
      acval(iclas) = min(max(inval(iclas),0.0_rp),table_cfi(1)%ivcfi(iclas,table_cfi(1)%nvcfi(iclas)))
      !
      ! get edge points and rib lengths of hypercube
      !
      do iedge = 1,table_cfi(1)%nvcfi(iclas)-1
        if ((table_cfi(1)%ivcfi(iclas,iedge+1) >= acval(iclas)) .and. &
        (table_cfi(1)%ivcfi(iclas,iedge) <= acval(iclas))) then
          ivalu(iclas) = iedge
          dxval(iclas) = (acval(iclas) - table_cfi(1)%ivcfi(iclas,iedge)) &
          / (table_cfi(1)%ivcfi(iclas,iedge+1) - table_cfi(1)%ivcfi(iclas,iedge))
          exit
        end if
      end do
      do icube = 1,ncube         
         rtemp = (-1.0_rp)**( real(icube-1_ip,rp)/(2.0_rp**real(iclas-1,rp)) + 1.0_rp )
        if (rtemp < 0.0_rp) then
          vcube(iclas,icube) = ivalu(iclas)
        else if (rtemp > 0.0_rp) then
          vcube(iclas,icube) = ivalu(iclas) + 1_ip
        end if
      end do
    end do
    !
    ! order distance in magnitude from high to low
    !
    do iclas = 1,table_cfi(1)%ndcfi
      order(iclas) = iclas
    end do
    iclas = 1_ip
    do while (iclas < table_cfi(1)%ndcfi)
      rtemp = 0.0_rp
      do jclas = iclas,table_cfi(1)%ndcfi
        if (dxval(order(jclas)) >= rtemp) then
          rtemp = dxval(order(jclas))
          itemp = jclas
        end if
      end do
      idumi = order(iclas)
      order(iclas) = order(itemp)
      order(itemp) = idumi
      iclas = iclas + 1_ip
    end do
    !
    ! Get edge points of line elements DX over which is interpolated
    !
    do iclas = 1,table_cfi(1)%ndcfi
      dimdx(1,iclas) = 0_ip
    end do
    do iclas = 2,table_cfi(1)%ndcfi+1
      do jclas = 1,table_cfi(1)%ndcfi
        dimdx(iclas,jclas) = dimdx(iclas-1,jclas)
        if (order(iclas-1) == jclas) dimdx(iclas,jclas) = 1_ip
      end do
    end do
    !
    ! Get point numbers
    !
    numpo(1) = 1_ip
    do iclas = 1,table_cfi(1)%ndcfi
      numpo(iclas+1) = 1_ip
      do jclas = 1,table_cfi(1)%ndcfi
        numpo(iclas+1) = numpo(iclas+1) &
                       + int(real(dimdx(iclas+1,jclas),rp)*2.0_rp**real(jclas-1,rp))
      end do
    end do
    !
    ! Obtain line numbers of edge points
    ! 
    do iclas = 1,table_cfi(1)%ndcfi+1
      nline(numpo(iclas)) = 1_ip
      do jclas = 1,table_cfi(1)%ndcfi
        nmult = 1_ip
        do iedge = jclas,table_cfi(1)%ndcfi
          nmult = nmult*nvalu(iedge+1)
        end do
        nline(numpo(iclas)) = nline(numpo(iclas)) &
                            + (vcube(jclas,numpo(iclas))-1_ip)*nmult
      end do
    end do
    !
    ! Get Function values of corresponding edge points
    !
    do iclas = 1,table_cfi(1)%ndcfi+1
      do jclas = 1,table_cfi(1)%nccfi
        fuval(iclas,jclas) = table_cfi(1)%table(nline(numpo(iclas)),jclas)
      end do
    end do
    !
    ! Get interpolated value
    !
    do jclas = 1,table_cfi(1)%nccfi
      reval(jclas) = fuval(1,jclas)
    end do
    do iclas = 1,table_cfi(1)%ndcfi 
      do jclas = 1,table_cfi(1)%nccfi
        reval(jclas) = reval(jclas) + (fuval(iclas+1,jclas) - fuval(iclas,jclas))*dxval(order(iclas))
      end do
    end do

  !
  ! non-premixed and mean mixture fraction is outside the flamability limit
  !  
  case(2_ip)

    !
    ! clip values to within allowable range 
    !
    acval(3) = min(max(inval_global(3),0.0_rp),1.0_rp)

    reval(1) = 0.0_rp
    reval(17) = 0.0_rp

    do iclas = 1,table_cfi(1)%nfcfi
      reval(iclas+1) = acval(3) * table_cfi(1)%inval(2,iclas) + (1-acval(3)) * table_cfi(1)%inval(1,iclas)
    end do 

  !
  ! non-premixed outside flamability limit
  !
  case(3_ip)

    ncube = int(2.0_rp**real(table_cfi(1)%ndcfi,rp),ip)

    if (inval_global(3) > table_cfi(1)%fmima(1)) then

    else if(inval_global(3) < table_cfi(1)%fmima(2)) then

    end if

    do iclas = 1, table_cfi(1)%ndcfi
      nvalu(iclas) =  table_cfi(1)%nvcfi(iclas)
      inval(iclas) = inval_global(iclas)
    end do
    if (inval(3) < table_cfi(1)%fmima(1)) then
      local_mf = inval(3)
      inval(3) = table_cfi(1)%fmima(1)
    else if(inval(3) > table_cfi(1)%fmima(2)) then
      local_mf = inval(3)
      inval(3) = table_cfi(1)%fmima(2)
    end if

    nvalu(table_cfi(1)%ndcfi+1) = 1_ip

    if (table_cfi(1)%nclas == 4_ip) inval(3) = inval(3) - table_cfi(1)%fmima(1)
    do iclas = 1,table_cfi(1)%ndcfi
      !
      ! clip values to within allowable range 
      !
      acval(iclas) = min(max(inval(iclas),0.0_rp),table_cfi(1)%ivcfi(iclas,table_cfi(1)%nvcfi(iclas)))
      !
      ! get edge points and rib lengths of hypercube
      !
      do iedge = 1,table_cfi(1)%nvcfi(iclas)-1
        if ((table_cfi(1)%ivcfi(iclas,iedge+1) >= acval(iclas)) .and. &
        (table_cfi(1)%ivcfi(iclas,iedge) <= acval(iclas))) then
          ivalu(iclas) = iedge
          dxval(iclas) = (acval(iclas) - table_cfi(1)%ivcfi(iclas,iedge)) &
          / (table_cfi(1)%ivcfi(iclas,iedge+1) - table_cfi(1)%ivcfi(iclas,iedge))
          exit
        end if
      end do
      do icube = 1,ncube
         rtemp = (-1.0_rp)**(real(icube-1_ip,rp)/(2.0_rp**real(iclas-1,rp))+1.0_rp)
        if (rtemp < 0.0_rp) then
          vcube(iclas,icube) = ivalu(iclas)
        else if (rtemp > 0.0_rp) then
          vcube(iclas,icube) = ivalu(iclas) + 1_ip
        end if
      end do
    end do
    !
    ! order distance in magnitude from high to low
    !
    do iclas = 1,table_cfi(1)%ndcfi
      order(iclas) = iclas
    end do
    iclas = 1_ip
    do while (iclas < table_cfi(1)%ndcfi)
      rtemp = 0.0_rp
      do jclas = iclas,table_cfi(1)%ndcfi
        if (dxval(order(jclas)) >= rtemp) then
          rtemp = dxval(order(jclas))
          itemp = jclas
        end if
      end do
      idumi = order(iclas)
      order(iclas) = order(itemp)
      order(itemp) = idumi
      iclas = iclas + 1_ip
    end do
    !
    ! Get edge points of line elements DX over which is interpolated
    !
    do iclas = 1,table_cfi(1)%ndcfi
      dimdx(1,iclas) = 0_ip
    end do
    do iclas = 2,table_cfi(1)%ndcfi+1
      do jclas = 1,table_cfi(1)%ndcfi
        dimdx(iclas,jclas) = dimdx(iclas-1,jclas)
        if (order(iclas-1) == jclas) dimdx(iclas,jclas) = 1_ip
      end do
    end do
    !
    ! Get point numbers
    !
    numpo(1) = 1_ip
    do iclas = 1,table_cfi(1)%ndcfi
      numpo(iclas+1) = 1_ip
      do jclas = 1,table_cfi(1)%ndcfi
        numpo(iclas+1) = numpo(iclas+1) &
                       + int(real(dimdx(iclas+1,jclas),rp)*2.0_rp**real(jclas-1,rp))
      end do
    end do
    !
    ! Obtain line numbers of edge points
    ! 
    do iclas = 1,table_cfi(1)%ndcfi+1
      nline(numpo(iclas)) = 1_ip
      do jclas = 1,table_cfi(1)%ndcfi
        nmult = 1_ip
        do iedge = jclas,table_cfi(1)%ndcfi
          nmult = nmult*nvalu(iedge+1)
        end do
        nline(numpo(iclas)) = nline(numpo(iclas)) &
                            + (vcube(jclas,numpo(iclas))-1_ip)*nmult
      end do
    end do
    !
    ! Get Function values of corresponding edge points
    !
    do iclas = 1,table_cfi(1)%ndcfi+1
      do jclas = 1,table_cfi(1)%nccfi
        fuval(iclas,jclas) = table_cfi(1)%table(nline(numpo(iclas)),jclas)
      end do
    end do
    !
    ! Get interpolated value
    !
    do jclas = 1,table_cfi(1)%nccfi
      reval(jclas) = fuval(1,jclas)
    end do
    do iclas = 1,table_cfi(1)%ndcfi
      do jclas = 1,table_cfi(1)%nccfi
        reval(jclas) = reval(jclas) + (fuval(iclas+1,jclas) - fuval(iclas,jclas))*dxval(order(iclas))
      end do
    end do

    reval(1) = 0.0_rp
    reval(17) = 0.0_rp

    if (local_mf < table_cfi(1)%fmima(1)) then
      do iclas = 1,table_cfi(1)%nfcfi
        reval(iclas+1) = table_cfi(1)%inval(1,iclas) * (1.0_rp-(local_mf / table_cfi(1)%fmima(1))) & 
                       + reval(iclas+1) * (local_mf / table_cfi(1)%fmima(1))
      end do
    else if (local_mf > table_cfi(1)%fmima(2))  then
      do iclas = 1,table_cfi(1)%nfcfi
        reval(iclas+1) = reval(iclas+1) * (1.0_rp-((local_mf-table_cfi(1)%fmima(2))/(1.0_rp-table_cfi(1)%fmima(2)))) &
                       + table_cfi(1)%inval(2,iclas) * ((local_mf - table_cfi(1)%fmima(2))/(1.0_rp - table_cfi(1)%fmima(2)))
      end do
    end if

  end select
  
end subroutine tem_getval
