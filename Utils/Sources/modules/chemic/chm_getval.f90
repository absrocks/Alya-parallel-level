subroutine chm_getval(itask,inval_global,reval,yscale,cvar) 
  !-----------------------------------------------------------------------
  !**** 
  ! NAME 
  !    chm_getval
  ! DESCRIPTION
  !    Returns properties from table for a given combination of reaction 
  !    progress variable, mixture fraction, their variances (stored
  !    in conce) and the scale reaction progress for the unscale model
  ! USES
  ! USED BY
  !   
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      : ip,rp
  use def_master, only      : wmean,conce,table_cfi
  use def_chemic, only      : kfl_uncfi_chm
  implicit none
  real(rp),     intent(in)  :: inval_global(table_cfi(1)%ndcfi)
  real(rp),     intent(out) :: reval(table_cfi(1)%nccfi)
  real(rp),     intent(out) :: yscale
  real(rp),     intent(out) :: cvar

  integer(ip)               :: iclas,iedge,icube,jclas,ncube,itask,ii,loc
  integer(ip)               :: ivalu(table_cfi(1)%ndcfi),vcube(table_cfi(1)%ndcfi,16),numpo(table_cfi(1)%ndcfi+1)
  integer(ip)               :: order(table_cfi(1)%ndcfi),itemp,idumi,dimdx(table_cfi(1)%ndcfi+1,table_cfi(1)%ndcfi)
  integer(ip)               :: nline(16),nmult,nvalu(table_cfi(1)%ndcfi+1)
  real(rp)                  :: local_mf
  real(rp)                  :: y_c_eq,delta_y_c,delta_f
  real(rp)                  :: dxval(table_cfi(1)%ndcfi),acval(table_cfi(1)%ndcfi),rtemp
  real(rp)                  :: fuval(16,table_cfi(1)%nccfi),fuval_temp(16,table_cfi(1)%nccfi)
  real(rp)                  :: inval(table_cfi(1)%ndcfi)
  real(rp)                  :: aux_value(table_cfi(1)%nvcfi(3))
  real(rp)                  :: ry_c_eq2

  select case (itask)

  !
  ! premixed or mean of mixture fraction is inside the flamability limit
  !
  case(1_ip)

    ncube = int(2.0_rp**real(table_cfi(1)%ndcfi,rp),ip)
  
    do iclas = 1, table_cfi(1)%ndcfi
      nvalu(iclas) = table_cfi(1)%nvcfi(iclas)   ! N_yc, N_yc*yc, N_f, N_fv
      inval(iclas) = inval_global(iclas)         ! Y_c, Y_c*Y_c, f, f_var
    end do
    !
    ! For unscale model: c = Y_c / Y_c_eq
    ! 
    if ( kfl_uncfi_chm == 1_ip ) then
       loc = minloc(abs(table_cfi(1) % ymass(:,1) - inval (3)), 1)

       if (inval(3) < table_cfi(1) % ymass(loc,1)) loc = loc - 1_ip
       !
       ! Linear interpolation: Y_c = Y_c_1 + (f - f_1) * Delta_Y_c / Delta_f 
       !
       delta_f   = table_cfi(1) % ymass(loc + 1,1) - table_cfi(1) % ymass(loc,1)
       delta_y_c = table_cfi(1) % ymass(loc + 1,2) - table_cfi(1) % ymass(loc,2)

       if (delta_f /= 0.0_rp ) then
          y_c_eq = table_cfi(1) % ymass(loc,2) + (inval(3) - table_cfi(1) % ymass(loc,1)) * delta_y_c / delta_f
       else
          y_c_eq = table_cfi(1) % ymass(loc,2)
       endif
       !
       !  c = Y_c / Y_c_eq
       !
       yscale   = inval(1) / y_c_eq 
       inval(1) = yscale
       !
       ! c_v = Y_c^2 / Y_c_eq^2 -  (Y_c / Y_c_eq)^2 
       !
       ry_c_eq2 = y_c_eq * y_c_eq
       ry_c_eq2 = 1.0_rp / ry_c_eq2
!!DMM       cvar     = inval(2) * ry_c_eq2 - yscale * yscale
       cvar     = inval(2) * ry_c_eq2
       inval(2) = cvar

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
        rtemp = (-1.0_rp)**(int((icube-1_ip)/(2.0_rp**real(iclas-1)))+1_ip) 
        if (rtemp < 0.0_rp) then
          vcube(iclas,icube) = ivalu(iclas)
        else if (rtemp > 0.0_rp) then
          vcube(iclas,icube) = ivalu(iclas) + 1_ip
        end if
      end do
    end do
    !
    ! get line numbers in table
    !
    do icube = 1,ncube
      nline(icube) = 1
      do jclas = 1,table_cfi(1)%ndcfi
        nmult = 1
        do iedge = jclas,table_cfi(1)%ndcfi
          nmult = nmult*nvalu(iedge+1)
        end do
        nline(icube) = nline(icube) &
                     + (vcube(jclas,icube)-1)*nmult
      end do
    end do
    !
    ! Get Function values of corresponding edge points
    !
    do iclas = 1,ncube
      do jclas = 1,table_cfi(1)%nccfi
        fuval(iclas,jclas) = table_cfi(1)%table(nline(iclas),jclas)
      end do
    end do
    !
    ! interpolate values
    !
    do jclas = 1,table_cfi(1)%ndcfi
      do iclas = 1,2**(table_cfi(1)%ndcfi-jclas)
        icube = iclas * 2
        do iedge = 1,table_cfi(1)%nccfi
          fuval_temp(iclas,iedge) = (1-dxval(jclas))*fuval(icube-1,iedge)+dxval(jclas)*fuval(icube,iedge)
        end do
      end do
      do iclas = 1,2**(table_cfi(1)%ndcfi-jclas)
        do iedge = 1,table_cfi(1)%nccfi
          fuval(iclas,iedge) = fuval_temp(iclas,iedge)
        end do
      end do
    end do
    !
    ! Get return value
    !
    do jclas = 1,table_cfi(1)%nccfi
      reval(jclas) = fuval(1,jclas)
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
      nvalu(iclas) = table_cfi(1)%nvcfi(iclas)
      inval(iclas) = inval_global(iclas)
    end do
    !
    ! Clipping mixture fraction to access database 
    !
    if (inval(3) < table_cfi(1)%fmima(1)) then
      local_mf = inval(3)
      inval(3) = table_cfi(1)%fmima(1)
    else if(inval(3) > table_cfi(1)%fmima(2)) then
      local_mf = inval(3)
      inval(3) = table_cfi(1)%fmima(2)
    end if
    !
    ! For unscale model: c = Y_c / Y_c_eq
    ! 
    if ( kfl_uncfi_chm == 1_ip ) then
       loc = minloc(abs(table_cfi(1) % ymass(:,1) - inval (3)), 1)

       if (inval(3) < table_cfi(1) % ymass(loc,1)) loc = loc - 1_ip
       !
       ! Linear interpolation: Y_c = Y_c_1 + (f - f_1) * Delta_Y_c / Delta_f 
       !
       delta_f   = table_cfi(1) % ymass(loc + 1,1) - table_cfi(1) % ymass(loc,1)
       delta_y_c = table_cfi(1) % ymass(loc + 1,2) - table_cfi(1) % ymass(loc,2)

       if (delta_f /= 0.0_rp ) then
          y_c_eq = table_cfi(1) % ymass(loc,2) + (inval(3) - table_cfi(1) % ymass(loc,1)) * delta_y_c / delta_f
       else
          y_c_eq = table_cfi(1) % ymass(loc,2)
       endif
       !
       !  c = Y_c / Y_c_eq
       !
       inval(1) = inval(1) / y_c_eq 
       yscale   = inval(1)

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
        rtemp = (-1.0_rp)**(int((icube-1_ip)/(2.0_rp**real(iclas-1)))+1_ip)
        if (rtemp < 0.0_rp) then
          vcube(iclas,icube) = ivalu(iclas)
        else if (rtemp > 0.0_rp) then
          vcube(iclas,icube) = ivalu(iclas) + 1_ip
        end if
      end do
    end do
    !
    ! get line numbers in table
    !
    do icube = 1,ncube
      nline(icube) = 1
      do jclas = 1,table_cfi(1)%ndcfi
        nmult = 1
        do iedge = jclas,table_cfi(1)%ndcfi
          nmult = nmult*nvalu(iedge+1)
        end do
        nline(icube) = nline(icube) &
                     + (vcube(jclas,icube)-1)*nmult
      end do
    end do
    !
    ! Get Function values of corresponding edge points
    !
    do iclas = 1,ncube
      do jclas = 1,table_cfi(1)%nccfi
        fuval(iclas,jclas) = table_cfi(1)%table(nline(iclas),jclas)
      end do
    end do
    !
    ! interpolate values
    !
    do jclas = 1,table_cfi(1)%ndcfi
      do iclas = 1,2**(table_cfi(1)%ndcfi-jclas)
        icube = iclas * 2
        do iedge = 1,table_cfi(1)%nccfi
          fuval_temp(iclas,iedge) = (1-dxval(jclas))*fuval(icube-1,iedge)+dxval(jclas)*fuval(icube,iedge)
        end do
      end do
      do iclas = 1,2**(table_cfi(1)%ndcfi-jclas)
        do iedge = 1,table_cfi(1)%nccfi
          fuval(iclas,iedge) = fuval_temp(iclas,iedge)
        end do
      end do
    end do
    !
    ! Get return value
    !
    do jclas = 1,table_cfi(1)%nccfi
      reval(jclas) = fuval(1,jclas)
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
  
end subroutine chm_getval
