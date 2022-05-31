subroutine chm_reatab()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_reatab
  ! NAME 
  !    chm_reatab
  ! DESCRIPTION
  !    Read table properties for flamelet combustion model
  ! USES
  ! USED BY
  !    chm_iniunk: initialize flow field
  !    chm_endite: update table properties to start doiter in temper with updated coefficients 
  !    chm_endste: properties available for all modules at the end of time-step
  !                (there is no update of properties during chemic iteration)   
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only          : ip,rp
  use def_master, only          : wmean,conce,inotmaster,visck,condk, &
                                  sphec,therm,lescl,massk 
  use def_domain, only          : npoin,ltypb,nnode,lnodb,nboun
  use def_chemic, only          : kfl_radia_chm,rspec_chm,kfl_wallc_chm, &
                                  kfl_premix_chm,kfl_varYc_chm,kfl_varZ_chm, &
                                  kfl_ufpv_chm, xZr_chm,xZs_chm,kfl_heat_loss_chm, &
                                  table_fw, nclas_chm 
  use mod_interp_tab,  only     : fw_lookup 

  implicit none
  integer(ip)               :: ipoin,iclas,ivalu
  integer(ip)               :: iboun,pblty,inodb,pnodb
  real(rp)                  :: retva(table_fw % main_table % nvar)      ! Values read in from flamelet table: 
                                                              ! (1) S_c (2) Wmean (3) Lambda (4) Mu (5-10) cp_low (11-16) cp_high (17) S_c*c
  real(rp)                  :: tab_scale_control(table_fw % main_table % ndim)  ! Concentration for each variable at each node
  real(rp)                  :: control(table_fw % main_table % ndim)


  if (INOTMASTER) then
    do ipoin = 1,npoin
      !
      ! Initialization
      !
      massk(ipoin,1:nclas_chm) = 0.0_rp

      control = 0.0_rp
      do iclas = 1, table_fw % main_table % ndim
         select case (table_fw % main_table % coords(iclas) % name)
         case ('CMEAN','C    ')
             control(iclas) = conce(ipoin,1,1)
         case ('CVAR ')
             control(iclas) = conce(ipoin,2,1)
         case ('CHIST')
             control(iclas) = xZr_chm(ipoin) + xZs_chm(ipoin)
         case ('ZMEAN','Z    ')
             control(iclas) = conce(ipoin,3,1)
         case ('ZVAR ')
             control(iclas) = conce(ipoin,4,1)
         case ('IMEAN','I    ')
             control(iclas) = therm(ipoin,1)
         end select
      enddo

      call fw_lookup( control, tab_scale_control, table_fw, retva )

      !
      ! update properties 
      !
      select case(kfl_ufpv_chm )
        case (0)
            massk(ipoin,1)  = retva(1)
            lescl(ipoin)    = retva(17)  ! filter{w_c * c}
        case (1)
            massk(ipoin,1)  = 0.0_rp
            lescl(ipoin)    = 0.0_rp
        case (2)
            massk(ipoin,1)  = retva(1)   ! omega_Yc
            lescl(ipoin)    = 0.0_rp
        case (3)
            massk(ipoin,1)  = retva(17)  ! rho * dYc/dt
            lescl(ipoin)    = 0.0_rp

      endselect

      wmean(ipoin,1) = retva(2)

      do iclas = 1,nclas_chm
        condk(ipoin,iclas) = retva(3)
        visck(ipoin,iclas) = retva(4)
      end do
      do ivalu = 1,6 
        sphec(ipoin,ivalu,1) = retva(4+ivalu)
        sphec(ipoin,ivalu,2) = retva(10+ivalu)
      end do

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
