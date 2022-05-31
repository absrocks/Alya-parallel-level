subroutine chm_post_lookup(auxvar)
  !------------------------------------------------------------------------
  ! lookup form postprocessingtable on nodes
  !------------------------------------------------------------------------    
  use def_chemic
  use def_master,     only  : conce, therm
  use def_kintyp,     only  : ip,rp
  use def_domain,     only  : npoin 
  use mod_interp_tab, only  : fw_lookup 
  implicit none

  real(rp),  intent(inout) :: auxvar(npoin,*)

  integer(ip)              :: ipoin,ipostvar,iclas 
  real(rp)                 :: retva(posttable_fw % main_table % nvar)
  real(rp)                 :: control(posttable_fw % main_table % ndim)   
  real(rp)                 :: scale_control(posttable_fw % main_table % ndim)


  do ipostvar=1,posttable_fw % main_table % nvar 
     do ipoin=1,npoin 
        auxvar(ipoin,ipostvar) = 0.0_rp
     enddo
  enddo

  do ipoin=1,npoin
      control = 0.0_rp
      do iclas = 1, posttable_fw % main_table % ndim
         select case (posttable_fw % main_table % coords(iclas) % name)
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

      call fw_lookup( control, scale_control, posttable_fw, retva )

      do ipostvar=1,posttable_fw % main_table % nvar
         auxvar(ipoin,ipostvar) = retva(ipostvar) 
      enddo
  end do

end subroutine chm_post_lookup

