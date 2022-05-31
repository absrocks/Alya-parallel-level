subroutine chm_fields
  !------------------------------------------------------------------------
  !****f* chemic/chm_fields
  ! NAME 
  !    chm_fields
  ! DESCRIPTION
  !    This routine:
  !    Reads the initial values for the cocentration
  !    It is important to keep the same order of the species in the fields
  !    as used in ALYA (check the order because it may differ from *chm.dat)    
  !
  ! USES
  ! USED BY
  !    chm_iniunk
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  use mod_postpr
  use mod_memchk
  implicit none
  integer(ip)             :: ipoin,iclas

  !----------------------------------------------------------------------
  !
  ! Write concentration fields
  !
  !----------------------------------------------------------------------

   !
   ! Set all fields to zero 
   !
   do ipoin=1,npoin
      do iclas = 1,nclas_chm
         if( kfl_fixno_chm(iclas,ipoin) == 0 ) then
            bvess_chm(iclas,ipoin) = 0.0_rp
         end if
      enddo
   enddo


   if (kfl_model_chm == 4) then
      if (kfl_rstar == 0 .and. kfl_weigh_in_eq_CMC_chm == 1) then
         do ipoin=1,npoin
            bvess_chm(nvar_CMC_chm+1,ipoin) = xfiel(kfl_field_chm(3))%a(1,ipoin,1)
         end do

         if (kfl_solve_enth_CMC_chm /= 0) then
            do ipoin = 1, npoin
               bvess_chm(nvar_CMC_chm,ipoin) = xfiel(kfl_field_chm(2))%a(1,ipoin,1)
            end do
         end if
      else
         call runend('CHEMIC REAPHY: providing conditional values for each field in CMC not implemented yet')
      end if

   else
      if( kfl_spec_name_chm > 1_ip ) then
         do ipoin=1,npoin
            do iclas = 1,kfl_spec_name_chm
               if( kfl_fixno_chm(iclas,ipoin) == 0 ) then
                  bvess_chm( (Field_ind_chm( iclas ) + 1_ip),ipoin) = xfiel(iclas+kfl_field_chm(1)-1_ip)%a(1,ipoin,1)
               end if
            enddo
         enddo
      else
         do ipoin=1,npoin
            do iclas = 1,nspec_chm
               if( kfl_fixno_chm(iclas,ipoin) == 0 ) then
                  bvess_chm(iclas,ipoin) = xfiel(iclas+kfl_field_chm(1)-1_ip)%a(1,ipoin,1)
               end if
            enddo
         enddo
      end if
   end if

end subroutine chm_fields
 
