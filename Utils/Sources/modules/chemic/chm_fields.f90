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

   do ipoin=1,npoin
      do iclas = 1,nspec_chm
         if( kfl_fixno_chm(iclas,ipoin) == 0 ) then
            bvess_chm(iclas,ipoin) = xfiel(iclas+kfl_field_chm(1)-1_ip)%a(1,ipoin,1)
         end if
      enddo
   enddo

end subroutine chm_fields
 
