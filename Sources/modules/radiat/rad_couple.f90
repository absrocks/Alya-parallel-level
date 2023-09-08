subroutine rad_couple()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_couple
  ! NAME 
  !    rad_couple
  ! DESCRIPTION
  !    This routine checks which other relevant modules are ON
  ! USES
  ! USED BY
  !    rad_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_radiat
  use mod_memchk

  implicit none
  integer(4)  :: istat
  integer(ip) :: ipoin,ispec

!!  print *,kfl_paral," says hello"
  !
  ! NASTAL has both TEMPE and DENSI
  !
  if (kfl_coupl(ID_NASTAL,ID_RADIAT)==1) then 
     tempe_rad => tempe(:,1)
  !
  ! NASTIN sometimes has DENSI, and it must be coupled to TEMPER to work
  !
  else if  ( kfl_coupl(ID_TEMPER,ID_RADIAT)==1) then
     tempe_rad => tempe(:,1)  !  
  else if ( kfl_atest_rad==0_ip) then
     call runend ("RADIAT NEEDS NASTAL OR NASTIN AND TEMPER TO PROVIDE TEMPERATURE")
  end if  
  !
  ! Finally we must know if we have species worked out by CHEMIC or not
  !
  if (kfl_coupl(ID_CHEMIC,ID_RADIAT)==1) then
     conce_rad => conce(:,:,1)
  else 
     if (INOTMASTER) then
        allocate(conce_rad(npoin,nspec_rad),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'CONCE_RAD','rad_couple',conce_rad) 
        do ipoin=1,npoin
           conce_rad(ipoin,1) = 1.0_rp
        enddo
        do ispec=2,nspec_rad
           do ipoin=1,npoin
              conce_rad(ipoin,ispec) = 0.0_rp
           enddo
        enddo
     else
        allocate(conce_rad(1,nspec_rad),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'CONCE_RAD','rad_couple',conce_rad) 
        conce_rad(1,1) = 1.0_rp
        do ispec=2,nspec_rad
           conce_rad(1,ispec) = 0.0_rp
        enddo
     endif
  end if

end subroutine rad_couple
      
