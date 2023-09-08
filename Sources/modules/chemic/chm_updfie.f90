subroutine chm_updfie()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_updfie
  ! NAME 
  !    chm_updfie
  ! DESCRIPTION
  !    This routine updates:
  !    1. Advection
  !    2. Density
  !    3. Temperature
  ! USES
  !    chm_updunk
  ! USED BY
  !    chemic
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  implicit none
  integer(ip) :: ipoin,idime

  if( INOTMASTER ) then

     !-------------------------------------------------------------------
     !
     ! VELOC_CHM: Velocity
     !
     !-------------------------------------------------------------------
 

     if( kfl_advec_chm == -2 ) then
        !
        ! Use VELOC
        !
        veloc_chm => veloc(:,:,1)

     else if( kfl_advec_chm == -1 ) then
        !
        ! From METEO file
        !
        continue

     else if( kfl_advec_chm > 0 ) then
        !
        ! From function
        !
        call chm_velfun(npoin,coord,veloc_chm)

     end if

     !-------------------------------------------------------------------
     !
     ! DENSI_CHM: Density
     !
     !-------------------------------------------------------------------
 
     if( lawde_chm == -2 ) then
        !
        ! Use DENSI
        !
        densi_chm => densi(:,1)

     else if( lawde_chm == -1 ) then
        !
        ! From METEO file
        !
        continue

     end if

     !-------------------------------------------------------------------
     !
     ! TEMPE_CHM: Temperature
     !
     !-------------------------------------------------------------------
 
     if( lawte_chm == -2 ) then
        !
        ! Use TEMPE
        !
        tempe_chm => tempe(:,1)

     else if( lawte_chm == -1 ) then
        !
        ! From METEO file
        !
        continue

     end if

  end if

end subroutine chm_updfie
