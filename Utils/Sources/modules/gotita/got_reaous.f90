subroutine got_reaous()
  !-----------------------------------------------------------------------
  !****f* Gotita/got_reaous
  ! NAME 
  !    got_reaous
  ! DESCRIPTION
  !    This routine reads the output strategy for the incomcdropible NS
  !    equations.
  !    The variable numbers are:
  !    1. Droplet velocity
  !    2. Droplet concentration
  !    3. Streamlines
  ! USES
  !    ecoute
  ! USED BY
  !    got_turnon
  !***
  !-----------------------------------------------------------------------
  use def_inpout
  use def_master
  use def_gotita
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: dummi

  if( INOTSLAVE ) then
     !
     ! Initializations
     !
     npp_bound_got = 0            ! Boundary conditions
     kfl_exacs_got = 0            ! Exact solution
     kfl_resid_got = 0
     pos_cutof_got = 0.0_rp       ! Cut-off for postprocess
     !
     ! Reach the section
     !
     call ecoute('got_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('got_reaous') 
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDOU')
        call ecoute('got_reaous')

        call posdef(2_ip,dummi)

        if(words(1)=='CUTOF') then
           pos_cutof_got=getrea('CUTOF',0.0_rp,'#Postprocess cut-off')

        else if(words(1)=='OUTPU') then
           !
           ! Error
           !
           if(words(2)=='ERROR') then
              !
              ! Exact solution
              !
              kfl_exacs_got=getint('SOLUT',1_ip,'#Exact solution')
              expar_got=param(1:nexap_got)
           end if

        end if
     end do

  end if

end subroutine got_reaous
    
