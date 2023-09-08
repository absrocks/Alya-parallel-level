subroutine Endste()
  !-----------------------------------------------------------------------
  !****f* master/Endste
  ! NAME
  !    Endste
  ! DESCRIPTION
  !    This routine closes a time step.
  ! USES
  !    Nastin
  !    Temper
  !    Codire
  !    Alefor
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use mod_bourgogne_pinotnoir
  use def_coupli
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: imodu,icoup

  call bourgogne(1_ip)
  !
  ! Initializations
  !
  kfl_gotim = 0
  !
  ! End a time step for each module
  !
  do iblok = 1,nblok
    call moduls(ITASK_ENDSTE)
  end do
  call Kermod(ITASK_ENDSTE)
  !
  ! Postprocess ppm
  !
  call posppm()
  !
  ! Service adapti
  !
  call Adapti(ITASK_ENDSTE)
  !
  ! Live information
  !
  call livinf(nine,' ',zero)
  !
  ! Check if the time evolution has to be stopped or not
  !  
  if( kfl_gotim /= 0 ) then !Some module is still running, then wake up all modules
     do imodu = 1,mmodu-1
        if( kfl_modul(imodu) /= 0 ) then
           momod(imodu) % kfl_stead = 0
        end if
     end do
  end if

  if( cutim >= timef-epsilon(1.0_rp) )  kfl_gotim = 0
  if( ittim >= mitim )                  kfl_gotim = 0
  !
  ! Close postprocess file if necessary
  !
  call openfi(6_ip)
  !
  ! Write restart files
  !
  call restar(two) ! General data
  !
  ! Output memory evolution
  !
  call output_memory_evolution()
  
end subroutine Endste
