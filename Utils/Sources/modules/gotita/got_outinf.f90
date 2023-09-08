subroutine got_outinf
  !-----------------------------------------------------------------------
  !****f* Gotita/got_outinf
  ! NAME 
  !    got_outinf
  ! DESCRIPTION
  !    This routine writes on the incomcdropible Navier-Stokes files
  ! USED BY
  !    got_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_gotita
  implicit none
  character(60) :: equat

  if(kfl_paral<=0) then
     !
     ! Write information in Result file
     !
     if(kfl_rstar/=2) then   
        equat=''
     end if
     write(momod(modul)%lun_outpu,200)&
          densi_got,ddrop_got,muair_got,deair_got,kfact_got

  end if
  !
  ! Formats
  !
110 format(/,10x,a)
200 format(//,&
       & 10x,'Physical parameters:',/,&
       & 10x,'--------------------',/&
       & 10x,'Density water=    ',e12.6,/&
       & 10x,'Droplet diameter= ',e12.6,/&
       & 10x,'Air viscosity=    ',e12.6,/&
       & 10x,'Air density=      ',e12.6,/&
       & 10x,'Inertia (K)=      ',e12.6)
111 format(&
       & 10x,a,10(1x,e12.6))
113 format(&
       & 10x,a)

end subroutine got_outinf
      
