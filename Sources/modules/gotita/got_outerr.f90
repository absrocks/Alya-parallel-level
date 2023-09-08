subroutine got_outerr
  !------------------------------------------------------------------------
  !****f* Gotita/got_outerr
  ! NAME 
  !    got_outerr
  ! DESCRIPTION
  !    This routine checks if there are errros and warnings
  ! USES
  ! USED BY
  !    got_turnon
  !***
  !------------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_gotita
  use mod_outfor, only : outfor
  implicit none
  integer(ip)   :: ierro=0,iwarn=0
  character(20) :: messa
  !
  ! Postprocess
  !
  if(kfl_probl_got==2.and.postp(1) % npp_stepi (1)>0) then
     postp(1) % npp_stepi (1)=0
     iwarn=iwarn+1
     call outfor(2_ip,momod(modul) % lun_outpu,&
          'SOLVING ONLY DROPLET VELOCITY: CANNOT POSTPROCESS '//postp(1) % wopos (1, 1))
  end if
  if(kfl_probl_got==2.and.postp(1) % npp_stepi (14)>0) then
     postp(1) % npp_stepi (14)=0
     iwarn=iwarn+1
     call outfor(2_ip,momod(modul) % lun_outpu,&
          'SOLVING ONLY DROPLET VELOCITY: CANNOT POSTPROCESS '//postp(1) % wopos (1,14))
  end if
  if(kfl_probl_got==2.and.postp(1) % npp_stepi (9)>0) then
     postp(1) % npp_stepi (9)=0
     iwarn=iwarn+1
     call outfor(2_ip,momod(modul) % lun_outpu,&
          'SOLVING ONLY DROPLET VELOCITY: CANNOT POSTPROCESS '//postp(1) % wopos (1, 9))
  end if
  if(kfl_probl_got==2.and.postp(1) % npp_stepi (4)>0) then
     postp(1) % npp_stepi (4)=0
     iwarn=iwarn+1
     call outfor(2_ip,momod(modul) % lun_outpu,&
          'SOLVING ONLY DROPLET VELOCITY: CANNOT POSTPROCESS '//postp(1) % wopos (1, 4))
  end if
  !
  ! Check the transient evolution
  !
  if(kfl_timei/=0) then
     if(kfl_timei_got == 0) then
        iwarn=iwarn+1
        call outfor(2_ip,momod(modul) % lun_outpu,&
             'STEADY NAVIER STOKES EQUATIONS IN A TRANSIENT CALCULATION')
     end if
  end if
  !
  ! Time integration scheme and accuracy
  !
  if(kfl_timei_got==1.and.kfl_tiacc_got>2) then
     ierro=ierro+1
     call outfor(1_ip,momod(modul) % lun_outpu,&
          'WRONG TIME INTEGRATION ORDER USING TRAPEZOIDAL RULE. MUST BE 1 OR 2.')
  end if
  !
  ! Time tracking of the subscales
  !
  if(kfl_timei_got==0.and.kfl_sgsti_got/=0) then
     iwarn=iwarn+1
     kfl_sgsti_got=0
     call outfor(2_ip,momod(modul) % lun_outpu,&
          'CANNOT TRACK THE SUBGRID SCALES IN TIME FOR STATIONARY PROBLEM')
  end if
  !
  ! Solver
  !
  if(solve(1)%nkryd==0.and.solve(1)%kfl_algso==8) then
     ierro=ierro+1
     call outfor(1_ip,momod(modul) % lun_outpu,'KRYLOV DIOMENSION MUST BE >0 WHEN USING GMRES SOLVER')
  end if

  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------

  call errors(3_ip,ierro,iwarn,' ')

end subroutine got_outerr
