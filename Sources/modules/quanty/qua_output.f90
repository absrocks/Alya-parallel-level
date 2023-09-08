subroutine qua_output(itask)
  !------------------------------------------------------------------------
  !****f* Quanty/qua_output
  ! NAME 
  !    qua_output
  ! DESCRIPTION
  ! End of a quanty time step 
  ! 
  ! ITASK = -1 ... Initial solution
  ! ITASK =  0 ... When timemarching is true. There is output or
  !                post-process of results if required.
  ! ITASK =  1 ... When timemarching is false. Output and/or 
  !                post-process of results is forced if they have not 
  !                been written previously. 
  ! USES
  !    qua_outvar 
  !    qua_radpos
  !    qua_outbcs
  !    qua_outwit
  !    qua_exaerr
  !    qua_openfi
  !    suplot
  !    pspltm
  ! USED BY
  !    qua_turnof
  !    qua_timste
  !    qua_iniunk
  !    qua_endite
  !    qua_endste
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_quanty
  use def_solver
  use mod_postpr
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: itime,ivari
  integer(ip), save       :: pos_alrea(nvarp_qua)

  select case(itask)

  case(-1_ip)
     !
     ! Initial solution
     !
     do ivari=1,nvarp_qua
        if(npp_iniso_qua==1.and.npp_stepi_qua(ivari)>0) then
           call qua_outvar(ivari)
        end if
     end do

  case(0_ip)
     !
     ! End of a time step and at a given time step
     !
     if(cutim>=pos_tinit_qua) then
        do ivari=1,nvarp_qua
           pos_alrea(ivari)=0
           !
           ! At a given time step
           !           
           if(ittim>=npp_inits_qua.and.pos_alrea(ivari)==0.and.npp_stepi_qua(ivari)>0) then     
              if(mod(ittim,npp_stepi_qua(ivari))==0) then
                 call qua_outvar(ivari)
                 pos_alrea(ivari)=1
              end if
           end if
           !
           ! At a given time
           !
           if(pos_alrea(ivari)==0) then  
              do itime=1,nvart_qua
                 if(   abs(pos_times_qua(itime,ivari)-cutim)<(0.5_rp*dtime).and.&
                      &    pos_times_qua(itime,ivari)>0.0_rp) then
                    call qua_outvar(ivari)
                 end if
              end do
           end if
        end do
     end if

  case(1_ip)
     !
     ! End of the run
     !
     do ivari=1,nvarp_qua   
        if(npp_stepi_qua(ivari)/=0.and.pos_alrea(ivari)==0) then
           call qua_outvar(ivari)
        end if
     end do

  case(2_ip)
     !
     ! End of an inner iteration
     !
     call qua_outvar(1_ip)

  case(3_ip)
     !
     ! End of the total iteration
     !
     call qua_outvar(2_ip)

  end select

end subroutine qua_output
