subroutine qua_reaous()
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_reaous
  ! NAME 
  !    qua_reaous
  ! DESCRIPTION
  !    This routine reads the output strategy for the Schrodinger equation.
  ! INPUT
  ! OUTPUT
  
  ! USES
  !    ecoute
  ! USED BY
  !    qua_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_quanty
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: ivarp

  if(kfl_paral<=0) then
     !
     ! Initializations
     !
     npp_inits_qua          = 0           ! Postprocess initial step
     npp_iniso_qua          = 0           ! Postprocess initial solution
     npp_stepi_qua          = 0           ! Postprocess step interval
     pos_tinit_qua          = 0.0_rp      ! Postprocess initial time
     pos_times_qua          = 0.0_rp      ! Postprocess times   
     !solve_qua(1)%kfl_cvgso = 0           ! Solver convergence
     !
     ! Reach the section
     !
     call ecoute('qua_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('qua_reaous')
     end do
     !
     ! Begin to read data.
     !
     do while(words(1)/='ENDOU')
        call ecoute('qua_reaous')
        if(words(1)=='POSTP') then
           !
           ! Post-process
           !    
           do ivarp=1,nvarp_qua
              if(words(2)==wopos_qua(1,ivarp).and.words(3)=='STEPS') then
                 npp_stepi_qua(ivarp) = &
                      getint('STEPS',1_ip,'#Postprocess step interval for '//wopos_qua(1,ivarp))

                 if(exists('ATTIM')) pos_times_qua(1:nvart_qua,ivarp)=param(4:nvart_qua+3)
              end if
           end do             

        else if(words(1)=='START') then
           !
           ! Starting time and step of post-process
           !
           if(exists('STEP ')) then
              npp_inits_qua = getint('STEP ',0_ip,'#Initial step to start postprocess')  
           end if
           if(exists('TIME ')) then
              pos_tinit_qua = getrea('TIME ',0.0_rp,'#Initial step to start postprocess')
           end if

        else if(words(1)=='OUTPU') then

           if(words(2)=='SOLVE') then
              !
              ! Solver convergence
              !
              solve_qua(1)%kfl_cvgso=1
           end if

        end if

     end do

  end if

end subroutine qua_reaous
    
