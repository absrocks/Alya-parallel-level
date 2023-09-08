subroutine got_reaphy()
  !-----------------------------------------------------------------------
  !****f* Gotita/got_reaphy
  ! NAME 
  !    got_reaphy
  ! DESCRIPTION
  !    This routine reads the physical problem definition for the
  !    incomcdropible NS equations.
  ! USES
  !    ecoute
  !    memchk
  !    runend
  ! USED BY
  !    got_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_gotita
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: idime,ipoin

  if(kfl_paral<=0) then
     !
     ! Initializations (defaults)
     !
     kfl_diffu_got = 0                                    ! Diffusion
     kfl_difun_got = 1                                    ! Diffusion function
     kfl_forme_got = 1                                    ! Conservative/non-conservative
     kfl_probl_got = 1                                    ! Problem to be solved
     kfl_timei_got = 0                                    ! Time integration off
     kfl_timec_got = 0                                    ! Continuity time integration off
     kfl_timem_got = 0                                    ! Momentum time integration off
     kfl_velfu_got = 0                                    ! Velocity function
     ddrop_got     = 1.0_rp
     densi_got     = 1.0_rp
     diffu_got(1)  = 1.0e-6
     diffu_got(2)  = 1.0
     gravi_got     = 0.0_rp                               ! Gravity vector
     grnor_got     = 0.0_rp                               ! Gravity norm
     leinf_got     = 1.0_rp
     muair_got     = 1.0_rp
     veair_got     = 1.0_rp
     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute('got_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('got_reaphy')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDPH')
        call ecoute('got_reaphy')
        if(words(1)=='PROBL') then
           !
           ! Problem definition data
           !          
           call ecoute('got_reaphy')
           do while(words(1)/='ENDPR')

              if(words(1)=='TEMPO') then                    ! Temporal evolution
                 if(words(2)=='ON   ') then
                    kfl_timei_got = 1 
                    kfl_timem_got = 1 
                    kfl_timec_got = 1 
                 else if(words(2)=='OFF  ') then
                    kfl_timei_got = 0
                    kfl_timem_got = 0 
                    kfl_timec_got = 0 
                 else if(exists('MOMEN')) then
                    if(getcha('MOMEN','ON   ','#du/dt')=='ON   ') then
                       kfl_timei_got = 1 
                       kfl_timem_got = 1                        
                    end if
                 else if(exists('CONTI')) then
                    if(getcha('CONTI','ON   ','#du/dt')=='ON   ') then
                       kfl_timei_got = 1 
                       kfl_timec_got = 1                        
                    end if
                 end if

              else if(words(1)=='PROBL') then               ! Problem to be solved
                 if(words(2)=='BOTH ') then
                    kfl_probl_got=1
                 else if(words(2)=='DROPL') then
                    kfl_probl_got=2
                 else if(words(2)=='WATER') then
                    kfl_probl_got=3
                 end if

              else if(words(1)=='DIFFU') then               ! Diffusion on
                 if(exists('ON   ')) kfl_diffu_got = 1
                 if(exists('SMOOT')) kfl_diffu_got = 2

              else if(words(1)=='FORM ') then               ! Conservative/non-conservative
                 if(words(2)=='NONCO') then
                    kfl_forme_got=1
                 else if(words(2)=='CONSE') then
                    kfl_forme_got=2
                 end if

              else if(words(1)=='AIRVE'.and.kfl_probl_got/=3) then  ! Air velocity
                 if(words(2)=='VELOC') then
                    kfl_velfu_got=0
                 else if(words(2)=='GIVEN') then
                    kfl_velfu_got=-1
                    call got_memphy(1_ip)
                    do ipoin=1,npoin
                       call ecoute('got_reaphy')
                       do idime=1,ndime
                          veloc_got(idime,ipoin)=param(1+idime)
                       end do
                    end do
                 else if(words(2)=='FUNCT') then
                    kfl_velfu_got=getint('FUNCT',0_ip,'#Air Velocity function')
                 end if
                 
              else if(words(1)=='DROPL'.and.kfl_probl_got==3) then  ! Droplet velocity velocity
                 call got_memphy(2_ip)
                 do ipoin=1,npoin
                    call ecoute('got_reaphy')
                    do idime=1,ndime
                       vdrop(idime,ipoin,1)=param(1+idime)
                    end do
                 end do
                 
              end if
              call ecoute('got_reaphy')
           end do
        else if(words(1)=='PROPE') then
           !
           ! Properties
           !           
           call ecoute('got_reaphy')
           do while(words(1)/='ENDPR')

              if(words(1)=='DIFFU') then                    ! Diffusion (k)
                 diffu_got(1) = param(1)
                 diffu_got(2) = param(2)
                 diffu_got(3) = param(3)
              else if(words(1)=='FUNCT') then               ! Diffusion function
                 if(words(2)=='EXPLO') then
                    kfl_difun_got=1
                 else if(words(2)=='TANH ') then
                    kfl_difun_got=2
                 else if(words(2)=='HEAVI') then
                    kfl_difun_got=3
                 else if(words(2)=='ELEME') then
                    kfl_difun_got=4
                 else if(words(2)=='CONST') then
                    kfl_difun_got=5
                 end if
              else if(words(1)=='DENSI') then               ! Density (rho)
                 densi_got = param(1)
              else if(words(1)=='AIRDE') then               ! Air density
                 deair_got = param(1)
              else if(words(1)=='AIRVI') then               ! Air viscosity
                 muair_got = param(1)
              else if(words(1)=='VELOC') then               ! Characteristic air velocity
                 veair_got = param(1)
              else if(words(1)=='LENGT') then               ! Characteristic length
                 leinf_got = param(1)
              else if(words(1)=='DROPL') then               ! Droplet diameter
                 ddrop_got = param(1)
              end if
              call ecoute('got_reaphy')
           end do
        end if
     end do

  end if

end subroutine got_reaphy
