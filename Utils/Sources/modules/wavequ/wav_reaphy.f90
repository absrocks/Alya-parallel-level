subroutine wav_reaphy()
  !------------------------------------------------------------------------
  !****f* Wavequ/wav_reaphy
  ! NAME 
  !    wav_rwaphy
  ! DESCRIPTION
  !    Read physical problem
  ! USES
  ! USED BY
  !    wav_turnon
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_wavequ
  use def_domain
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: imate,ielem,ii

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     nmate_wav     = 1
     sourc_wav     = 0.0_rp                               ! Source term  
     !
     ! Reach the section
     !
     call ecoute('wav_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('wav_reaphy')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDPH')
        call ecoute('wav_reaphy')
        if(words(1)=='PROBL') then
           !
           ! Problem definition data
           !
           call ecoute('wav_reaphy')
           do while(words(1)/='ENDPR')
              if(words(1)=='SOURC') then               ! Source term
                 if(exists('CONST')) then
                    kfl_sourc_wav = -1
                    sourc_wav(1)  = getrea('VALUE',0.0_rp,'#Source term')
                 else if(exists('FUNCT')) then
                    do ii=1,nsour_wav
                       sourc_wav(ii)=param(2+ii)
                    end do
                    kfl_sourc_wav = getint('FUNCT',1_ip,  '#Source term function')
                 end if
              end if
              call ecoute('wav_reaphy')
           end do
        else if(words(1)=='PROPE') then
           !
           ! Properties
           !               
           if(words(2)=='MATER') &
                nmate_wav=getint('MATER',1_ip,'#Temperature number of materials')

           call wav_memphy

           call ecoute('wav_reaphy')
           do while(words(1)/='ENDPR')
              if(words(1)=='DENSI') then                    ! Density (rho)
                 if(nmate_wav==1) then
                    densi_wav(1:10,1)=param(1:10)
                 else
                    call ecoute('wav_reaphy')
                    do while(words(1)/='ENDDE')
                       imate=int(param(1))
                       densi_wav(1:10,imate)=param(2:11)
                       call ecoute('wav_reaphy')
                    end do
                 end if
              else if(words(1)=='KAPPA') then               ! Kappa (k)
                 if(nmate_wav==1) then
                    kappa_wav(1:10,1)=param(1:10)
                 else
                    call ecoute('wav_reaphy')
                    do while(words(1)/='ENDKA')
                       imate=int(param(1))
                       kappa_wav(1:10,imate)=param(2:11)
                       call ecoute('wav_reaphy')
                    end do
                 end if
              else if(words(1)=='LAWDE') then               ! Law Density (rho)
                 if(nmate_wav==1) then
                    lawde_wav(1)=0
                 else
                    call ecoute('wav_reaphy')
                    do while(words(1)/='ENDLA')
                       imate=int(param(1))
                       if(words(1)=='CONST') then
                          lawde_wav(imate)=0
                       end if
                       call ecoute('wav_reaphy')
                    end do
                 end if
              else if(words(1)=='LAWKA') then               ! Law Kappa (k)
                 if(nmate_wav==1) then
                    lawka_wav(1)=0
                 else
                    call ecoute('wav_reaphy')
                    do while(words(1)/='ENDLA')
                       imate=int(param(1))
                       if(words(1)=='CONST') then
                          lawka_wav(imate)=0
                       end if
                       call ecoute('wav_reaphy')
                    end do
                 end if
              else if(words(1)=='MATER') then               ! List of materials
                 call ecoute('wav_reaphy')
                 do while(words(1)/='ENDMA')
                    ielem=int(param(1))
                    lmate_wav(ielem)=int(param(2))
                    call ecoute('wav_reaphy')
                 end do
              end if
              call ecoute('wav_reaphy')
           end do
        end if
     end do
     where(lmate_wav==0) lmate_wav=1
  end if

end subroutine wav_reaphy
