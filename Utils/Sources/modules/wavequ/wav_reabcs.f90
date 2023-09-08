subroutine wav_reabcs
  !-----------------------------------------------------------------------
  !****f* wavequ/wav_reabcs
  ! NAME 
  !    wav_reabcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions
  ! USES
  ! USED BY
  !    wav_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_wavequ
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: ipoin,iknbo,pnodb,iboun
  integer(ip) :: knodb(mnodb)

  if(kfl_paral<=0.and.kfl_ptask/=2) then
     !
     ! Initialization
     !
     kfl_onnod_wav = 0     ! If there are b.c.'s on nodes
     kfl_onbou_wav = 0     ! If there are b.c.'s on boundaries
     kfl_absor_wav = 0     ! If there are absorbing b.c.
     !
     ! Reach the nodal-wise section
     !
     call ecoute('wav_reabcs')
     do while(words(1)/='BOUND')
        call ecoute('wav_reabcs')
     end do
     if(exists('UNKNO')) then
        iknbo=1
     else
        iknbo=0 
     end if
     !
     ! Loop over nodes and or boundaries
     !
     do while(words(1)/='ENDBO')
        call ecoute('wav_reabcs')

        if(words(1)=='CODES'.and.exists('BOUND')) then
           !
           ! User-defined codes on boundaries
           !          
           call wav_membcs(2_ip)
           kfl_onbou_wav =  1
           kfl_fixbo     => kfl_fixbo_wav
           call reacod(2_ip)

        else if(words(1)=='ONNOD') then
           if(.not.exists('OFF  ')) then
              !
              ! There are b.c.'s on nodes
              !
              call wav_membcs(1_ip)
              kfl_onnod_wav = 1
              call ecoute('wav_reabcs')
              do while(words(1)/='ENDON')
                 ipoin                = int(param(1))
                 kfl_fixno_wav(ipoin) = int(param(2))
                 bvess_wav(ipoin)     = param(3)
                 call ecoute('wav_reabcs')
              end do
           end if

        else if(words(1)=='ONBOU') then
           if(.not.exists('OFF  ').and.nboun/=0) then
              !
              ! There are b.c.'s on boundaries
              !
              call wav_membcs(2_ip)
              kfl_onbou_wav = 1
              call ecoute('wav_reabcs')
              if(iknbo==0) then
                 do while(words(1)/='ENDON')
                    iboun                = int(param(1))
                    kfl_fixbo_wav(iboun) = int(param(2))
                    call ecoute('wav_reabcs')
                 end do
              else
                 do while(words(1)/='ENDON')
                    pnodb=int(param(2)) 
                    knodb(1:pnodb)=int(param(3:3+pnodb))
                    call finbou(pnodb,knodb,iboun)
                    kfl_fixbo_wav(iboun) = int(param(3+pnodb))
                    call ecoute('wav_reabcs')
                 end do
              end if
           end if
        end if
     end do
     if(kfl_onbou_wav==1) then
        iboun=0
        do while(iboun<nboun)
           iboun=iboun+1
           if(kfl_fixbo_wav(iboun)==4.or.kfl_fixbo_wav(iboun)==5) then
              iboun=nboun
              kfl_absor_wav=1
           end if
        end do
     end if
  end if

end subroutine wav_reabcs
