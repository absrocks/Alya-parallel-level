subroutine fildef(itask) 
  !-----------------------------------------------------------------------
  !****f* master/fildef
  ! NAME
  !    fildef
  ! DESCRIPTION
  !    This subroutine deals with filters
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master 
  use def_kermod
  use def_inpout
  use mod_memchk
  use mod_ecoute, only :  ecoute
  implicit none  
  integer(ip), intent(in)    :: itask
  integer(ip)                :: jfilt,imodu,ipara,ivarp,inumb,jcofi,modu1
  integer(ip)                :: icofi(mfilt),jmodu
  integer(4)                 :: istat
  character(20)              :: wfilt

  if( itask == -2_ip ) then

     !-------------------------------------------------------------------
     !
     ! Deallocate memory
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then
        do jfilt = 1,nfilt
           call memchk(two,istat,memor_dom,'LFILT','fildef',filte(jfilt) % lfilt)
           deallocate(filte(jfilt) % lfilt,stat=istat)
           if(istat/=0) call memerr(two,'LFILT','fildef',0_ip)
        end do
     end if

  else if( itask == -1_ip ) then

     !-------------------------------------------------------------------
     !
     ! Allocate memory
     !
     !-------------------------------------------------------------------

     !
     ! Take off unused filters
     !
     do jfilt = 1,nfilt
        ipara = 0
        inumb = filte(jfilt) % numbe
        do imodu = 1,mmodu
           if( kfl_modul(imodu) /= 0 ) then
              do ivarp = 1,nvarp
                 if( momod(imodu) % postp(1) % nfilt(ivarp) == inumb ) then
                    ipara = ipara + 1
                 end if
              end do
           end if
        end do
        if( ipara == 0 ) filte(jfilt) % numbe = 0
     end do
     !
     ! Error if an non-existing filter is used
     !
     do imodu = 1,mmodu
        if( kfl_modul(imodu) /= 0 ) then
           do ivarp = 1,nvarp
              inumb = momod(imodu) % postp(1) % nfilt(ivarp)
              if( inumb /= 0 ) then
                 jfilt = 0
                 ipara = 0
                 do while ( jfilt < nfilt ) 
                    jfilt = jfilt + 1
                    if( inumb == filte(jfilt) % numbe ) then
                       do jcofi = 1,ncofi 
                          jmodu = filte(jfilt) % modul(jcofi)
                          if( jmodu > 0 ) then
                             if( kfl_modul(jmodu) == 0 ) then
                                wfilt = intost(filte(jfilt) % numbe)
                                call runend('FILDEF: FILTER '//trim(wfilt)//' USES MODULE '&
                                     //namod(jmodu)//' WHICH IS NOT SOLVED')
                             end if
                          end if
                       end do
                       jfilt = nfilt + 1
                    end if
                 end do
                 if( jfilt /= nfilt + 1 ) then
                    wfilt = intost(inumb)
                    call runend('FILDEF: FILTER '//trim(wfilt)//' USED BY MODULE '&
                         //namod(imodu)//' WAS NOT DEFINED')
                 end if
              end if
           end do
        end if
     end do
     !
     ! Allocate memory
     !
     if( INOTMASTER ) then
        do jfilt = 1,nfilt
           if( filte(jfilt) % numbe /= 0 ) then
              allocate(filte(jfilt) % lfilt(npoin),stat=istat) 
              call memchk(zero,istat,memor_dom,'LFILT','fildef',filte(jfilt) % lfilt)
           end if
        end do
     end if

  else if( itask == 1_ip ) then

     !-------------------------------------------------------------------
     !
     ! Used for Parall service
     !
     !-------------------------------------------------------------------

     do jfilt = 1,mfilt
        call iexcha( filte(jfilt) % numbe )
        do jcofi = 1,ncofi
           call iexcha( filte(jfilt) % modul(jcofi) )
           call iexcha( filte(jfilt) % ifilt(jcofi) )
           do ipara = 1,10
              call rexcha( filte(jfilt) % param(jcofi,ipara) )
           end do
        end do
     end do

  else if( itask == 2_ip ) then

     !-------------------------------------------------------------------
     !
     ! Read filters
     !
     !-------------------------------------------------------------------

     if( words(1) == 'FILTE' ) then
        if( words(2) /= 'OFF  ' ) then
           call ecoute('fildef')

           do jfilt = 1,mfilt
              do jcofi = 1,ncofi
                 filte(jfilt) % modul(jcofi) = -1
              end do
           end do

           nfilt = 0
           do jfilt = 1,mfilt
              icofi(jfilt) = 1
           end do

           do while( words(1) /= 'ENDFI' )

              kfilt = nfilt
              inumb = int(param(1))
              if( inumb > mfilt ) call runend('FILDEF: WRING FILTER NUMBER')

              jfilt = 0
              do while( jfilt < nfilt )
                 jfilt = jfilt + 1
                 if( inumb == filte(jfilt) % numbe ) then
                    icofi(inumb) = icofi(inumb) + 1
                    nfilt = jfilt
                    jfilt = nfilt
                 end if
              end do
              if( icofi(inumb) == 1 )     nfilt = nfilt + 1
              if( icofi(inumb) >  ncofi ) call runend('FILDEF: A FILTER CAN ONLY HAVE TWO COMPONENTS')
              if( nfilt > mfilt-1 )       call runend('FILDEF: WRONG FILTER NUMBER')

              imodu = 0
              modu1 = 0
              do while( imodu < mmodu )
                 imodu = imodu + 1
                 if( words(2) == namod(imodu)(1:5) .and. words(2) /= '' ) then   
                    modu1 = imodu
                    imodu = mmodu + 1
                 end if
              end do
              if( words(2) == 'KERNE' ) then
                 modu1 = mmodu
              else if( imodu /= mmodu + 1 ) then
                 call runend('FILDEF: FILTER MUST BE COMPUTED BY THE KERNEL OR A MODULE')
              end if

              filte(nfilt) % numbe = inumb
              filte(nfilt) % ifilt(icofi(inumb)) = int(param(2)) 
              filte(nfilt) % modul(icofi(inumb)) = modu1
              do ipara = 1,10
                 filte(nfilt) % param(icofi(inumb),ipara) = param(2+ipara)
              end do
              if( icofi(inumb) > 1 ) nfilt = kfilt
              call ecoute('fildef')

           end do
        end if
     end if

  else if( itask == 3_ip ) then

     !-------------------------------------------------------------------
     !
     ! Check if current variable should be filtered and point to filter
     !
     !-------------------------------------------------------------------

     kfl_filte = 0
     if( ivapo /= 0 .and. kfl_outfo > 0 ) then
        if( postp(1) % nfilt(ivapo) /= 0 ) then
           do kfilt = 1,nfilt
              if( postp(1) % nfilt(ivapo) == filte(kfilt) % numbe ) then
                 kfl_filte =  kfilt
                 if( INOTMASTER ) gefil => filte(kfl_filte) % lfilt
              end if
           end do
        end if
     end if

  end if

end subroutine fildef
