subroutine Filter(itask)
  !-----------------------------------------------------------------------
  !****f* master/Filter
  ! NAME
  !    Filter
  ! DESCRIPTION
  !    This routine computes the filters to be used in postprocess
  ! USES
  !    modules
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use mod_moduls, only : moduls 
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: imodu,ipoin,icofi,ipara,iffil
  integer(ip)             :: ivari,ivarp,npoin_tmp
  !
  ! Where we are
  ! 1 = initial solution
  ! 2 = end of a time step
  ! 3 = end of the run 
  !
!!$  if( itask == 1 ) then
!!$     ittyp = ITASK_INITIA
!!$  else if( itask == 2 ) then
!!$     ittyp = ITASK_ENDTIM
!!$  else if( itask == 3 ) then
!!$     ittyp = ITASK_ENDRUN
!!$  else if( itask == 4 ) then
!!$     ittyp = ITASK_ENDINN
!!$  end if
  ittyp =itask
  if( IMASTER ) then
     npoin_tmp = npoin
     npoin     = 0
  end if

  do kfilt = 1,nfilt

     if( filte(kfilt) % numbe /= 0 ) then
        !
        ! Check if filter must be used for next output
        !
        iffil = 0
        do imodu = 0,mmodu
           if( kfl_modul(imodu) /= 0 ) then
              postp => momod(imodu) % postp
              ivarp =  0
              do while( ivarp < nvarp )
                 ivarp = ivarp + 1
                 ivari = ivarp
                 call posdef(111_ip,ivari)
                 if( ivari /= 0 ) then
                    if( momod(imodu) % postp(1) % nfilt(ivari) == filte(kfilt) % numbe ) then
                       iffil = 1
                       ivarp = nvarp
                    end if
                 end if
              end do
           end if
        end do

        if( iffil == 1 ) then
           !
           ! Filter is used: initialize filter
           !
           if( INOTMASTER ) then
              do ipoin = 1,npoin           
                 filte(kfilt) % lfilt(ipoin) = 0
              end do
           end if
           do icofi = 1,ncofi
              !
              ! Tell modules to compute their filters
              !
              kfl_filte = filte(kfilt) % ifilt(icofi)

              if( kfl_filte /= 0 ) then
                 do imodu = 0,mmodu

                    if( filte(kfilt) % modul(icofi) == imodu ) then 

                       if( INOTMASTER ) then
                          modul =  imodu
                          gefil => filte(kfilt) % lfilt
                          do ipara = 1,10
                             pafil(ipara) = filte(kfilt) % param(icofi,ipara)
                          end do
                       end if

                       if( imodu == mmodu ) then
                          call Kermod(-ITASK_FILTER)
                       else
                          call moduls( ITASK_FILTER)     ! Filter computed by module MODUL
                       end if

                    end if

                 end do
              end if
           end do
        end if
     end if
  end do

  if( IMASTER ) npoin = npoin_tmp 

end subroutine Filter
