subroutine modser
  !------------------------------------------------------------------------
  !****f* master/outerr
  ! NAME 
  !    outerr
  ! DESCRIPTION
  !    This routine checks the order of modules iterations, and compute
  !    number of modules and services
  ! OUTPUT
  !    LMORD(IORDE,IBLOK) ... Module number to run at IORDE position
  !                           within block IBLOK
  !    NMODU ................ Number of modules
  !    NSERV ................ Number of services
  ! USES
  ! USED BY
  !    Reapro
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_master
  implicit none
  integer(ip) :: imodu,iorde,iserv,jblok
  !
  ! Automatic ordering of the modules: LMORD
  !
  if( lmord(1,1) == 0 ) then
     nblok = 1     
     iorde = 0
     !imodu = mmodu
     !if( kfl_modul(imodu) /= 0 ) then
     !   iorde = iorde + 1
     !   lmord(iorde,1) = imodu        
     !end if
     do imodu = 1,mmodu-1
        if( kfl_modul(imodu) /= 0 ) then
           iorde = iorde + 1
           lmord(iorde,1) = imodu
        end if
     end do
  else
     !
     ! Add KERMOD to block 1 if it was not previously included
     !
     !iorde = 0
     !do iblok = 1,nblok
     !   do imodu = 1,mmodu
     !      if( lmord(imodu,iblok) == mmodu ) iorde = 1
     !   end do
     !end do
     !if( iorde == 0 ) then
     !   nblok = nblok + 1
     !   do jblok = nblok,2,-1
     !      iblok = jblok - 1
     !      do imodu = 1,mmodu
     !         lmord(imodu,jblok) = lmord(imodu,iblok) 
     !      end do
     !   end do
     !   iblok = 1
     !   do imodu = 1,mmodu
     !      lmord(imodu,iblok) = 0
     !   end do
     !   lmord(1,iblok) = mmodu
     !end if
     
  end if
  !
  ! Count the number of modules used: NMODU
  !
  nmodu = 0
  do imodu = 1,mmodu
     if( kfl_modul(imodu) /= 0 ) nmodu = nmodu + 1
  end do
  !
  ! Count the number of services used: NSERV
  ! 
  nserv = 0
  do iserv = 1,mserv
     if( kfl_servi(iserv) /= 0 ) nserv = nserv + 1
  end do

end subroutine modser
