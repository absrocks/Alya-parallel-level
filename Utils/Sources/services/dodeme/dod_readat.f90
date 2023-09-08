!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_readat.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Read data
!> @details Read input data
!> @} 
!-----------------------------------------------------------------------
subroutine dod_readat
  use def_parame
  use def_master
  use def_dodeme 
  use def_inpout
  use def_domain
  use mod_ecoute, only :  ecoute
  use mod_messages, only : livinf
  implicit none
  integer(ip)           :: isubd,jsubd,intij,iboun
  integer(ip)           :: pelty,kelem,ielem,pnodb
  integer(ip)           :: knodb(mnodb)

  if( INOTSLAVE ) then
     !
     ! Initialization
     !
     nsubd = 1                            ! Number of subdomains

     !-------------------------------------------------------------------
     !
     ! Numerical data
     !
     !-------------------------------------------------------------------

     call ecoute('dod_readat')
     do while(words(1)/='NUMER')
        call ecoute('dod_readat')
     end do
     do while(words(1)/='ENDNU')
        call ecoute('dod_reanut')
        if( words(1) == 'SUBDO' ) then
           nsubd = getint('SUBDO',1_ip,'#Number of subdomains')
        end if
     end do

     !-------------------------------------------------------------------
     !
     ! Subdomain data
     !
     !-------------------------------------------------------------------

     call dod_memall(1_ip)

     call ecoute('dod_readat')
     do while( words(1) /= 'SUBDO' )
        call ecoute('dod_readat')
     end do

     do while( words(1) /= 'ENDSU' )

        if( words(1) == 'ELEME' ) then
           call livinf(0_ip,'READ ELEMENT SUBDOMAINS',0_ip)
           call ecoute('dod_readat') 
           kelem = 0
           do while(words(1)/='ENDEL')
              kelem = kelem+1
              ielem = int(param(1))
              isubd = int(param(2))
              lsubd_nelem(ielem) = isubd
              pelty = ltype(ielem)
              call ecoute('dod_readat') 
           end do
           if( kelem /= nelem ) call runend('READAT: WRONG LIST OF ELEMENT SUBDOMAIN')
        end if
        call ecoute('dod_readat')
     end do

     !-------------------------------------------------------------------
     !
     ! Interface definition
     !
     !-------------------------------------------------------------------

     call ecoute('dod_readat')
     do while( words(1) /= 'INTER' )
        call ecoute('dod_readat')
     end do
     do while( words(1) /= 'ENDIN' )

        if( words(1) == 'CONNE' ) then
           !
           ! Boundary connectivity: LNSUB_NBOUN(IBOUN)
           !
           if( words(2) == 'UNKNO' ) then

              call ecoute('dod_readat')
              do while( words(1) /= 'ENDCO')
                 pnodb = int(param(2))
                 knodb(1:pnodb) = int(param(3:2+pnodb))
                 call finbou(pnodb,knodb,iboun)
                 if( iboun == 0 ) then
                    call runend('DOD_READAT: COULD NOT FIND BOUNDARY')
                 else
                    prescribed_boundaries(iboun) = int(param(3+pnodb))
                 end if
                 call ecoute('reabcs')
              end do

           else

              call ecoute('dod_readat')
              do while( words(1) /= 'ENDCO' )
                 iboun = int(param(1))
                 if( iboun < 1 .or. iboun > nboun ) call runend('DOD_READAT: WRONG BOUNDARY CONNECTIVITY')
                 prescribed_boundaries(iboun) = int(param(2),ip)
                 call ecoute('dod_readat')
              end do

           end if

        else if( words(1) == 'BOUND' ) then
           ! 
           ! Interface types: INTYP(ISUBD,JSUBD)
           !
           do isubd = 1,nsubd 
              call ecoute('dod_readat')         
              do jsubd = 1,nsubd
                 intyp_dod(isubd,jsubd) = int(param(jsubd))
              end do
           end do
           call ecoute('dod_readat')         
           if( words(1) /= 'ENDBO' ) call runend('READAT: WRONG BOUNDARY-TYPE CARD')

        else if( words(1) == 'DEFIN' ) then
           !
           ! Interface definition: ICTOP_DID(INTIJ)
           !
           intij = getint('DEFIN',1_ip,'#INTERFACE DEFINITION')
           do while( words(1) /= 'ENDDE' )
              if( words(1) == 'TYPET' ) then
                 if(      words(2) == 'CHIME' ) then
                    ictop_dod(intij) = DOD_CHIMERA
                 else if( words(2) == 'PATCH' ) then
                    ictop_dod(intij) = DOD_PATCH
                 else if( words(2) == 'HOLED' ) then
                    ictop_dod(intij) = DOD_HOLED_PATCH
                 else if( words(2) == 'PRESC' ) then
                    ictop_dod(intij) = DOD_PRESCRIBED
                 end if
              end if
              call ecoute('dod_readat')
           end do
        end if
        call ecoute('dod_readat')

     end do

  end if
end subroutine dod_readat
