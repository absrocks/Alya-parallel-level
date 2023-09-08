subroutine exm_reaous
!-----------------------------------------------------------------------
!****f* Exmedi/exm_reaous
! NAME 
!    exm_reaous
! DESCRIPTION
!    This routine reads the output strategy 
! USES
!    listen
! USED BY
!    exm_turnon
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk

  use      def_exmedi
  use mod_ecoute, only :  ecoute

  implicit none
  integer(ip) jdrep,kdrep,iwopo,ivars,idumy
  character(35)  :: messa

  if( INOTSLAVE ) then
     
     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute('exm_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('exm_reaous')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDOU')
        call ecoute('exm_reaous')

        call posdef(2_ip,idumy)
       
        if(words(1)=='PARAM') then
        ! Postprocess general parameters
           call ecoute('exm_reaous')
           do while(words(1)/='ENDPA')
              if(words(1)=='ISOCH') then
                 call runend('EXM_REAOUS: ISOCHRONES THRESHOLD MUST BE WRITTEN IN PROPERTIES SECTION')
!                 if(words(2)=='HIGHE') then
!                    thiso_exm(1)=  1.0_rp         ! isochrones are taken when value becomes higher than threshold
!                    thiso_exm(2)= param(2)        ! isochrones threshold
!                 else if(words(2)=='LOWER') then
!                    thiso_exm(1)= -1.0_rp         ! isochrones are taken when value becomes lower than threshold
!                    thiso_exm(2)= param(2)        ! isochrones threshold
!                 end if                 
              end if
              call ecoute('exm_reaous')
           end do

        
        end if
     end do
  end if

end subroutine exm_reaous
    
