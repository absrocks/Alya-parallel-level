subroutine nsa_reaous
!-----------------------------------------------------------------------
!****f* Nastal/nsa_reaous
! NAME 
!    nsa_reaous
! DESCRIPTION
!    This routine reads the output strategy 
! USES
!    ecoute
! USED BY
!    nsa_turnon
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk
  use      mod_iofile

  use      def_nastal
use mod_lodi_nsa, only: CHARAs
    use mod_ecoute, only :  ecoute
    
  implicit none
  
  integer(ip)    :: iwopo,ivars,dummi, coso 
  character(35)  :: messa
  


  if(kfl_paral<=0) then
     !
     ! Initializations.
     !

     kfl_chkpo_nsa(1)  =  0        ! Don't read a restart file
     kfl_chkpo_nsa(2)  =  1        ! Write a read a restart file
     vtkrestart_nsa    =  0        ! Output restart interval for the VTK Kessler output
     kfl_pro2d_nsa     =  0        ! 2D profile distributions
     npopr_nsa         =  0
 
     bodyr_nsa =     0.0_rp
     kfl_bodyf_nsa =   0

     !
     ! Reach the section.
     !
     rewind(lisda)
     call ecoute('nsa_reaous')
     do while(words(1)/='OUTPU')
        call ecoute('nsa_reaous')
     end do
     !
     ! Begin to read data.
     !
     !-----------------------------------------------------------------------
     ! ADOC[0]> $-----------------------------------
     ! ADOC[0]> $-- Output and postprocess
     ! ADOC[0]> $-----------------------------------
     ! ADOC[0]> OUTPUT_POSTPROCESS
     !-----------------------------------------------------------------------
     do while(words(1)/='ENDOU')
        call ecoute('nsa_reaous')

        call posdef(2_ip,dummi)
        
        if(words(1)=='VTKRE' .or. words(1)=='VTKOU') then
           vtkrestart_nsa = param(1)
        end if

        if(words(1)=='POSTP') then
          if(words(2)=='CHARA') then
            if(words(3)=='STEPS') then
CHARAs%idofn = getrea('STEPS', 1.0_rp, 'LODI, sigma parameter')
            endif 
          end if
        end if


     end do

     !-----------------------------------------------------------------------
     ! ADOC[0]> END_OUTPUT_POSTPROCESS
     !-----------------------------------------------------------------------

  end if


50 format(2a)
  
end subroutine nsa_reaous
    
