subroutine ibm_reanut()
  !-----------------------------------------------------------------------
  !****f* Immbou/ibm_reanut
  ! NAME
  !    ibm_reanut
  ! DESCRIPTION
  !    Read IB numerical problem
  ! OUTPUT
  ! USED BY
  !    ibm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  use def_inpout
  use mod_ecoute, only :  ecoute
  implicit none

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_inter_ibm = 2             ! Interpolation
     kfl_nforc_ibm = 2             ! Numerical force calculation
     kfl_advec     = 0             ! Always mesh velocity
     safet_ibm     = 1.0_rp        ! Safety factor for time step    
     beta_ibm      = 0.25_rp       ! Beta Newmark
     gamma_ibm     = 0.5_rp        ! Gamma Newmark
     !
     ! Reach the section 
     !
     call ecoute('ibm_reanut')
     do while(words(1)/='NUMER')
        call ecoute('ibm_reanut')
     end do
     call ecoute('ibm_reanut') 

     do while( words(1) /= 'ENDNU' ) 

        if(words(1)=='FRING') then
           if( words(2) == 'ELEME' ) then
              kfl_inter_ibm = 1
           elseif( words(2) == 'KRIGI' ) then
              kfl_inter_ibm = 2  
           else if( words(2) == 'ALE  ' ) then
              kfl_inter_ibm = 3  
           end if
        else if(words(1)=='FORCE') then
           if( words(2) == 'ORDIN' ) then
              kfl_nforc_ibm = 1
           elseif( words(2) == 'AVERA' ) then
              kfl_nforc_ibm = 2
           end if
        else if(words(1)=='NEWMA') then
           if(exists('BETA ')) then
              beta_ibm = getrea('BETA ',0.25_rp,'#Beta for Newmark')
           end if
           if(exists('GAMMA')) then
              gamma_ibm = getrea('GAMMA',0.5_rp,'#Gamma for Newmark')
           end if

        end if

        call ecoute('ibm_reanut')

     end do
     
     !if( kfl_aleib_ibm >= 1 ) kfl_advec = 1

  end if

end subroutine ibm_reanut
