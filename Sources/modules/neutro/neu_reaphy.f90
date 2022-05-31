!-----------------------------------------------------------------------
!> @addtogroup NeutroInput
!> @{
!> @file    neu_reaphy.f90
!> @author  Guillaume Houzeaux
!> @date    29/03/2016
!> @brief   Read physical data
!> @details Read physical data
!> @} 
!-----------------------------------------------------------------------
subroutine neu_reaphy()

  use def_parame
  use def_inpout
  use def_master 
  use def_neutro
  use def_domain
  use mod_ADR, only : ADR_read_data
  use mod_ecoute, only :  ecoute
  implicit none

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     ! MODULEX
     !
     num_energies_neu   = 1       ! Number of energy groups
     num_directions_neu = 64      ! Number of directions
     kfl_icosa_neu      = 0
     kfl_snord_neu      = 0
     aniso_neu          = 0.0_rp
     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute('neu_reaphy') 
     do while( words(1) /= 'PHYSI' )
        call ecoute('neu_reaphy')
     end do
     call ecoute('neu_reaphy')
     !--><group>
     !-->    <groupName>PHYSICAL_PROBLEM</groupName>
     !
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> $ Physical properties definition
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> PHYSICAL_PROBLEM
     !
     do while( words(1) /= 'ENDPH' )
        !
        ! Read ADR data
        !
        call ADR_read_data(ADR_neu) 

        if(      words(1) == 'ENERG' ) then 
           !
           ! Energy groups
           !
           num_energies_neu = getint('ENERG',1_ip,'#Number of energy groups')

        else if( words(1) == 'DIREC' ) then 
           !
           ! Direction 
           !
           num_directions_neu = getint('DIREC',1_ip,'#Number of directions')

        else if( words(1) == 'ANISO' ) then                    
           !
           ! Linear anisotropy factor
           !
           aniso_neu = param(1)    
           if( aniso_neu < -1.0_rp ) call runend('neu_reaphy: anisotropy factor cannot be lower than -1.0')     
           if( aniso_neu >  1.0_rp ) call runend('neu_reaphy: anisotropy factor cannot be greater than 1.0')     

        else if( words(1)=='METHO')   then       
           !             
           ! Numerical method
           !
           if( words(2) =='DOM  ') then
              if(      words(3) == 'ICO20')  then 
                 num_directions_neu = 20
                 kfl_icosa_neu = 1  
              else if( words(3) == 'ICO80')  then 
                 num_directions_neu = 80         
                 if( ndime == 2 ) num_directions_neu = 40 
                 kfl_icosa_neu = 1
              else if( words(3) == 'SN8  ')   then 
                 num_directions_neu = 80
                 if( ndime == 2 ) num_directions_neu = 40 
                 kfl_snord_neu = 8
              else if( words(3) == 'SN6  ')   then 
                 num_directions_neu = 48
                 if( ndime == 2 ) num_directions_neu = 24        
                 kfl_snord_neu = 6
              else if( words(3) == 'SN4  ')   then 
                 num_directions_neu = 24
                 if( ndime == 2 ) num_directions_neu = 12        
                 kfl_snord_neu = 4
              else if( words(3) == 'SN10 ')   then 
                 num_directions_neu = 120
                 if( ndime == 2 ) num_directions_neu = 60 
                 kfl_snord_neu = 10
              else if( words(3) == 'SN12 ')   then 
                 num_directions_neu = 168
                 if( ndime == 2 ) num_directions_neu = 84 
                 kfl_snord_neu = 12
              else 
                 call runend('neu_reaphy: error reading number of directions')                  
              end if
              
           end if
           
        end if

        call ecoute('neu_reaphy')
        
     end do
     !--></group>
     !
     ! ADOC[0]> END_PHYSICAL_PROBLEM
     !
     
  end if
  
end subroutine neu_reaphy
