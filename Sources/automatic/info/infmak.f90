subroutine infmak() 
  !------------------------------------------------------------------------
  !****f* info/infmak 
  ! NAME  
  !    infmak 
  ! DESCRIPTION 
  !    Output makefile info 
  ! USES 
  !    outfor
  ! USED BY 
  !*** 
  !------------------------------------------------------------------------ 
  use def_master 
  implicit none 
  
  coutp(1)  = 'COMPILED BY:      am2455'   
  coutp(2)  = 'DATE:             Thu Sep  7 23:49:31 EDT 2023'   
  coutp(3)  = 'HOSTNAME:         n0004'   
  coutp(4)  = 'ARCHITECTURE:     x86_64'   
  coutp(5)  = 'OPERATING SYSTEM: Linux 4.18.0-372.26.1.el8_6.x86_64'   
  coutp(6)  = 'END'                           
  coutp(10) = 'F90:             mpif90'      
  coutp(11) = 'FPPFLAGS:        -fpp'    
  coutp(12) = 'FCFLAGS:         -module $O -c' 
  coutp(13) = 'FOPT:            -O3' 
  coutp(14) = 'CSALYA:          -DNDIMEPAR -DVECTOR_SIZE=16 -DMETIS' 
  coutp(15) =  'END'                          
  
end subroutine infmak 
