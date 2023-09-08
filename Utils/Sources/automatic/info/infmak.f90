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
  
  coutp(1)  = 'COMPILED BY:      juancarlos'   
  coutp(2)  = 'DATE:             Wed 11 Sep 2019 05:01:13 PM EDT'   
  coutp(3)  = 'HOSTNAME:         user-X11DAi-N'   
  coutp(4)  = 'ARCHITECTURE:     x86_64'   
  coutp(5)  = 'OPERATING SYSTEM: Linux 5.0.0-27-generic'   
  coutp(6)  = 'END'                           
  coutp(10) = 'F90:             $(MPICH)/mpif90'      
  coutp(11) = 'FPPFLAGS:        -fpp'    
  coutp(12) = 'FCFLAGS:         -module $O -c' 
  coutp(13) = 'FOPT:            -O3' 
  coutp(14) = 'CSALYA:          -DNDIMEPAR -DVECTOR_SIZE=16 -DMETIS' 
  coutp(15) =  'END'                          
  
end subroutine infmak 
