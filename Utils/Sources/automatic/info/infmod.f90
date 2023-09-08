subroutine infmod() 
  !------------------------------------------------------------------------
  !****f* info/infmod 
  ! NAME  
  !    infmod 
  ! DESCRIPTION 
  !    Output module info 
  ! USES 
  !    outfor
  ! USED BY 
  !*** 
  !------------------------------------------------------------------------ 
  use def_master 
  implicit none 
  
  coutp(min(1,size(coutp)))= 'sh: 1: module: not found' 
  coutp(min(2,size(coutp)))= 'END'      
  
end subroutine infmod 
