subroutine infsvn() 
  !------------------------------------------------------------------------
  !****f* info/infsvn 
  ! NAME  
  !    infsvn 
  ! DESCRIPTION 
  !    Output svn info 
  ! USES 
  !    outfor
  ! USED BY 
  !*** 
  !------------------------------------------------------------------------ 
  use def_master 
  implicit none 
  
  coutp(min(1,size(coutp)))= 'END'      
  
end subroutine infsvn 
