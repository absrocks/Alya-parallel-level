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
  
  coutp(min(1,size(coutp)))= 'Path: .' 
  coutp(min(2,size(coutp)))= 'Working Copy Root Path: /home/juancarlos/svn/Morphology' 
  coutp(min(3,size(coutp)))= 'URL: svn+ssh://bsc21735@dt01.bsc.es/gpfs/projects/bsc21/svnroot/Alya/Trunk' 
  coutp(min(4,size(coutp)))= 'Relative URL: ^/Trunk' 
  coutp(min(5,size(coutp)))= 'Repository Root: svn+ssh://bsc21735@dt01.bsc.es/gpfs/projects/bsc21/svnroot/Alya' 
  coutp(min(6,size(coutp)))= 'Repository UUID: 2cd113d9-5f18-0410-b867-96f721f0e0d9' 
  coutp(min(7,size(coutp)))= 'Revision: 13557' 
  coutp(min(8,size(coutp)))= 'Node Kind: directory' 
  coutp(min(9,size(coutp)))= 'Schedule: normal' 
  coutp(min(10,size(coutp)))= 'Last Changed Author: bsc21759' 
  coutp(min(11,size(coutp)))= 'Last Changed Rev: 13505' 
  coutp(min(12,size(coutp)))= 'Last Changed Date: 2019-08-30 06:15:00 -0400 (Fri, 30 Aug 2019)' 
  coutp(min(13,size(coutp)))= 'END'      
  
end subroutine infsvn 
