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
  
  coutp(min(1,size(coutp)))= 'Currently Loaded Modules:' 
  coutp(min(2,size(coutp)))= '1) easybuild' 
  coutp(min(3,size(coutp)))= '3) GCCcore/12.2.0' 
  coutp(min(4,size(coutp)))= '5) binutils/2.39' 
  coutp(min(5,size(coutp)))= '' 
  coutp(min(6,size(coutp)))= '' 
  coutp(min(7,size(coutp)))= '' 
  coutp(min(8,size(coutp)))= '' 
  coutp(min(9,size(coutp)))= '' 
  coutp(min(10,size(coutp)))= '' 
  coutp(min(11,size(coutp)))= '7) numactl/2.0.16' 
  coutp(min(12,size(coutp)))= '9) impi/2021.7.1' 
  coutp(min(13,size(coutp)))= '11) imkl-FFTW/2022.2.1' 
  coutp(min(14,size(coutp)))= '2) null' 
  coutp(min(15,size(coutp)))= '' 
  coutp(min(16,size(coutp)))= '' 
  coutp(min(17,size(coutp)))= '' 
  coutp(min(18,size(coutp)))= '4) zlib/1.2.12' 
  coutp(min(19,size(coutp)))= '' 
  coutp(min(20,size(coutp)))= '' 
  coutp(min(21,size(coutp)))= '6) intel-compilers/2022.2.1' 
  coutp(min(22,size(coutp)))= '8) UCX/1.13.1' 
  coutp(min(23,size(coutp)))= '' 
  coutp(min(24,size(coutp)))= '' 
  coutp(min(25,size(coutp)))= '10) imkl/2022.2.1' 
  coutp(min(26,size(coutp)))= '12) intel/2022b' 
  coutp(min(27,size(coutp)))= 'END'      
  
end subroutine infmod 
