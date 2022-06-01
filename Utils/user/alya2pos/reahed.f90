!-----------------------------------------------------------------------
!> @addtogroup Alya2pos
!> @{
!> @file    reahed.f90
!> @author  Mariano Vazquez
!> @brief   Header's reader
!> @date    18/02/2013
!> @details Header's reader
!> @} 
!-----------------------------------------------------------------------
subroutine reahed(ii,npart_par,wwww8,iiiii,rrrrr)
  use def_kintyp, only : ip,rp
  implicit none
  integer(ip),  intent(in)  :: ii
  character(8), intent(out) :: wwww8(10)
  integer(4),   intent(out) :: iiiii(10)
  real(8),      intent(out) :: rrrrr(10)
  integer(ip)               :: npart_par
  integer(4)                :: ihead

  read(ii) ihead    ! Header: 1234
  read(ii) wwww8(1) ! AlyaPost
  read(ii) wwww8(2) ! Version
  read(ii) wwww8(3) ! NAME
  read(ii) wwww8(4) ! SCALA/VECTO        
  read(ii) wwww8(5) ! NELEM/NPOIN/NBOUN  
  read(ii) wwww8(6) ! INTEG/REAL         
  read(ii) wwww8(7) ! 4BYTE/8BYTE        
  read(ii) wwww8(8) ! SEQUE/PARAL        
  read(ii) wwww8(9) ! NOFIL/FILTE
  read(ii) iiiii(1) 
  read(ii) iiiii(2) 
  read(ii) iiiii(3) ! # Subdomains
  
  if( wwww8(2)(1:5) /= 'V0001' ) then
     read(ii) iiiii(4) ! Time step
  end if
  if( wwww8(2)(1:5) /= 'V0001' .and. wwww8(2)(1:5) /= 'V0002' ) then
     read(ii) iiiii(5) ! TAG1
     read(ii) iiiii(6) ! TAG2
  end if
  read(ii) rrrrr(1)    ! Time
  
  if( wwww8(2)(1:5) == 'V0001' ) then
     iiiii(4) = int(rrrrr(1),4)
  end if

  if( int(iiiii(3),ip) /= npart_par ) stop
  
end subroutine reahed

