subroutine hlm_fivout()

  !-----------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_fivout.f90
  ! NAME 
  !    hlm_fivout
  ! DESCRIPTION
  !    This routine outputs values of field vectors at given sites 
  !    in an output file 'field_vectors.out'.
  ! USES
  ! USED BY
  !    hlm_fivecs
  !-----------------------------------------------------------------

  use def_parame
  use def_master
  use def_domain
  use def_helmoz
  use mod_postpr
  use mod_iofile

  implicit none

  integer(ip) :: ii,ierror
  character(len=20) :: filename

  filename = 'field_vectors.out'
  open(unit=1,file=filename,iostat=ierror)
  if (ierror == 0_ip) then
	write(1,*) 'ELECTRIC FIELD INTENSITY'
	write(1,*)
	write(1,*) '  SITE       Re{Ex}          Re{Ey}          Re{Ez}  '  
	write(1,*)
	do ii = 1,nsite_hlm
		write(1,1409) ii,real(elefi_hlm(1,ii)),real(elefi_hlm(2,ii)),real(elefi_hlm(3,ii))
	enddo
	write(1,*)
	write(1,*) '  SITE       Im{Ex}          Im{Ey}          Im{Ez}  '  
	write(1,*)
	do ii = 1,nsite_hlm
		write(1,1409) ii,aimag(elefi_hlm(1,ii)),aimag(elefi_hlm(2,ii)),aimag(elefi_hlm(3,ii))
	enddo
	write(1,*)
	write(1,*)
	write(1,*) 'MAGNETIC FIELD INTENSITY'
	write(1,*)
	write(1,*) '  SITE       Re{Hx}          Re{Hy}          Re{Hz}  '  
	write(1,*)
	do ii = 1,nsite_hlm
		write(1,1409) ii,real(magfi_hlm(1,ii)),real(magfi_hlm(2,ii)),real(magfi_hlm(3,ii))
	enddo
	write(1,*)
	write(1,*) '  SITE       Im{Hx}          Im{Hy}          Im{Hz}  '  
	write(1,*)
	do ii = 1,nsite_hlm
		write(1,1409) ii,aimag(magfi_hlm(1,ii)),aimag(magfi_hlm(2,ii)),aimag(magfi_hlm(3,ii))
	enddo
	1409 format (i5,'    ',3e16.8)
	close (unit=1)
  else
	write(*,*) 'ERROR OPENING FILE: ',filename
  endif

end subroutine hlm_fivout
