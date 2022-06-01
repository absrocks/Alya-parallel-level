subroutine read_johans
!!$  
  use KindType
  use def_master
  use InpOut
  implicit none
  integer(ip)   :: i
  real(rp)      :: rdumm, direction, inidir, z_exp(10)
  
  open(90,FILE='practice_data3.txt',STATUS='old')   !!,ERR=101)
!  open(91,FILE='./Ug.txt',STATUS='old') !!,ERR=101)
!  open(92,FILE='./Vg.txt',STATUS='old') !!,ERR=101)
!  read(90,*)  ! read header
!  read(91,*)  ! read header
!  read(92,*)  ! read header
  do i =1,216
     print *, i
     read(90,*) johans(i,1:56)
!     read(91,*) rdumm, u_geos(i,1)
!     read(92,*) rdumm, u_geos(i,2)
     rad_time(i) = float(i-1)*3600.0_rp
     ! print *, i , rad_time(i), rad_heat(i), u_geos(i,1), u_geos(i,2)
     wall_tempe(i) = johans(i, 30)
     rad_heat(i) = johans(i,49)
     ug_time(i) = rad_time(i)
     inidir =    (270.0_rp - johans(1,28))*pi/180
     direction = (270.0_rp-johans(i,28))*pi/180.0 -inidir
     u_geos(i,1) = johans(i,17)*cos(direction)
     u_geos(i,2) = johans(i,17)*sin(direction)
  end do
  
  close(90)
  open(90,FILE='Ugeos.txt') !,STATUS='new')
  open(91,FILE='qR.txt') !,STATUS='new')
  write(90,'(a)') '# time [h]  Ugx  Ugy'
  write(91,'(a)') '# time [h]  Radheat [W/m2]' 
  do i =1, 216
     write(90,'(3(f15.7,2x))')  float(i-1), u_geos(i,1), u_geos(i,2)
     write(91,'(2(f15.7,2x))')  float(i-1), rad_heat(i)
  end do
  close(90)

  close(91)

  !postprocess exoerimental data
  
  open(90,FILE='Uexp_z.txt')

  z_exp(1:10)=(/ 0.5, 1.5, 5.0, 18.0, 40.0, 59.0, 80.0, 98.0, 120.0, 138.0 /)
  write(90,'(a)') '# z  U(time)'
  do i =1, 10
     write (90,'(218(f15.7,2x))' ) z_exp(i), johans(1:216, 6 + i) 
  end do
  close(90)

  open(90,FILE='Direxp_z.txt')
  write(90,'(a)') '# z  Dir(time)'
  do i =1, 10
     write (90,'(218(f15.7,2x))' ) z_exp(i), johans(1:216, i +17 ) 
  end do
  close(90)

  open(90,FILE='Tempexp_z.txt')
  write(90,'(a)') '# z  Temp(time)'
  do i =1, 10
     write (90,'(218(f15.7,2x))' ) z_exp(i), johans(1:216, i +29 ) 
  end do
  close(90)

  open(90,FILE='TKEexp_z.txt')
  write(90,'(a)') '# z  Temp(time)'
  do i =5, 10
     write (90,'(218(f15.7,2x))' ) z_exp(i), johans(1:216, i +35 )*0.5 
  end do
  close(90)

  open(90,FILE='Heatexp_z.txt')
  write(90,'(a)') '# z  Temp(time)'
  do i =5, 10
     write (90,'(218(f15.7,2x))' ) z_exp(i), johans(1:216, i +46 ) 
  end do
  close(90)
  


  
end subroutine read_johans


subroutine read_radiat
!!$  
  use KindType
  use def_master
  use InpOut
  implicit none
  integer(ip)   :: i
  real(rp)      :: rdumm
  open(90,FILE=rad_file,STATUS='old')   !!,ERR=101)
!  open(91,FILE='./Ug.txt',STATUS='old') !!,ERR=101)
!  open(92,FILE='./Vg.txt',STATUS='old') !!,ERR=101)
  read(90,*)  ! read header
!  read(91,*)  ! read header
!  read(92,*)  ! read header
  do i =1,66
     read(90,*) rad_time(i), rad_heat(i)
!     read(91,*) rdumm, u_geos(i,1)
!     read(92,*) rdumm, u_geos(i,2)
     
     rad_time(i) = rad_time(i)*3600.0
     
    ! print *, i , rad_time(i), rad_heat(i), u_geos(i,1), u_geos(i,2)
     
  end do
  close(90)
!  close(91)
!  close(92)
end subroutine read_radiat

subroutine read_vegeo
!!$  
  use KindType
  use def_master
  use InpOut
  implicit none
  integer(ip)   :: i
  real(rp)      :: rdumm
  open(91,FILE='./Ug.txt',STATUS='old') !!,ERR=101)
  read(91,*)  ! read header
  do i =1,66
     read(91,*) ug_time(i), u_geos(i,1), u_geos(i,2)
     
     ug_time(i) = ug_time(i)*3600.0
     
    ! print *, i , rad_time(i), rad_heat(i), u_geos(i,1), u_geos(i,2)
     
  end do
  
  close(91)
end subroutine read_vegeo

subroutine interp_radiat(qrad)
!!$  
  use KindType
  use def_master
  use InpOut
  implicit none
  real(rp), intent(out) :: qrad
  real(rp) ::   facto
  integer(ip)   :: i
  i=1
  do while (ctime.gt.rad_time(i))
     i = i + 1
  end do
!  print *, 'itera', i
  facto = (ctime- rad_time(i-1))/(rad_time(i)-rad_time(i-1))
  qrad = facto*rad_heat(i) + (1.0-facto)*rad_heat(i-1)
  
end subroutine interp_radiat

subroutine interp_wallte(walltempe)
!!$  
  use KindType
  use def_master
  use InpOut
  implicit none
  real(rp), intent(out) :: walltempe
  real(rp) ::   facto
  integer(ip)   :: i
  i=1
  do while (ctime.gt.rad_time(i))
     i = i + 1
  end do
!  print *, 'itera', i
  facto = (ctime- rad_time(i-1))/(rad_time(i)-rad_time(i-1))
  walltempe = facto*wall_tempe(i) + (1.0-facto)*wall_tempe(i-1)
  
end subroutine interp_wallte



subroutine interp_ugeos(velgeo)
!!$  
  use KindType
  use def_master
  use InpOut
  implicit none
  real(rp), intent(out) :: velgeo(2)
  real(rp) ::   facto
  integer(ip)   :: i
  i=1
  do while (ctime.gt.ug_time(i))
     i = i + 1
  end do
!  print *, 'itera', i
  facto = (ctime- ug_time(i-1))/(ug_time(i)-ug_time(i-1))
  velgeo(1:2) = facto*u_geos(i,1:2) + (1.0-facto)*u_geos(i-1,1:2)
  
end subroutine interp_ugeos
