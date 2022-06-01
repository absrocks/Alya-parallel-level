program rotacio
  
  implicit none
  
  integer, parameter :: ip = selected_int_kind(9)   ! 4-byte integer
  integer, parameter :: rp = kind(1.00)            ! double precision
 
  integer(ip)       ::  i, ipoin,npoin, nelem, nnode, point(2), ndime, istat, aux
  real(rp)          ::  &
       angle, new_z(3), old_z(3), rot_ax(3),vecaux(3),rotvec(3), &
       dotpro,&
       modve, sinang, cosang
  character (88)    ::  linea
  integer(ip), pointer    :: elem(:,:)
  real(rp), pointer    :: coord(:,:)


  open(10,file='coord.dom.geo',status='unknown')
  open(20,file='rot-coord.dom.geo',status='unknown')
 ! open(30,file='',status='unknown')

!!  npoin=       763302
!!  nelem=          3941541
!!  point(1) = 33965  ! Node del pla que volem orientar perpendicular a l'eix z
!!  point(2) = 34581  ! Amb el punt(1) i el punt (2) formarem el nou eix z

  npoin= 1150092  
  nelem= 6023281

  ndime = 3
  nnode = 4
  

  allocate(coord(npoin,ndime), stat=istat)
  allocate(elem(nelem,nnode), stat=istat)
 
  read(10, '(a)') linea
  write(20, '(a)') linea

  do i=1,npoin
!!!     read(10,*) aux, coord(i,1:ndime)                   ! node i coordenades
     read(10,*) ipoin, coord(ipoin,1:ndime)                   ! esto es para cuando es NOTSORTED
  enddo

!!  point(1) = 168673  ! Node del pla que volem orientar perpendicular a l'eix z
!!  point(2) = 890120 ! Amb el punt(1) i el punt (2) formarem el nou eix z
!!  new_z(1:ndime) = coord(point(2),1:ndime) - coord(point(1),1:ndime)
 

!! uso estas coordenadas para definir el nuevo eje en vez de unos nodos
  new_z(1)= -166.13  + 127.554
  new_z(2)= -68.9136 - 110.858
  new_z(3)= -146.659 + 8.4514

  modve = sqrt(new_z(1)*new_z(1) + new_z(2)*new_z(2) + new_z(3)*new_z(3))

  new_z = - new_z/modve
  
  old_z(1) = 0.0_rp
  old_z(2) = 0.0_rp
  old_z(3) = 1.0_rp

  call cropro(new_z,old_z,rot_ax)

  cosang = sqrt(dotpro(new_z,new_z)) * sqrt(dotpro(old_z,old_z))
  cosang = dotpro(new_z,old_z) / cosang

  angle  = acos(angle)
  sinang = sin(angle)

!!!! test
  vecaux= new_z
  call cropro(rot_ax,vecaux,rotvec)  
  rotvec = vecaux * cosang + rotvec * sinang + rot_ax * dotpro(rot_ax,vecaux) * (1.0_rp - cosang)
!!!
  
  do i=1,npoin
     vecaux(1:ndime) = coord(i,1:ndime)  
     call cropro(rot_ax,vecaux,rotvec)  
     ! la formula de olinde rodriguez
     rotvec = vecaux * cosang + rotvec * sinang + rot_ax * dotpro(rot_ax,vecaux) * (1.0_rp - cosang)

     write(20,*) i, rotvec
     
  end do

  read(10, '(a)') linea
  write(20, '(a)') linea

  
!  do i=1,2
!     read(10, '(a)') linea
!     write(20, '(a)') linea
!  end do 
!  do i=1,nelem
!     read(10,*) aux, elem(i,1:nnode)
!     write(20,*) aux, elem(i,1:nnode)
!  enddo
!
!  write(20,*)'END_ELEMENTS'

end program rotacio

subroutine cropro(a,b,c)
  implicit none

  integer, parameter :: ip = selected_int_kind(9)   ! 4-byte integer
  integer, parameter :: rp = kind(1.00)            ! double precision

  real(rp)  :: c(3),a(3),b(3)

  c(1)=   a(2)*b(3) - a(3)*b(2) 
  c(2)= - a(1)*b(3) + a(3)*b(1) 
  c(3)=   a(1)*b(2) - a(2)*b(1) 

end subroutine  cropro

function dotpro(a,b)
  implicit none

  integer, parameter :: ip = selected_int_kind(9)   ! 4-byte integer
  integer, parameter :: rp = kind(1.00)            ! double precision

  real(rp)  :: dotpro,a(3),b(3)

  dotpro= a(1) * b(1) + a(2) * b(2) + a(3) * b(3) 

end function dotpro
