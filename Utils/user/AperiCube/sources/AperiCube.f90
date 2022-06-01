program apericube
  !
  ! AperiCube: cube mesh generator
  ! Mariano Vazquez - Barcelona Supercomputing Center
  !

  use def_cputim

  implicit none
  integer(4)    :: istat
  integer       :: &
       ndime,npoin,nelem,nboun,ndofn,npost,nstep,ntior
  integer       :: &
       npord(3), neord(3)
  real          :: &
       ddref(3), diffu(3) , xmaxi,ymaxi,zmaxi,timax
  real          :: &
       copoi(3)
  real, pointer ::  & 
       coord(:,:),  &    ! nodal coordinates
       unkno(:,:,:)      ! nodal unknown
  integer, pointer ::  & 
       mperm(:,:,:)     ! points
  character(80) :: namda

  real :: time1,time2

  ! parameters:
  
  call cputim(time1)

  xmaxi = 1.0
  ymaxi = 1.0
  zmaxi = 1.0

!  npord(1) = 501
!  npord(2) = 501
!  npord(3) = 501

!  npord(1) = 2
!  npord(2) = 2
!  npord(3) = 2
  ! derived parameters:
!  neord(1)= npord(1)-1
!  neord(2)= npord(2)-1
!  neord(3)= npord(3)-1

  neord(1) = 2 !x
  neord(2) = 2 !y
  neord(3) = 2 !z
  ! derived parameters:
  npord(1)= neord(1)+1
  npord(2)= neord(2)+1
  npord(3)= neord(3)+1

  ndime = 3       ! space dimension (always 3!!!)

  write(6,*) ' ---|'
  write(6,*) ' ---| AperiCube '
  write(6,*) ' ---| -----------------'

  nelem= neord(1)*neord(2)*neord(3)
  npoin= npord(1)*npord(2)*npord(3)

  nboun= 2 * (neord(1)*neord(2) + neord(2)*neord(3) + neord(3)*neord(1))
  
  ddref(1)= xmaxi / real(neord(1))
  ddref(2)= ymaxi / real(neord(2))
  ddref(3)= zmaxi / real(neord(3))

  ! allocatable vectors:
  allocate(coord(ndime,npoin),stat=istat)
  allocate(mperm(npord(1),npord(2),npord(3)),stat=istat)

  write(6,*) ' ---|'
  write(6,*) ' ---|   Problem size:'
  write(6,*) ' ---|        npoin:',npoin
  write(6,*) ' ---|        nelem:',nelem
  write(6,*) ' ---|        nboun:',nboun
  write(6,*) ' ---|'
  write(6,*) ' ---|   Defining coord, elem, bound...'
  call acubeDB(&
       ndime,npoin,nelem,nboun,coord,ddref,xmaxi,ymaxi,zmaxi,npord,neord,mperm)

  call cputim(time2)
  time_total = time2-time1
  
  write(6,*) ' ---|'
  write(6,*) ' ---|   Files written:'
  write(6,*) ' ---|   fort.50  -> dom.dat file'
  write(6,*) ' ---|   fort.100 -> coordinates'
  write(6,*) ' ---|   fort.200 -> elements'
  write(6,*) ' ---|   fort.300 -> boundaries'
  write(6,*) ' ---|   fort.400 -> boundary condition codes'
  write(6,*) ' ---|   fort.500 -> field (generic)'
  write(6,*) ' ---|   Codes: '
  write(6,*) ' ---|       1 -> x=0'
  write(6,*) ' ---|       2 -> y=0'
  write(6,*) ' ---|       3 -> z=0'
  write(6,*) ' ---|       4 -> x=1'
  write(6,*) ' ---|       5 -> y=1'
  write(6,*) ' ---|       6 -> z=1'
  write(6,*) ' ---|'
  write(6,*) ' ---|   VERY IMPORTANT: '
  write(6,*) ' ---|   Always remember to include all the code combinations... they are: '
  write(6,*) ' ---|'
  write(6,*) ' ---|    1   ' 
  write(6,*) ' ---|    2     '
  write(6,*) ' ---|    3     '
  write(6,*) ' ---|    4     '
  write(6,*) ' ---|    5     '
  write(6,*) ' ---|    6     '
  write(6,*) ' ---|    1 & 3   ' 
  write(6,*) ' ---|    2 & 3   ' 
  write(6,*) ' ---|    4 & 3   ' 
  write(6,*) ' ---|    5 & 3    '
  write(6,*) ' ---|    1 & 6    ' 
  write(6,*) ' ---|    2 & 6    '
  write(6,*) ' ---|    4 & 6    '
  write(6,*) ' ---|    5 & 6    '
  write(6,*) ' ---|    1  &  2  &  3  ' 
  write(6,*) ' ---|    1  &  2  &  6  ' 
  write(6,*) ' ---|    1  &  3  &  5  '
  write(6,*) ' ---|    1  &  5  &  6  ' 
  write(6,*) ' ---|    2  &  3  &  4  '
  write(6,*) ' ---|    2  &  4  &  6  ' 
  write(6,*) ' ---|    3  &  4  &  5  '
  write(6,*) ' ---|    4  &  5  &  6  ' 
  write(6,*) ' ---|'
  write(6,*) ' ---|'
  write(6,*) ' ---|   CPU time: '
  write(6,*) ' ---|       Total=    ',time_total
  write(6,*) ' ---|'
  write(6,*) ' ---|   Done.'
  write(6,*) ' ---|'




end program apericube
