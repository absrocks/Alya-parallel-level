program depo_initial 

  implicit none
  integer(4)            :: one4=1,two4=2
  integer(4)            :: i,ii,jj,kk,j,nline,modu,nlinep
  real(8)               :: xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,time,r
  integer(4)            :: subdom,part,max_part,max_particule,part_id,part_type,max_depo
  character(150)        :: fil_name,name,fil_depo
  character(4)          :: cpart,num
  integer(4), pointer   :: family(:)
  logical               :: dir_e

  nline=0
  nlinep=0
  part=0
  max_depo=0

  call GETARG(one4,name)
  
  fil_name    = trim(name)//'.pts.res'
  fil_depo    = trim(name)//'-deposition.pts.res'

  if(len(trim(name))==0)then
     write(6,*) &
          '--| Usage: alya-particule [name]'
     write(6,*) '--|'
     write(6,*) '--|'
     write(6,*) '--| Try again motherfucker !'
     write(6,*) '--|'
     stop
  end if

  write(6,*) '--|'
  write(6,*) '--| Deposition and initial position '
  write(6,*) '--|'
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  open(unit=10,file=fil_depo,status='old')
  do 
     read(10,*,end=106)
     nline= nline+1
  end do
106 continue
  close(10)
  nline=nline-5
  write(6,*) '--|'
  write(6,*) '--| nline of depo file = ',nline
  write(6,*) '--|'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  inquire(file='paraview_depo_initial.csv', exist=dir_e)
  if ( dir_e ) then
     !write(*,*)"exist file"
  else
     open (12, file="paraview_depo_initial.csv", status="new")
     close(12)
  end if

  allocate (family(nline))
  do i=1,nline
     family(i) = 0
  enddo
  
  max_depo=nline

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  write(6,*) '--|'
  write(6,*) '--| number of particule in depo file = ',max_depo
  write(6,*) '--|'
  write(6,*) '--| we take only the particule deposition above Z=0.04388'
  write(6,*) '--|' 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(12, file="paraview_depo_initial.csv", status="old")
  write(12,*)'time,part_id,part_type,xcoord,ycoord,zcoord,initial_or_final'
 
  open(unit=10,file=fil_depo,status='old')
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  do ii =1,nline
     read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord
     if (zcoord >= 0.04388) then
        write(12,84) time,part_id,part_type,xcoord,ycoord,zcoord,'2'
        family(ii) = part_id
     end if
  end do
  close (10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit=10,file=fil_name,status='old')
  do 
     read(10,*,end=107)
     nlinep= nlinep+1
  end do
107 continue
  close(10)
  nlinep=nlinep-13
  !nlinep=16881450
  write(6,*) '--|'
  write(6,*) '--| nline of particule file == ',nlinep
  write(6,*) '--|'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(unit=10,file=fil_name,status='old')
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  do ii =1,nlinep
     read (10,*) time,part_id,xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,part_type,subdom,r,r,r,r,r
     if (time <= 0.170000E-03 ) then
        do i=1,max_depo
           if (family(i-1)==part_id) write(12,84) time,part_id,part_type,xcoord,ycoord,zcoord,'1'
           !if (family(i-1)==part_id) write(*,*) time,part_id,part_type,xcoord,ycoord,zcoord,'1'
        enddo 
     end if
  end do
  close (10)
  close (12)


84 format (e12.6,',',i5,',',i1,',',e12.6,',',e12.6,',',e12.6,',',A)

  deallocate (family)

end program depo_initial
 
