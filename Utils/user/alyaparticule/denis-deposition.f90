program alya_deposition

  implicit none
  integer(4)            :: one4=1,two4=2
  integer(4)            :: i,ii,jj,kk,j,nline,modu,nbr
  real(8)               :: xcoord,ycoord,zcoord,time,r,pastime,inter
  integer(4)            :: subdom,part,max_part,max_particule,part_id,part_type
  character(150)        :: fil_name,name
  character(4)          :: cpart,num
  logical               :: dir_e

  nline=0
  part=0
  max_part=0
  nbr=0  

  call GETARG(one4,name)
  
  fil_name    = trim(name)//'-deposition.pts.res'

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
  write(6,*) '--| Denis-deposition converter '
  write(6,*) '--|'
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  open(unit=10,file=fil_name,status='old')
  do 
     read(10,*,end=106)
     nline= nline+1
  end do
106 continue
  close(10)
  nline=nline-5
  write(6,*) '--|'
  write(6,*) '--| nline = ',nline
  write(6,*) '--|'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  inquire(file='depo_time.dat', exist=dir_e)
  if ( dir_e ) then
     !write(*,*)"exist file"
  else
     open (12, file="depo_time.dat", status="new")
     close(12)
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(12, file="depo_time.dat", status="old")
  open(unit=10,file=fil_name,status='old')
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  do ii =1,nline
     read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord
     if (zcoord >= 0.04388) then
        !do i=1,10
        inter=time-pastime
        !enddo
        !write(*,*)inter
        if (inter == 0.0000000000000000)then
           max_part= max_part + 1
        else 
           max_part= max_part + 1
           write(12,'(e12.6, I5)')time,max_part
        endif
     end if
     pastime=time

  end do
  close (10)
  close(12)
  write(6,*) '--|'
  write(6,*) '--| deposition number particule  = ',max_part
  write(6,*) '--|'



84 format (e12.6,',',i5,',',i1,',',e12.6,',',e12.6,',',e12.6)
  
end program alya_deposition
 
