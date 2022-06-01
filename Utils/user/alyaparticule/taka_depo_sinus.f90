program alya_target

  implicit none
  integer(4)            :: one4=1,two4=2
  integer(4)            :: i,ii,jj,kk,j,nline,modu
  real(8)               :: xcoord,ycoord,zcoord,time,r,per
  integer(4)            :: subdom,part,max_part,max_particule,part_id,part_type
  character(150)        :: fil_name,name
  character(4)          :: cpart,num
  logical               :: dir_e
  integer(4)            :: max_part_type,type,nline1,depo,out,depo1,numdepo,BC
  integer(4), pointer   :: family(:)
  character(150)        :: nunam_pos1,filsa,fils1,fils2

  nline=0
  part=0
  max_part_type=1

  call GETARG(one4,name)

  fil_name    = trim(name)//'-deposition.pts.csv'

  if(len(trim(name))==0)then
     write(6,*) &
          '--| Usage: alya-particle [name]'
     write(6,*) '--|'
     write(6,*) '--|'
     write(6,*) '--| Try again !'
     write(6,*) '--|'
     stop
  end if


  write(6,*) '--|'
  write(6,*) '--| Alya-deposition separator '
  write(6,*) '--|'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  open(unit=10,file=fil_name,status='old')
  do 
     read(10,*,end=106)
     nline= nline+1
  end do
106 continue
  close(10)
  nline=nline-1
  write(6,*) '--|'
  write(6,*) '--| nline = ',nline
  write(6,*) '--|'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  filsa = 'particle_sinus.csv'

  inquire(file=filsa, exist=dir_e)
  if ( dir_e ) then
     !write(*,*)"exist file"
  else
     open (12, file=filsa, status="new")
     close(12)
  end if

  open(12, file=filsa, status="old")
  write(12,*)'Current time,Particle number,x,y,z'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  fils1 = 'script1.sh'
  fils2 = 'script2.sh'

  inquire(file=fils1, exist=dir_e)
  if ( dir_e ) then
     !write(*,*)"exist file"
  else
     open (13, file=fils1, status="new")
     close(13)
  end if
  
  inquire(file=fils2, exist=dir_e)
  if ( dir_e ) then
     !write(*,*)"exist file"
  else
     open (14, file=fils2, status="new")
     close(14)
  end if
  

  open(13, file=fils1, status="old")
  open(14, file=fils2, status="old")
  write(13,*)'cp header.res reduced_fensap.pts.res' 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  
  open(unit=10,file=fil_name,status='old')
  read(10,*)
  do ii =1,nline
     read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord,type,BC
     if (BC== 5 .AND. type==-2) write(12,83) time,part_id,xcoord,ycoord,zcoord
     if (BC== 6 .AND. type==-2) write(12,83) time,part_id,xcoord,ycoord,zcoord
     if (BC== 7 .AND. type==-2) write(12,83) time,part_id,xcoord,ycoord,zcoord
     if (BC== 20 .AND. type==-2) write(12,83) time,part_id,xcoord,ycoord,zcoord
     if (BC== 21 .AND. type==-2) write(12,83) time,part_id,xcoord,ycoord,zcoord
     if (BC== 22 .AND. type==-2) write(12,83) time,part_id,xcoord,ycoord,zcoord
     !
     if (BC== 5 .AND. type==-2) write(13,84) part_id
     if (BC== 6 .AND. type==-2) write(13,84) part_id
     if (BC== 7 .AND. type==-2) write(13,84) part_id
     if (BC== 20 .AND. type==-2) write(13,84) part_id
     if (BC== 21 .AND. type==-2) write(13,84) part_id
     if (BC== 22 .AND. type==-2) write(13,84) part_id
     !
     if (BC== 5 .AND. type==-2) write(14,85) part_id
     if (BC== 6 .AND. type==-2) write(14,85) part_id
     if (BC== 7 .AND. type==-2) write(14,85) part_id
     if (BC== 20 .AND. type==-2) write(14,85) part_id
     if (BC== 21 .AND. type==-2) write(14,85) part_id
     if (BC== 22 .AND. type==-2) write(14,85) part_id

     
  enddo
  close(10)
  close(12)
  close(13)
  close(14)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

83 format (e12.6,',',i8,',',3(',',e12.6))
84 format ("awk '{if($2 ==",i8," && $5 >= -0.0067) print $0}' fensap.pts.res >> reduced_fensap.pts.res")
85 format ("bash getLine.sh",i8)  
  
end program alya_target
 
