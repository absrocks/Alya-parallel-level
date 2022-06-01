program alya_deposition

  implicit none
  integer(4)            :: one4=1,two4=2
  integer(4)            :: i,ii,jj,kk,j,nline,modu
  real(8)               :: xcoord,ycoord,zcoord,time,r,per
  real(8)               :: part_id,part_type,exist,set
  integer(4)            :: subdom,part,max_part,max_particule
  character(150)        :: fil_name,name
  character(4)          :: cpart,num
  logical               :: dir_e
  integer(4)            :: max_part_type,type,nline1,depo,out,depo1,numdepo
  integer(4), pointer   :: family(:)
  character(150)        :: nunam_pos1,filsa

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

  open(unit=10,file=fil_name,status='old')
  read(10,*)
  do ii =1,nline
     read (10,*) time,part_id,part_type,exist,set,xcoord,ycoord,zcoord
     if (int(part_type) > max_part_type) then
        max_part_type = int(part_type)
     endif
  enddo
  close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(6,*) '--|'
  do i=1,max_part_type
     
     depo=0
     out=0
     numdepo=0

     open(unit=10,file=fil_name,status='old')
     read(10,*)
     do ii =1,nline
        read (10,*) time,part_id,part_type,exist,set,xcoord,ycoord,zcoord
        if (int(part_type) == i) then
           numdepo=numdepo+1
           if (int(exist)==-2)then
              depo=depo+1
           else
              out=out+1
           endif
        endif


     end do
     close(10)
     write(*,*) (real(depo)/real(numdepo))*100
     write(6,*) '--|'
     write(6,*) '--| type particle = ',i
     write(6,*) '--| number of deposited particle = ',depo
     write(6,*) '--| number of particle going throught the outlet/plugg= ',out
     write(6,*) '--| total number of particles= ',numdepo
     write(6,*) '--|'

  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end program alya_deposition
 
