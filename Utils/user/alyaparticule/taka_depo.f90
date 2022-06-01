program alya_deposition

  implicit none
  integer(4)            :: one4=1,two4=2
  integer(4)            :: i,ii,jj,kk,j,nline,modu
  real(8)               :: xcoord,ycoord,zcoord,time,r,per,part_id,part_type,exist,BC
  integer(4)            :: subdom,part,max_part,max_particule
  character(150)        :: fil_name,name
  character(4)          :: cpart,num
  logical               :: dir_e
  integer(4)            :: max_part_type,type,nline1,depo,out
  integer(4)            :: depo1,depo2,depo3,depo4
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  depo=0
  out=0
  
  open(unit=10,file=fil_name,status='old')
  read(10,*)
  do ii =1,nline
     read (10,*) time,part_id,part_type,exist,BC,xcoord,ycoord,zcoord
     if (part_type > max_part_type) then
        max_part_type = part_type
     endif

     if (exist==-2)then
        depo=depo+1
     else if (BC==7) then
        out=out+1
     endif
   
     
  end do
  close(10)
  
  write(6,*) '--|'
  write(6,*) '--| different number of type particle = ',max_part_type
  write(6,*) '--| number of deposited particle in the nasal cavity      = ',depo
  write(6,*) '--| number of particle going throught the outlet or plugin= ',out
  write(6,*) '--| total number of particle deposited                    = ',out+depo
  write(6,*) '--|'
  
  write(6,*) '--| Overall Deposition efficiency= ',(real(depo)/real(out+depo))*100,'%' 
  write(6,*) '--|'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  depo1=0
  depo2=0
  depo3=0
  depo4=0
  out=0
  
  open(unit=10,file=fil_name,status='old')
  read(10,*)
  do ii =1,nline
     read (10,*) time,part_id,part_type,exist,BC,xcoord,ycoord,zcoord
     if (part_type > max_part_type) then
        max_part_type = part_type
     endif
    
     if (BC==3)then
        depo1=depo1+1
     else if (BC==4)then
        depo2=depo2+1
     else if (BC==2.OR.BC==6)then
        depo3=depo3+1
     else if (BC==7.OR.BC==5)then   
        out=out+1   
     endif
   
     
  end do
  close(10)
  
  write(6,*) '--|'
  write(6,*) '--| different number of type particle = ',max_part_type
  write(6,*) '--| number of deposited particle in the common nasal meatus     = ',depo1
  write(6,*) '--| number of deposited particle in the nasopharynx             = ',depo2
  write(6,*) '--| number of deposited particle in the Ethmoid sinuses         = ',depo3
  write(6,*) '--| number of particle going throught the outlet or plugin      = ',out
  write(6,*) '--| total number of particle deposited                          = ',out+depo1+depo2+depo3
  write(6,*) '--|'
  
  write(6,*) '--| common nasal meatus Deposition efficiency= ',(real(depo1)/real(out+depo1+depo2+depo3))*100,'%' 
  write(6,*) '--|'
  write(6,*) '--| nasopharynx Deposition efficiency= ',(real(depo2)/real(out+depo1+depo2+depo3))*100,'%' 
  write(6,*) '--|'
  write(6,*) '--| Ethmoid sinuses Deposition efficiency= ',(real(depo3)/real(out+depo1+depo2+depo3))*100,'%' 
  write(6,*) '--|'

  

end program alya_deposition
 
