program alya_deposition

  implicit none
  integer(4)            :: one4=1,two4=2
  integer(4)            :: i,ii,jj,kk,j,nline,modu
  real(8)               :: xcoord,ycoord,zcoord,time,r,per
  integer(4)            :: subdom,part,max_part,max_particule,part_id,part_type
  character(150)        :: fil_name,name
  character(4)          :: cpart,num
  logical               :: dir_e
  integer(4)            :: max_part_type,type,nline1,depo,out
  integer(4), pointer   :: family(:)
  character(150)        :: nunam_pos1,filsa
  integer(4)            :: depo1,depo2,depo3,depo4,depo5,depo6,depo7
  real(8)               :: slice1,slice2,slice3,slice4,slice5,slice6



  nline=0
  part=0


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
     read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord,type
     if (part_type > max_part_type) then
        max_part_type = part_type
     endif
  enddo
  close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  write(6,*) '--|'
     do i=1,max_part_type
        write(6,*) '--|--------------------------------------'
        write(6,*) '--| for particle type=',i

  depo=0
  out=0
  
  open(unit=10,file=fil_name,status='old')
  read(10,*)
  do ii =1,nline
     read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord,type
     if (part_type == i) then

      if (type==-2)then
         depo=depo+1
      else
         out=out+1
      endif 
     endif 
   
     
  end do
  close(10)
  
  write(6,*) '--|'
  write(6,*) '--| different number of type particle = ',max_part_type
  write(6,*) '--| number of deposited particle in the nasal cavity      = ',depo
  write(6,*) '--| number of particle going throught the outlet or plugin= ',out
  write(6,*) '--| total number of particle deposited                    = ',out+depo
  write(6,*) '--|'
  
  write(6,*) '--| Deposition efficiency= ',(real(depo)/real(out+depo))*100,'%' 
  write(6,*) '--|'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  depo1=0
  depo2=0
  depo3=0
  depo4=0
  depo5=0
  depo6=0
  depo7=0

  slice1=0.0338
  slice2=0.0197
  slice3=0.0035
  slice4=-0.006
  slice5=-0.020
  slice6=-0.036

  open(unit=10,file=fil_name,status='old')
  read(10,*)
  do ii =1,nline
     read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord,type
     if (part_type == i) then

     if (type==-2 .AND. zcoord >= slice1) then
        depo1=depo1+1
     elseif (type==-2 .AND.zcoord <= slice1 .AND.zcoord >= slice2 )then
        depo2=depo2+1
     elseif (type==-2 .AND.zcoord <= slice2 .AND.zcoord >= slice3 )then
        depo3=depo3+1
     elseif (type==-2 .AND.zcoord <= slice3 .AND.zcoord >= slice4 )then
        depo4=depo4+1
     elseif (type==-2 .AND.zcoord <= slice4 .AND.zcoord >= slice5 )then
        depo5=depo5+1
     elseif (type==-2 .AND.zcoord <= slice5 .AND.zcoord >=slice6 )then
        depo6=depo6+1   
     elseif (type==-2 .AND.zcoord <= slice6) then
        depo7=depo7+1
     endif
    endif           
     
  end do
  close(10)

write(6,*) '--| !!!!!!!! VALID ONLY FOR SUBJECT 3 !!!!!!!!!!!!! '
write(6,*) '--|'
write(6,*) '--| number of deposited particle betwenn the nostrill and slice1  = ',depo1
write(6,*) '--| Deposition efficiency=',(real(depo1)/real(out+depo))*100,'%'
write(6,*) '--|'
write(6,*) '--| number of deposited particle betwenn slice1 and slice2        = ',depo2
write(6,*) '--| Deposition efficiency=',(real(depo2)/real(out+depo))*100,'%' 
write(6,*) '--|'
write(6,*) '--| number of deposited particle betwenn slice2 and slice3        = ',depo3
write(6,*) '--| Deposition efficiency=',(real(depo3)/real(out+depo))*100,'%' 
write(6,*) '--|'
write(6,*) '--| number of deposited particle betwenn slice3 and slice4        = ',depo4
write(6,*) '--| Deposition efficiency=',(real(depo4)/real(out+depo))*100,'%'
write(6,*) '--|'
write(6,*) '--| number of deposited particle betwenn slice4 and slice5        = ',depo5 
write(6,*) '--| Deposition efficiency=',(real(depo5)/real(out+depo))*100,'%'
write(6,*) '--|'
write(6,*) '--| number of deposited particle betwenn slice5 and slice6        = ',depo6 
write(6,*) '--| Deposition efficiency=',(real(depo6)/real(out+depo))*100,'%'
write(6,*) '--|'
write(6,*) '--| number of deposited particle betwenn slice6 and output        = ',depo7
write(6,*) '--| Deposition efficiency=',(real(depo7)/real(out+depo))*100,'%' 
write(6,*) '--|'

enddo

end program alya_deposition
 
