program alya_deposition

  implicit none
  integer(4)            :: one4=1,two4=2
  integer(4)            :: i,ii,jj,kk,j,nline,modu
  real(8)               :: xcoord,ycoord,zcoord,time,r,per,stk_inst,stk_eff
  integer(4)            :: subdom,part,max_part,max_particule,part_id,part_type
  character(150)        :: fil_name,name
  character(4)          :: cpart,num
  logical               :: dir_e
  integer(4)            :: max_part_type,type,nline1,depo,out
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
     read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord,stk_inst,stk_eff,type
     if (part_type > max_part_type) then
        max_part_type = part_type
     endif
     if (type==-2)then
        depo=depo+1
     else
        out=out+1
     endif
     
  end do
  close(10)
  write(6,*) '--|'
  write(6,*) '--| different number of type particle = ',max_part_type
  write(6,*) '--| number of deposited particle = ',depo
  write(6,*) '--| number of particle going throught the outlet= ',out
  write(6,*) '--|'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nullify (family)
  allocate (family(max_part_type))
  do i=1,max_part_type
     family(i) = 0
  enddo
  nline1=0
  open(unit=10,file=fil_name,status='old')
  read(10,*)
  do ii =1,nline
     read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord,stk_inst,stk_eff,type
     i=1
     if (type==-2)then
        nline1=nline1+1
        do while ( i <= max_part_type)
           if (part_type == i) family(i) = family(i) + 1
           i=i+1
        end do
     endif
  end do
  close(10)

  do i=1,max_part_type
     per= (real(family(i),8)/real(nline1,8))*100.0
     write(*,*)'type particle',i,',number of particle',family(i),'percentage',per,'%'
  enddo
   
  deallocate (family)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
i=1
do while ( i <= max_part_type)
   if      (i<10)then
      write(nunam_pos1,'(i1)') i
   else  
      write(nunam_pos1,'(i2)') i
   end if
   filsa = 'particle_'//trim(nunam_pos1)//'.csv'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   inquire(file=filsa, exist=dir_e)
   if ( dir_e ) then
      !write(*,*)"exist file"
   else
      open (12, file=filsa, status="new")
      close(12)
   end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   open(12, file=filsa, status="old")
   write(12,*)'Current time,Particle number,Particle type,x,y,z,Stk_inst,Stk_eff'

   open(unit=10,file=fil_name,status='old')
   read(10,*)
   do ii =1,nline
      read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord,stk_inst,stk_eff
      if (part_type == i) write(12,83) time,part_id,part_type,xcoord,ycoord,zcoord,stk_inst,stk_eff
   end do

  close(10)   
  close(12)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
   i=i+1
end do


83 format (e12.6,',',i6,',',i2,5(',',e12.6))
end program alya_deposition
 
