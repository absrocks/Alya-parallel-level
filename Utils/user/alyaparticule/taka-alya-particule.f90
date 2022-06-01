program alya_particule

  implicit none
  integer(4)            :: one4=1,two4=2
  integer(4)            :: i,ii,jj,kk,j,nline,modu,filter
  real(8)               :: xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,time,ax,ay,az,cd,re,r
  real(8)               :: modvel,modacc
  integer(4)            :: subdom,part,max_part,max_particule,part_type,part_id
  character(150)        :: fil_name,name
  character(4)          :: cpart,num
  integer(4), pointer   :: family(:)
  logical               :: dir_e

  nline=0
  part=0

  call GETARG(one4,name)
  call GETARG(two4,num)

  fil_name    = trim(name)//'.pts.res'

  if(len(trim(name))==0)then
     write(6,*) &
          '--| Usage: alya-particule [name] [module]'
     write(6,*) '--|'
     write(6,*) '--|'
     write(6,*) '--| Try again !'
     write(6,*) '--|'
     stop
  end if

  write(6,*) '--|'
  write(6,*) '--| Alya-particule converter '
  write(6,*) '--|'


  read (num,'(I4)') modu

  if(len(trim(num))==0)then
     modu=1
     write(6,*) '--|'
     write(6,*) '--| You are processing all the particules !!!'
     write(6,*) '--|'
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  open(unit=10,file=fil_name,status='old')
  do 
     read(10,*,end=106)
     nline= nline+1
  end do
106 continue
  close(10)
  nline=nline-13
  write(6,*) '--|'
  write(6,*) '--| nline = ',nline
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
  do ii =1,nline
     read (10,*) time,part_id,xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,part_type,subdom,ax,ay,az,cd,re
     max_part=max(part_id,part)
     part=max_part
  end do
  close (10)
  write(6,*) '--|'
  write(6,*) '--| max_part = ',max_part
  write(6,*) '--|'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    inquire(file='paraview_particule.csv', exist=dir_e)
    if ( dir_e ) then
       !write(*,*)"exist file"
    else
       open (12, file="paraview_particule.csv", status="new")
       close(12)
    end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(12, file="paraview_particule.csv", status="old")
  open(unit=10,file=fil_name,status='old')
  !write(12,*)'time,part_id,xcoord,ycoord,zcoord,magvel,xvelo,yvelo,zvelo,part_type,subdom,xacc,yacc,zacc,Cd,Re'
  write(12,*)'time,part_id,xcoord,ycoord,zcoord,modvel,modacc,part_type,subdom,Cd,Stokes nb'
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
  do ii =1,nline
     read (10,*) time,part_id,xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,part_type,subdom,ax,ay,az,cd,re
     modvel=sqrt(xvelo*xvelo+yvelo*yvelo+zvelo*zvelo)
     modacc=sqrt(ax*ax+ay*ay+az*az)
     !if (mod(part_id,modu)==0) write(12,83) time,part_id,xcoord,ycoord,zcoord,&
     !     xvelo,yvelo,zvelo,part_type,subdom,ax,ay,az,cd,re
     if (mod(part_id,modu)==0) write(12,84) time,part_id,xcoord,ycoord,zcoord,modvel,modacc,part_type,subdom,cd,re
  end do
  close (10)
  close (12)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! if you want create family with geometrical criteria
  ! (for example at t=0 2 family one above Z=0.006 and one under)
  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  filter=0

  if (filter ==1 ) then

     allocate (family(max_part))
     do i=1,max_part
        family(i) = 0
     enddo

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
     do ii =1,nline
        read (10,*) time,part_id,xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,part_type,subdom,r,r,r,r,r
        if (time <= 0.100000E-02 .and. zcoord  >= 0.00672873) then
           !write(*,*)" TYPE 1 ",part_id,family(part_id)
           if (family(part_id) == 0) family(part_id) = 1  
        else
           !write(*,*)" TYPE 2 ",part_id,family(part_id)
           if (family(part_id) == 0) family(part_id) = 2
        end if
     end do
     close (10)
     write(6,*) '--|'
     write(6,*) '--| we create 2 family one above Z=0.006 and one under'
     write(6,*) '--|' 


     open(12, file="paraview_particule.csv", status="old")
     open(unit=10,file=fil_name,status='old')
     write(12,*)'time,part_id,xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,part_type,subdom,family'
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
     do ii =1,nline
        read (10,*) time,part_id,xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,part_type,subdom,r,r,r,r,r
        if (mod(part_id,modu)==0) write(12,82) time,part_id,xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,part_type,subdom,family(part_id)
     end do
     close (10)
     close (12)
  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


82 format (e12.6,',',i6,',',e12.6,',',e12.6,',',e12.6,',',e12.6,',',e12.6,',',e12.6,',',i2,',',i6',',i1)
83 format (e12.6,',',i6,6(',',e12.6)',',i2,',',i6,5(',',e12.6))
84 format (e12.6,',',i6,5(',',e12.6)',',i2,',',i6,2(',',e12.6))


  if (filter ==1 ) then
     deallocate (family)
  end if

end program alya_particule
