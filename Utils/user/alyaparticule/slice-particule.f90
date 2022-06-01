program slice_particule

  implicit none
  integer(4)            :: one4=1,two4=2
  integer(4)            :: i,ii,jj,kk,j,nline,modu,nbr_part
  real(8)               :: xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,time,r,oldtime,initime
  integer(4)            :: subdom,part,max_part,max_particule,part_id,part_type,islice
  character(150)        :: fil_name,name
  character(4)          :: cpart,num
  character(4)          :: cslice
  integer(4), pointer   :: family(:)
  logical               :: dir_e
  real(8)               :: a(8),b(8),c(8),d(8)
  

  nline=0
  part=0

! 
!slice 1
a(1)=1.72627997E-02; b(1)=0.94291675; c(1)=0.33258078; d(1)=-5.10640517E-02
!
!slice 2
a(2)=0.0; b(2)=1.0; c(2)=0.0; d(2)=-8.25544968E-02
!  
!slice 3
a(3)=9.55353398E-03; b(3)=0.99355137; c(3)=-0.11297952; d(3)=-0.11639848
!
!slice 4
a(4)=0.0; b(4)=0.0; c(4)=-1.0; d(4)=-6.01993315E-03
!
!slice 5
a(5)=0.0; b(5)=0.0; c(5)=-1.0; d(5)=-5.3874608E-02
!
!slice 6
a(6)=0.0; b(6)=0.0; c(6)=-1.0; d(6)=-7.43472204E-02
!
!slice 7
a(7)=0.0; b(7)=0.0; c(7)=-1.0; d(7)=-0.11639445
!
!slice 8
a(8)=0.0; b(8)=0.0; c(8)=-1.0; d(8)=-0.18005240 
!

  call GETARG(one4,name)
  call GETARG(two4,num)
  
  fil_name    = trim(name)//'.pts.res'

  if(len(trim(name))==0)then
     write(6,*) &
          '--| Usage: alya-particule [name] '
     write(6,*) '--|'
     write(6,*) '--|'
     write(6,*) '--| Try again motherfucker !'
     write(6,*) '--|'
     stop
  end if

  write(6,*) '--|'
  write(6,*) '--| Process number of particules crossing 8 slices'
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
  write(6,*) '--| nline =',nline
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
     read (10,*) time,part_id,xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,part_type,subdom,r,r,r,r,r
     max_part=max(part_id,part)
     part=max_part
  end do
  close (10)
  write(6,*) '--|'
  write(6,*) '--| max_part =',max_part
  write(6,*) '--|'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
     read (10,*) time,part_id,xcoord,ycoord,zcoord,xvelo,yvelo,zvelo,part_type,subdom,r,r,r,r,r
     initime=time
  close (10)

  allocate (family(max_part))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do islice=1,8
 
   
   write(cslice,'(i1)')islice
   open  (12, file="slice"//trim(cslice)//".dat", status="new")
   oldtime=initime
   write(6,*) '--|'
   write(6,*) '--| Creation of slice =',islice
   write(6,*) '--|'
   
   do i=1,max_part
      family(i) = 0
   enddo
   
   nbr_part=0
   
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
      if (oldtime /= time)then
         !write(*,*)oldtime,nbr_part
         write(12,84)oldtime,nbr_part
         oldtime=time
      endif
      if (islice < 4) then
         if (a(islice)*xcoord+b(islice)*ycoord+c(islice)*zcoord+d(islice)>0.0.and.family(part_id)==0) then
            family(part_id) = 1
            nbr_part=nbr_part+1
         end if
      else if (islice >= 4) then
         if (a(islice)*xcoord+b(islice)*ycoord+c(islice)*zcoord+d(islice)>0.0.and.family(part_id)==0.and.ycoord>0.04) then
            family(part_id) = 1
            nbr_part=nbr_part+1
         end if
      end if
   end do
   close (10)
   close (12)
   
end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



83 format (e12.6,',',i5,',',e12.6,',',e12.6,',',e12.6,',',e12.6,',',e12.6,',',e12.6,',',i1,',',i4',',i1)
84 format (e12.6,'  ',i5)  
  deallocate (family)

end program slice_particule
 
