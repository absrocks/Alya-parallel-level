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
  integer(4)            :: BC,depo1,depo2,depo3,depo4
  integer(4)            :: depo5,depo6,depo7,depo8,depo9,numdepo
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
     read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord,type,BC
     if (part_type > max_part_type) then
        max_part_type = part_type
     endif
  enddo
  close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     write(6,*) '--|'
     do i=1,max_part_type
        write(6,*) '--|--------------------------------------'
        write(6,*) '--| for particle type=',i
        numdepo=0
        depo1=0
        depo2=0
        depo3=0
        depo4=0
        depo5=0
        depo6=0
        depo7=0
        depo8=0
        depo9=0
        out=0

        open(unit=10,file=fil_name,status='old')
        read(10,*)
        do ii =1,nline
           read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord,type,BC
           if (part_type == i) then
              numdepo=numdepo+1

              if (BC==2)then
                 depo1=depo1+1
              else if (BC==3)then
                 depo2=depo2+1
              else if (BC==4)then
                 depo3=depo3+1
              else if (BC==5)then
                 depo4=depo4+1
              else if (BC==6)then
                 depo5=depo5+1
              else if (BC==7)then
                 depo6=depo6+1   
              else if (BC==0)then   
                 out=out+1   
              endif
           endif

        end do
        close(10)

        write(6,*) '--|'
        write(6,*) '--| number of deposited particle in the anterior part         = ',depo1
        write(6,*) '--| number of deposited particle in the interior meatus       = ',depo2
        write(6,*) '--| number of deposited particle in the middle meatus         = ',depo3
        write(6,*) '--| number of deposited particle in the olfactory             = ',depo4
        write(6,*) '--| number of deposited particle in the rest of nasal cavity  = ',depo5
        write(6,*) '--| number of deposited particle in the pluggin               = ',depo6
        write(6,*) '--| number of particle going throught the outlet               = ',out
        write(6,*) '--| total number of particle deposited                         = ', &
             out+depo1+depo2+depo3+depo4+depo5+depo6
        write(6,*) '--|'
        write(6,*) '--| nasal cavity Deposition efficiency= ',(real(depo1+depo2+depo3+depo4+depo5)/real(numdepo))*100,'%'
        write(6,*) '--|'
        write(6,*) '--| anterior Deposition efficiency= ',(real(depo1)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| interior meatus Deposition efficiency= ',(real(depo2)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| middle meatus Deposition efficiency= ',(real(depo3)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| olfactory Deposition efficiency= ',(real(depo4)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| rest of nasal cavity Deposition efficiency= ',(real(depo5)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| pluggin Deposition efficiency= ',(real(depo6)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--------------------------------------------------------------------------------'    
        write(6,*) '--|',(real(depo1)/real(numdepo))*100,' ',(real(depo2)/real(numdepo))*100,'  ',&
(real(depo3)/real(numdepo))*100,' ',(real(depo4)/real(numdepo))*100

       

     end do
end  program alya_deposition
