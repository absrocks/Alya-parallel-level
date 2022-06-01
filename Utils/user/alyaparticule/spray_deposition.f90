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
  write(6,*) '--| Alya-deposition spray separator '
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

  depo1=0
  depo2=0
  depo3=0
  depo4=0
  depo5=0
  out=0

  open(unit=10,file=fil_name,status='old')
  read(10,*)
  do ii =1,nline
     read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord,type,BC
     if (BC==2)then ! anterior part
        depo1=depo1+1
     else if (BC==3)then ! inf meatus
        depo2=depo2+1
     else if (BC==4)then ! mid meatus
        depo3=depo3+1
     else if (BC==5)then !  olfa
        depo4=depo4+1
     else if (BC==6)then ! rest nasal
        depo5=depo5+1  
     else if (BC==0)then   
        out=out+1   
     endif            
  enddo
  close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
        write(6,*) '--|'
        write(6,*) '--| number of deposited particle in anterior part          = ',depo1
        write(6,*) '--| number of deposited particle in inferior meatus        = ',depo2
        write(6,*) '--| number of deposited particle in middle meatus          = ',depo3
        write(6,*) '--| number of deposited particle in olfactory part         = ',depo4
        write(6,*) '--| number of deposited particle in rest of nasal          = ',depo5
        write(6,*) '--| number of particle going throught the outlel, trachea...= ',out
        write(6,*) '--| total number of particle deposited                          = ', &
             out+depo1+depo2+depo3+depo4+depo5
        write(6,*) '--|'
        write(6,*) '--| total spray Deposition efficiency= ',(real(depo1+depo2+depo3+depo4+depo5)/real(nline))*100,'%'
        write(6,*) '--|'
        write(6,*) '--| anterior Deposition efficiency= ',(real(depo1)/real(nline))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| interior meatus Deposition efficiency= ',(real(depo2)/real(nline))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| middle meatus Deposition efficiency= ',(real(depo3)/real(nline))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| olfactory Deposition efficiency= ',(real(depo4)/real(nline))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| rest of nasal cavity Deposition efficiency= ',(real(depo5)/real(nline))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--|'
       
   end program alya_deposition

