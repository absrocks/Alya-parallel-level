program alya_deposition

  implicit none
  integer(4)            :: one4=1,two4=2
  integer(4)            :: i,ii,jj,kk,j,nline,modu
  real(8)               :: xcoord,ycoord,zcoord,time,r,per
  integer(4)            :: subdom,part,max_part,max_particule
  character(150)        :: fil_name,name
  character(4)          :: cpart,num
  logical               :: dir_e
  integer(4)            :: max_part_type
  real(8)               :: part_id,part_type,kfl_exist,BC
  integer(4)            :: out,depo1,depo2,depo3,depo4
  integer(4)            :: depo5,depo6,depo7,depo8,depo9,numdepo
  integer(4)            :: depo10,depo11,depo12,depo13,depo14
  integer(4), pointer   :: family(:)
  character(150)        :: nunam_pos1,filsa

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
  do ii =1,nline
     read (10,*) time,part_id,part_type,kfl_exist,BC,xcoord,ycoord,zcoord
     if (part_type > max_part_type) then
        max_part_type = part_type
     endif
  enddo
  close(10)
  print*,'MAX_PART_TYPE',max_part_type
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
        depo10=0
        depo11=0
        depo12=0
        depo13=0
        depo14=0
        out=0

        open(unit=10,file=fil_name,status='old')
        read(10,*)
        do ii =1,nline
          read (10,*) time,part_id,part_type,kfl_exist,BC,xcoord,ycoord,zcoord
         !!!! read (10,*) time,part_id,part_type,xcoord,ycoord,zcoord,type,BC
           if (part_type == i) then
              numdepo=numdepo+1

              if (BC==132.0)then !central
                 depo1=depo1+1
              else if (BC==352.0)then !inhaler
                 depo2=depo2+1
              else if (BC==133.0 .or. BC==137.0 .or. BC==138.0 .or. BC==141.0 .or. BC==142.0 .or. BC==149.0   &
                   .or. BC==150.0 .or. BC==151.0 .or. BC==152.0 .or. BC==165.0 .or. BC==166.0 .or. BC==167.0  &
                   .or. BC==168.0 .or. BC==169.0 .or. BC==170.0 .or. BC==171.0 .or. BC==172.0 .or. BC==193.0  &
                   .or. BC==194.0 .or. BC==195.0 .or. BC==196.0 .or. BC==197.0 .or. BC==198.0 .or. BC==199.0  &
                   .or. BC==200.0 .or. BC==201.0 .or. BC==202.0 .or. BC==203.0 .or. BC==204.0 .or. BC==205.0  &
                   .or. BC==206.0 .or. BC==244.0 .or. BC==245.0 .or. BC==246.0 .or. BC==247.0 .or. BC==248.0  &
                   .or. BC==249.0 .or. BC==250.0 .or. BC==251.0 .or. BC==252.0 .or. BC==253.0 .or. BC==254.0  &
                   .or. BC==255.0 .or. BC==294.0 .or. BC==295.0 .or. BC==296.0 .or. BC==297.0 .or. BC==298.0  &
                   .or. BC==299.0 .or. BC==300.0  .or. BC==301.0 .or. BC==302.0 .or. BC==303.0 )then 
                   depo3=depo3+1        !lll
                else if (BC==138.0 .or. BC==143.0 .or. BC==144.0 .or. BC==153.0 .or. BC==154.0 .or. BC==155.0  &
                   .or. BC==156.0 .or. BC==173.0 .or. BC==174.0 .or. BC==175.0 .or. BC==176.0 .or. BC==177.0   &
                   .or. BC==178.0 .or. BC==179.0 .or. BC==180.0 .or. BC==207.0 .or. BC==208.0 .or. BC==209.0   &
                   .or. BC==210.0 .or. BC==211.0 .or. BC==212.0 .or. BC==213.0 .or. BC==214.0 .or. BC==215.0   &
                   .or. BC==216.0 .or. BC==217.0 .or. BC==218.0 .or. BC==219.0 .or. BC==220.0 .or. BC==256.0   &
                   .or. BC==257.0 .or. BC==258.0 .or. BC==259.0 .or. BC==260.0 .or. BC==261.0 .or. BC==262.0   &
                   .or. BC==263.0 .or. BC==264.0 .or. BC==265.0 .or. BC==266.0 .or. BC==267.0 .or. BC==304.0 .or. BC==305.0 )then 
                   depo4=depo4+1        !lul
                else if (BC==353.0)then 
                   depo5=depo5+1       !mouth
                else if (BC==139.0 .or. BC==145.0 .or. BC==146.0 .or. BC==159.0 .or. BC==160.0 .or. BC==185.0 .or. BC==186.0  &
                    .or. BC==187.0 .or. BC==188.0 .or. BC==227.0 .or. BC==228.0 .or. BC==229.0 .or. BC==230.0 .or. BC==231.0  &
                    .or. BC==232.0 .or. BC==233.0 .or. BC==234.0 .or. BC==235.0 .or. BC==270.0 .or. BC==271.0 .or. BC==272.0  &
                    .or. BC==273.0 .or. BC==274.0 .or. BC==275.0 .or. BC==276.0 .or. BC==277.0 .or. BC==278.0 .or. BC==279.0  &
                    .or. BC==306.0 .or. BC==307.0 .or. BC==308.0 .or. BC==309.0 .or. BC==310.0 .or. BC==311.0 .or. BC==312.0  &
                    .or. BC==313.0 .or. BC==314.0 .or. BC==315.0 .or. BC==316.0 .or. BC==317.0 .or. BC==328.0 .or. BC==329.0  &
                    .or. BC==330.0 .or. BC==331.0 .or. BC==332.0 .or. BC==333.0 .or. BC==334.0 .or. BC==335.0 .or. BC==336.0  &
                    .or. BC==337.0 .or. BC==338.0 .or. BC==339.0 .or. BC==346.0 .or. BC==347.0 .or. BC==348.0 .or. BC==349.0  &
                    .or. BC==350.0 .or. BC==351.0 )then
                   depo6=depo6+1       !rll
                else if (BC==145.0 .or. BC==146.0 .or. BC==159.0 .or. BC==185.0 .or. BC==186.0 .or. BC==227.0 .or. BC==228.0 & 
                   .or. BC==229.0 .or. BC==230.0 .or. BC==270.0 .or. BC==271.0 .or. BC==272.0 .or. BC==273.0 .or. BC==306.0  &
                   .or. BC==307.0 .or. BC==308.0 .or. BC==309.0 .or. BC==310.0 .or. BC==311.0 .or. BC==312.0 .or. BC==313.0  &
                   .or. BC==328.0 .or. BC==329.0 .or. BC==330.0 .or. BC==331.0 .or. BC==332.0 .or. BC==333.0 .or. BC==334.0  &
                   .or. BC==335.0 .or. BC==336.0 .or. BC==337.0 .or. BC==346.0 .or. BC==347.0 .or. BC==348.0 .or. BC==349.0  &
                   .or. BC==350.0 .or. BC==351.0 )then
                   depo7=depo7+1       !rml
                else if (BC==134.0 .or. BC==139.0 .or. BC==140.0 .or. BC==147.0 .or. BC==148.0 .or. BC==161.0 .or. BC==162.0 &
                   .or. BC==163.0 .or. BC==164.0 .or. BC==189.0 .or. BC==190.0 .or. BC==191.0 .or. BC==192.0 .or.  BC==236.0 &
                   .or. BC==237.0 .or. BC==238.0 .or. BC==239.0 .or. BC==240.0 .or. BC==241.0 .or. BC==242.0 .or. BC==243.0  &
                   .or. BC==280.0 .or. BC==281.0 .or. BC==282.0 .or. BC==283.0 .or. BC==284.0 .or. BC==285.0 .or. BC==286.0  &
                   .or. BC==287.0 .or. BC==288.0 .or. BC==289.0 .or. BC==290.0 .or. BC==291.0 .or. BC==292.0 .or. BC==293.0  &
                   .or. BC==318.0 .or. BC==319.0 .or. BC==320.0 .or. BC==321.0 .or. BC==322.0 .or. BC==323.0 .or. BC==324.0  &
                   .or. BC==325.0 .or. BC==326.0 .or. BC==327.0 .or. BC==340.0 .or. BC==341.0 .or. BC==342.0 .or. BC==343.0  &
                   .or. BC==344.0 .or. BC==345.0  )then
                 depo8=depo8+1       !rul
              else if (BC==354.0)then
                 depo9=depo9+1       !upper airways
              else if (BC<132.0)then
                 if(  BC== 123 .or. BC== 124 .or. BC==125 .or. BC==126 .or. BC==127 .or. BC==128 .or. BC==129 .or. BC==130   .or. &
                      BC==77  .or. BC== 78.or. BC==79.or. BC==80 .or. BC==81 .or. BC==82 .or. BC==83 .or. BC==84 .or. BC==85 .or. &
                      BC==86  .or. BC==87 .or. BC==88 .or. BC==89 .or. BC==90 .or. BC==91 .or. BC==131 )then
                    depo10=depo10 + 1 !!!L3, Lul
                    out = out + 1
                  else if( BC==92  .or. BC==93  .or. BC==94  .or. BC==95  .or. BC==96   .or. BC==97 .or. BC==98  .or. BC==99  .or. &
                           BC==100 .or. BC==101 .or. BC==102 .or. BC==103 .or. BC==104 .or. BC==105 .or. BC==106 .or. BC==107 .or. &
                           BC==108 .or. BC==109 .or. BC==110 .or. BC==111 .or. BC==112 .or. BC==113 .or. BC==114 .or. BC==115 .or. &
                           BC==116 .or. BC==117 .or. BC==118 .or. BC==119 .or. BC==120 .or. BC==121 .or. BC==122 )then
                     depo11 = depo11 + 1  !!!L6, lll
                     out = out + 1
                  else if( BC==34  .or. BC==35 .or. BC==36 .or. BC== 37.or. BC==38 .or. BC==39 .or. BC==40 .or. BC==42 .or. &
                       BC== 42 .or. BC==43 ) then
                     depo12 = depo12 +1  !!!R4 rml
                     out = out + 1
                  else if( BC==44  .or. BC==45 .or. BC==46 .or. BC==47 .or. BC==48 .or. BC==49 .or. BC==50 .or. BC==51 .or. &
                           BC==52  .or. BC==53 .or. BC==54 .or. BC==55 .or. BC==56 .or. BC==57 .or. BC==58 .or. BC==59 .or. &
                           BC==60  .or. BC==61 .or. BC==62 .or. BC==63 .or. BC==64 .or. BC==65 .or. BC==66 .or. BC==67 .or. &
                           BC==68  .or. BC==69 .or. BC==70 .or. BC==71 .or. BC==72 .or. BC==73 .or. BC==74 .or. BC==75 .or. &
                           BC==76) then
                     depo13 =depo13 + 1 !!!R6 rll
                     out=out+1
                  else
                     out = out + 1
                     print*,'otras tapitas'
                     depo14 =depo14 + 1 !! rul
                  end if
              endif
           endif

        end do
        close(10)

        write(6,*) '--|'
        write(6,*) '--| number of deposited particle in the central   = ',depo1
        write(6,*) '--| number of deposited particle in the inhaler   = ',depo2
        write(6,*) '--| number of deposited particle in the LLL       = ',depo3
        write(6,*) '--| number of deposited particle in the LUL       = ',depo4
        write(6,*) '--| number of deposited particle in the mouth     = ',depo5
        write(6,*) '--| number of deposited particle in the RLL       = ',depo6
        write(6,*) '--| number of deposited particle in the RML       = ',depo7
        write(6,*) '--| number of deposited particle in the RUL       = ',depo8
        write(6,*) '--| number of deposited particle in the UA        = ',depo9
        write(6,*) '--| number of particle going throught the outlet  = ',out
        write(6,*) '--| total number of particle deposited            = ', &
             depo1+depo2+depo3+depo4+depo5+depo6+depo7+depo8+depo9
        write(6,*) '--|'

        write(6,*) '--| central Deposition efficiency   = ',(real(depo1)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| inhaler Deposition efficiency   = ',(real(depo2)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| LLL Deposition efficiency       = ',(real(depo3)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| LUL Deposition efficiency       = ',(real(depo4)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| mouth Deposition efficiency     = ',(real(depo5)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| RLL Deposition efficiency       = ',(real(depo6)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| RML Deposition efficiency       = ',(real(depo7)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| RUL Deposition efficiency       = ',(real(depo8)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| UA Deposition efficiency        = ',(real(depo9)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--| output efficiency               = ',(real(out)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--|'
        write(6,*) '--| output lul                      = ',(real(depo10)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--|'
        write(6,*) '--| output lll                      = ',(real(depo11)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--|'
        write(6,*) '--| output rml                      = ',(real(depo12)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--|'
        write(6,*) '--| output rll                      = ',(real(depo13)/real(numdepo))*100,'%' 
        write(6,*) '--|'
        write(6,*) '--|'
        write(6,*) '--| output rul                      = ',(real(depo14)/real(numdepo))*100,'%' 
        write(6,*) '--|'

     end do
   end program alya_deposition



   !mouth
   !ua
   !central
   !rul
   !rml
   !rll
   !lul
   !lll
   !outlets going to rul:
     !b1r
     !b2r
     !b3r

   !outlets going to rml
     !b4r
     !b5r

   !outlets going to rll
     !b6r
     !b7r
     !b8r
     !b9r
     !b10r

   !outlets going to lul
     !b12l
     !b3l
     !b4l
     !b5l

   !outlets going to lll
     !b6l
     !b8l
     !b9l
     !b10l
