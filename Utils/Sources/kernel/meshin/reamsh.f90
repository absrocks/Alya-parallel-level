subroutine reamsh()
  !-----------------------------------------------------------------------
  !****f* Domain/reamsh
  ! NAME
  !    reamsh
  ! DESCRIPTION
  !    Read cartesian mesh
  ! OUTPUT
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_inpout
  use def_meshin
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip)          :: isour

  call ecoute('reamsh')
  rsuni      = 0.0_rp
  rscal      = 1.0_rp
  mleve      = 8
  nsour      = 0
  nvoxx      = 0
  nvoxy      = 0
  nvoxz      = 0
  kfl_ifbox  = 0
  boxin(1,1) = 0.0_rp
  boxin(2,1) = 0.0_rp
  boxin(3,1) = 0.0_rp
  boxin(1,2) = 0.0_rp
  boxin(2,2) = 0.0_rp
  boxin(3,2) = 0.0_rp

  do while(words(1)/='ENDME')

     if(words(1)=='PARAM') then
        !
        ! Parameters
        !
        call ecoute('reamsh')
        do while(words(1)/='ENDPA')

           if(words(1)=='UNIFO') then
              rsuni=getrea('UNIFO',0.0_rp,'#Uniform size')

           else if(words(1)=='SCALI') then
              rscal=getrea('SCALI',1.0_rp,'#Scaling')

           else if(words(1)=='MAXIM') then
              mleve=getint('MAXIM',0_ip,'#Maximum nunber of levels')

           else if(words(1)=='SOURC') then
              nsour=getint('SOURC',0_ip,'#Number of sources')
              call memmsh(1_ip)

           else if(words(1)=='NUMBE') then
              nvoxx = int(param(1))
              nvoxy = int(param(2))
              nvoxz = int(param(3))

           else if(words(1)=='BOXSI') then
              kfl_ifbox    = 1
              boxin(1,1) = param(1)
              boxin(2,1) = param(2)
              boxin(3,1) = param(3)
              boxin(1,2) = param(4)
              boxin(2,2) = param(5)
              boxin(3,2) = param(6)

           end if
           call ecoute('reamsh')
        end do

     else if(words(1)=='SOURC') then
        !
        ! Sources
        !
        call ecoute('reamsh')
        if(nsour>0) then
           isour=0
           do while(words(1)/='ENDSO')
              isour=isour+1
              if(isour>nsour) &
                   call runend('REAMSH: TOO MANY SOURCES IN SOURCE FIELD')
              rsgeo(1,1,isour)=param(1)
              rsgeo(2,1,isour)=param(2)
              rsgeo(3,1,isour)=param(3)
              call ecoute('reamsh')
              rsgeo(1,2,isour)=param(1)
              rsgeo(2,2,isour)=param(2)
              rsgeo(3,2,isour)=param(3)
              call ecoute('reamsh')
              rsgeo(1,3,isour)=param(1)
              rsgeo(2,3,isour)=param(2)
              rsgeo(3,3,isour)=param(3)
              call ecoute('reamsh')
              rsour(1,isour)=param(1)
              rsour(2,isour)=param(2)
              rsour(3,isour)=param(3)
              call ecoute('reamsh')
           end do
        else
           do while(words(1)/='ENDSO')
              call ecoute('reamsh')
           end do

        end if

     end if

     call ecoute('reamsh')

  end do
  !
  ! Check errors
  !
  !if(nvoxx==0.or.nvoxy==0.or.nvoxz==0)&
  !     call runend('REAMSH: WRONG NUMBER OF BINS')
  if(  kfl_ifbox==1.and.( &
       boxin(1,2)<=boxin(1,1).or.&
       boxin(2,2)<=boxin(2,1).or.&
       boxin(3,2)<=boxin(3,1)))&
       call runend('REAMSH: WRONG BOUNDING BOX')

end subroutine reamsh
