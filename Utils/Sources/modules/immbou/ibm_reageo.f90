subroutine ibm_reageo()
  !-----------------------------------------------------------------------
  !****f* Immbou/ibm_reageo
  ! NAME
  !    ibm_reageo
  ! DESCRIPTION
  !    Read IB
  ! OUTPUT
  ! USED BY
  !    ibm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  use def_inpout
  use mod_memchk
  use mod_iofile
  use mod_ecoute,             only : ecoute
  use mod_read_domain_arrays, only : read_domain_arrays_types
  implicit none
  integer(ip)   :: ipoin,ipara,idime,iimbo,dummi,kfl_icgns,knode,inodb
  integer(ip)   :: ktyib,kboun,iboun,iblty,ielem,ielty,inode,kelem,ktype
  integer(ip)   :: kboel
  character(20) :: mess1

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_icgns = 0
     ktyib     = 0                   ! IB # of nodes is not calculated
     ktype     = 0                   ! IB # of nodes is not calculated
     kboel     = 0                   ! Connected elements are given in boundary list
     !
     ! Reach the section 
     !
     call ecoute('ibm_reageo')
     do while( words(1) /= 'GEOME' )
        call ecoute('ibm_reageo')
     end do
     call ecoute('reageo')
     !
     ! Allocate memory
     !
     call ibm_memall(1_ip)
     !
     ! Allocate model memory
     !
     do iimbo = 1,nimbo
        if( imbou(iimbo) % kfl_model /= 0 ) then
           igene = iimbo
           call ibm_memall(4_ip)
        end if
     end do

     do while( words(1) /= 'ENDGE') 
        !
        ! Postprocess
        !
        if( words(1) == 'IB   ' .and. nimbo > 0 ) then

           iimbo = getint('NUMBE',0_ip,'#IB NUMBER')
           if( iimbo < 1 .or. iimbo > nimbo ) call runend('REAGEO: WRONG IB NUMBER')

           do while( words(1) /= 'ENDIB')

              if( words(1) == 'MASS ' ) then
                 imbou(iimbo) % massa = param(1)
                 if( param(1) <= 0.0_rp ) call runend('REAGEO: WRONG PARTICLE MASS')
              else if( words(1) == 'DENSI' ) then
                 imbou(iimbo) % densi = param(1)
                 if( param(1) <= 0.0_rp ) call runend('REAGEO: WRONG PARTICLE DENSITY')
              else if( words(1) == 'MOMEN' ) then
                 if( param(1) < 0.0_rp ) call runend('REAGEO: WRONG PARTICLE MOMENTUM OF INERTIA')
                 imbou(iimbo) % momin      = 0.0_rp  ! initialize to zero in case not all are defined
                 do ipara = 1,nnpar
                    imbou(iimbo) % momin(ipara) = param(ipara)
                 end do
              else if( words(1) == 'POSGR' ) then
                 if( param(1) < -0.5e12_rp ) call runend('REAGEO: WRONG PARTICLE CENTER OF GRAVITY')
                 imbou(iimbo) % posgr      = 0.0_rp  ! initialize to zero in case not all are defined
                 do ipara = 1,nnpar
                    imbou(iimbo) % posgr(ipara) = param(ipara)
                 end do
              else if( words(1) == 'POSIL' ) then
                 do ipara = 1,3
                    imbou(iimbo) % posil(ipara,1) = param(ipara)
                 end do
              else if( words(1) == 'VELOL' ) then
                 do ipara = 1,3
                    imbou(iimbo) % velol(ipara,1) = param(ipara)
                 end do
              else if( words(1) == 'ACCEL' ) then
                 do ipara = 1,3
                    imbou(iimbo) % accel(ipara,1) = param(ipara)
                 end do
              else if( words(1) == 'POSIA' ) then
                 do ipara = 1,3
                    imbou(iimbo) % posia(ipara,1) = param(ipara)
                 end do
              else if( words(1) == 'VELOA' ) then
                 do ipara = 1,3
                    imbou(iimbo) % veloa(ipara,1) = param(ipara)
                 end do
              else if( words(1) == 'ACCEA' ) then
                 do ipara = 1,3
                    imbou(iimbo) % accea(ipara,1) = param(ipara)
                 end do

              else if( words(1) == 'COORD' ) then 
                 !
                 ! COOIB: IB coordinates
                 !                     
                 if( imbou(iimbo) % kfl_typeb == - 1 ) then
                    do ipoin = 1,imbou(iimbo) % npoin
                       read(nunit,*,err=3) dummi,(imbou(iimbo) % coord(idime,ipoin),idime=1,ndime)                    
                    end do
                 else
                    do ipoin = 1,imbou(iimbo) % npoib
                       read(nunit,*,err=3) dummi,(imbou(iimbo) % cooib(idime,ipoin),idime=1,ndime)                    
                    end do
                 end if
                 call ecoute('reageo')
                 if( words(1) /= 'ENDCO'.and.words(1) /= 'END  ')&
                      call runend('REAGEO: WRONG IB COORDINATE FIELD')

              else if( words(1) == 'TYPES' ) then

                 if( imbou(iimbo) % kfl_typeb == - 1 ) then
                    !
                    ! LTYPE: IB element types 
                    !           
                    call read_domain_arrays_types(kfl_icgns,imbou(iimbo) % nelem,ktype,imbou(iimbo) % ltype,lexib)
                 else
                    !
                    ! LTYIB: IB boundary types
                    !           
                    call read_domain_arrays_types(kfl_icgns,imbou(iimbo) % nboib,ktyib,imbou(iimbo) % ltyib,lexib)
                 end if

              else if( words(1) == 'ELEME' ) then
                 !
                 ! LNOIB: IB boundaries
                 !
                 call ecoute('reageo')
                 if( ktype == imbou(iimbo) % nelem ) then
                    kelem = 0
                    do while( words(1) /= 'ENDEL'.and.words(1) /= 'END  ')
                       ielem = int(param(1))
                       kelem = kelem+1
                       if( ielem < 0 .or. ielem > imbou(iimbo) % nelem) then
                          mess1 = intost(iimbo)
                          call runend('REAGEO: WRONG ELEMENT NUMBER READING LNOIB FOR PARTICLE '//trim(mess1))
                       end if
                       do inode = 1,nnode(imbou(iimbo) % ltype(ielem))
                          imbou(iimbo) % lnods(inode,ielem) = int(param(inode+1))
                       end do
                       call ecoute('reageo')
                    end do
                    if( kelem /= imbou(iimbo) % nelem ) &
                         call runend('REAGEO: WRONG TOTAL NUMBER OF ELEMENTS')                
                 else
                    kelem = 0
                    do while( words(1) /= 'ENDEL'.and.words(1) /= 'END  ')
                       ielem = int(param(1))
                       kelem = kelem+1
                       if( ielem < 0 .or. ielem > imbou(iimbo) % nelem ) then
                          mess1 = intost(iimbo)
                          call runend('REAGEO: WRONG ELEMENT NUMBER READING LNOIB FOR PARTICLE '//trim(mess1))
                       end if
                       knode = nnpar-1
                       call fintyp(ndime,knode,ielty)
                       imbou(iimbo) % ltype(ielem) = ielty
                       do inode = 1,nnode(ielty)
                          imbou(iimbo) % lnods(inode,ielem) = int(param(inode+1))
                       end do
                       call ecoute('reageo')
                    end do
                    if(kelem/=imbou(iimbo) % nelem) &
                         call runend('REAGEO: WRONG TOTAL NUMBER OF ELEMENTS') 

                 end if

              else if( words(1) == 'BOUND' ) then
                 !
                 ! LNOIB: IB boundaries
                 !
                 if( exists('ELEME' ) ) kboel = -1
                 call ecoute('reageo')
                 if( ktyib == imbou(iimbo) % nboib ) then
                    kboun = 0
                    do while( words(1) /= 'ENDBO'.and.words(1) /= 'END  ')
                       iboun = int(param(1))
                       kboun = kboun+1
                       if( iboun < 0 .or. iboun > imbou(iimbo) % nboib ) then
                          mess1 = intost(iimbo)
                          call runend('REAGEO: WRONG BOUNDARY NUMBER READING LNOIB FOR PARTICLE '//trim(mess1))
                       end if
                       do inodb = 1,nnode(imbou(iimbo) % ltyib(iboun))
                          imbou(iimbo) % lnoib(inodb,iboun) = int(param(inodb+1))
                       end do
                       call ecoute('reageo')
                    end do
                    if( kboun /= imbou(iimbo) % nboib ) &
                         call runend('REAGEO: WRONG TOTAL NUMBER OF BOUNDARIES')                     
                 else
                    kboun = 0
                    do while( words(1) /= 'ENDBO' .and. words(1) /= 'END  ' )
                       iboun = int(param(1))
                       kboun = kboun+1
                       if( iboun < 0 .or. iboun > imbou(iimbo) % nboib ) then
                          mess1 = intost(iimbo)
                          call runend('REAGEO: WRONG BOUNDARY NUMBER READING LNOIB FOR PARTICLE '//trim(mess1))
                       end if
                       knode = nnpar-1+kboel
                       call fintyp(ndimb,knode,iblty)
                       imbou(iimbo) % ltyib(iboun)=iblty
                       do inodb = 1,nnode(iblty)
                          imbou(iimbo) % lnoib(inodb,iboun) = int(param(inodb+1))
                       end do
                       call ecoute('reageo')
                    end do
                    if(kboun/=imbou(iimbo) % nboib) &
                         call runend('REAGEO: WRONG TOTAL NUMBER OF BOUNDARIES') 
                 end if

              else if( words(1) == 'MODEL' ) then
                 !
                 ! Read model geometry and definition
                 !
                 call ibm_model_wind_turbine(2_ip,iimbo)
                 call ibm_model_wind_turbine(3_ip,iimbo)

              end if
              call ecoute('reageo')
           end do

        end if
        call ecoute('reageo')
     end do

  end if

  return

1 call runend('REAGEO: WRONG IB NUMBER OF NODES FOR ELEMENT '//trim(intost(ielem)))
2 call runend('REAGEO: WRONG IB ELEMENT CONNECTIVITY FOR ELEMENT '//trim(intost(ielem)))
3 call runend('REAGEO: WRONG IB COORDINATES FOR NODE '//trim(intost(ipoin)))

end subroutine ibm_reageo
