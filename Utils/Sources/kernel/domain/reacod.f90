subroutine reacod(itask)
  !-----------------------------------------------------------------------
  !****f* Domain/reacod
  ! NAME
  !    reacod
  ! DESCRIPTION
  !    Reads codes on nodes and boundaries
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_kermod
  use def_inpout
  use mod_iofile
  use mod_chktyp, only : check_type
  use mod_ecoute, only : ecoute
  use mod_outfor, only : outfor
  use mod_memory, only : memory_size
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,iboun,icode,idofn,iparb,ibopo,dummi,iauxi
  integer(ip)             :: ndofn,ncode,ibves,kpoin,ierro,nnand,iword
  integer(ip)             :: jcode,icono,nmcod,kcode(mcono),ivcod,ivcob
  integer(ip)             :: kfl_funno_tmp,kfl_funbo_tmp,ifunc,ktest_size
  integer(ip)             :: ntotn,mcono_tmp
  character(20)           :: messa
  integer(ip)             :: kfl_fixrs_tmp
  character(5)            :: wfixrs,wfname,wtag
  integer(ip), pointer    :: kfl_codes_tmp(:,:)

  ierro = 0

  if( itask == 1 .or. itask == -1 .or. itask == 100 ) then

     !-------------------------------------------------------------------
     !
     ! Master reads node codes
     !
     !-------------------------------------------------------------------

     !if( kfl_icodn == 0 ) call runend('REACOD: CODE CONDITIONS ARE ASSIGNED BUT THERE IS NO NODAL CODE')

     call ecoute('reacod')

     ktest_size = size(tncod,KIND=ip)
!     if (ktest_size == 0) &
!          call runend('REACOD: NODAL CODES BUT CONDITIONS ON BOUNDARIES. MISSING EXTRAPOLATE MAYBE...')
     if ( .not. associated(tncod) ) call runend('REACOD: NODAL CODES BUT CONDITIONS ON BOUNDARIES. MISSING EXTRAPOLATE MAYBE...')

     ncode = 0
     ndofn = tncod(1) % ndofn

     if( itask == 100 ) call outfor(84_ip,momod(modul)%lun_outpu,' ')

     do while( words(1) /= 'ENDCO' )

        nnand = 0
        if(exists('&    ')) then
           !
           ! Count number of multiple codes
           !
           nnand=0
           do iword=1,maxwp
              if(words(iword)=='&    ') nnand=nnand+1
           end do
           if(mcono==1) then
              call runend('REACOD: MULTIPLE CODES FOR NODES NOT POSSIBLE')
           end if
        end if

        if( words(2) /= 'OFF  ' ) then

           ncode = ncode + 1
           nmcod = nnand+1
           !Check that ncode <= mcodc or die 
           if  ( ncode > size(tncod(1) % l,1,ip) ) call runend('REACOD: NUMBER OF LINES IN BOUNDARY CONDITIONS SECTION EXCEEDS THE VALUE GIVEN  BY mcodc. MODIFY reacod.f90 AND RECOMPILE.')


           !
           ! Axes: local or global
           !
           kfl_fixrs_tmp = 0
           if( exists('AXES ') ) then
              wfixrs = getcha('AXES ','     ','#Axes')
              if( wfixrs == 'LOCAL' ) then
                 kfl_fixrs_tmp = -1
              else
                 kfl_fixrs_tmp =  0
              end if
           else
              if( nnand == 0 ) then
                 if( param(1) < 0 ) then
                    kfl_fixrs_tmp = int( param(5) )
                 else
                    kfl_fixrs_tmp = int( param(4+ndofn) )
                 end if
              else
                 if( param(1) < 0 ) then
                    kfl_fixrs_tmp = int( param(3+nmcod+1) )
                 else
                    kfl_fixrs_tmp = int( param(3+nmcod+ndofn) )
                 end if
              end if
           end if
           !
           ! KFL_FUNNO_TMP > 0: Time function
           ! KFL_FUNNO_TMP < 0: Space/time function
           !
           ! Time function
           !
           kfl_funno_tmp = 0
           wfname = '     '
           if( exists('TIMEF') ) then
              wfname  = getcha('TIMEF','     ','#Time Function name')
              do ifunc = 1,number_time_function
                 if( trim(wfname) == trim(time_function(ifunc) % name) ) then
                    kfl_funno_tmp = ifunc
                 end if
              end do
              if( kfl_funno_tmp == 0 ) call runend('REACOD: TIME FUNCTION '//trim(wfname)//' IS NOT DEFINED')
           else if( exists('FIELD') ) then
              kfl_funno_tmp = 1000_ip + getint('FIELD',10000_ip,'#Field')   ! a default value makes no sense here
           else
              if( nnand == 0 ) then
                 if( param(1) < 0 ) then
                    kfl_funno_tmp = int( param(4) )
                 else
                    kfl_funno_tmp = int( param(3+ndofn) )
                 end if
              else
                 if( param(1) < 0 ) then
                    kfl_funno_tmp = int( param(2+nmcod+1) )
                 else
                    kfl_funno_tmp = int( param(2+nmcod+ndofn) )
                 end if
              end if
           end if
           !
           ! Time/Space function
           !
           if( exists('FUNCT') ) then
              wfname  = getcha('FUNCT','     ','#Space/time Function name')
           else if( exists('SPACE') ) then
              wfname  = getcha('SPACE','     ','#Space/time Function name')
           end if
           if( number_space_time_function > 0 ) then 
              if( space_time_function(number_space_time_function) % numfield == 0 ) then  
                 !
                 ! Do not this for space_time fields
                 !
                 if( exists('FUNCT') .or. exists('SPACE') ) then
                    if( kfl_funno_tmp /= 0 ) call runend('REACOD: CANNOT PRESCRIBE BOTH TIME AND SPACE/TIME FUNCTIONS')
                    do ifunc = 1,number_space_time_function
                       if( trim(wfname) == trim(space_time_function(ifunc) % name) ) then
                          kfl_funno_tmp = -ifunc
                       end if
                    end do
                    if( kfl_funno_tmp == 0 ) call runend('REACOD: SPACE/TIME FUNCTION '//trim(wfname)//' IS NOT DEFINED')
                 end if
              end if
           end if
           !
           ! Code tag
           !
           if( exists('TAG  ') ) then
              wtag = getcha('TAG  ','     ','#Function name')
           else
              wtag = ''
           end if
           !
           ! Values on nodes function
           !
           if( exists('VALUE') ) then
              param(1) = - param(1)
              param(3) = getrea('VALUE',1.0_rp,'#VALUES ON NODE FUNCTION')
              call runend('REACOD: VALUE MUST BE REPLACED BY FIELD')
           end if

           tncod(1) % l(ncode) % tag   = wtag
           tncod(1) % l(ncode) % fname = wfname   ! Function name

           if( itask == 100 ) ioutp(1:3) = -100

           if( nnand == 0 ) then
              !
              ! Single code per node
              !
              if( param(1) < 0.0_rp ) then

                 tncod(1) % l(ncode) % lcode(1)  = abs(int( param(1) ))
                 tncod(1) % l(ncode) % kfl_fixno = int( param(2) )
                 tncod(1) % l(ncode) % kfl_value = int( param(3) )
                 tncod(1) % l(ncode) % kfl_funno = kfl_funno_tmp
                 tncod(1) % l(ncode) % kfl_fixrs = kfl_fixrs_tmp

              else

                 if( exists('FIELD') ) then
                    tncod(1) % l(ncode) % lcode(1)  = int( param(1) )
                    tncod(1) % l(ncode) % kfl_fixno = int( param(2) )
                    tncod(1) % l(ncode) % kfl_value = getint('FIELD',1_ip,'#FIELD TO USE AS A BOUNDARY CONDITION')
                    tncod(1) % l(ncode) % kfl_funno = kfl_funno_tmp
                    tncod(1) % l(ncode) % kfl_fixrs = kfl_fixrs_tmp
                 else
                    tncod(1) % l(ncode) % lcode(1)  = int( param(1) )
                    iauxi = int( param(2) )
                    tncod(1) % l(ncode) % kfl_fixno = int( param(2) )
                    tncod(1) % l(ncode) % kfl_value = 0
                    do idofn = 1,ndofn
                       tncod(1) % l(ncode) % bvess(idofn) = param(2+idofn)
                    end do
                    tncod(1) % l(ncode) % kfl_funno = kfl_funno_tmp
                    tncod(1) % l(ncode) % kfl_fixrs = kfl_fixrs_tmp
                 end if
              end if

              if( itask == 100 ) then
                 ioutp(1) = tncod(1) % l(ncode) % lcode(1)
              end if

           else
              !
              ! Multiple code per node
              !
              do jcode = 1,nmcod
                 kcode(jcode) = abs(int(param(jcode)))
              end do
              call heapsorti1(2_ip,nmcod,kcode)

              if( param(1) < 0.0_rp ) then
                 do jcode = 1,nmcod
                    tncod(1) % l(ncode) % lcode(jcode)  = kcode(jcode)
                 end do
                 tncod(1) % l(ncode) % kfl_fixno = int( param(nmcod+1)   )
                 tncod(1) % l(ncode) % kfl_value = int( param(nmcod+2)   )
                 tncod(1) % l(ncode) % kfl_funno = kfl_funno_tmp
                 tncod(1) % l(ncode) % kfl_fixrs = kfl_fixrs_tmp
              else
                 if( exists('FIELD') ) then
                    do jcode = 1,nmcod
                       tncod(1) % l(ncode) % lcode(jcode)  = kcode(jcode)
                    end do
                    tncod(1) % l(ncode) % kfl_fixno = int( param(nmcod+1)   )
                    tncod(1) % l(ncode) % kfl_value = getint('FIELD',1_ip,'#FIELD TO USE AS A BOUNDARY CONDITION')
                    tncod(1) % l(ncode) % kfl_funno = kfl_funno_tmp
                    tncod(1) % l(ncode) % kfl_fixrs = kfl_fixrs_tmp
                 else
                    do jcode = 1,nmcod
                       tncod(1) % l(ncode) % lcode(jcode)  = kcode(jcode)
                    end do
                    tncod(1) % l(ncode) % kfl_fixno = int( param(nmcod+1) )
                    tncod(1) % l(ncode) % kfl_value = 0
                    do idofn = 1,ndofn
                       tncod(1) % l(ncode) % bvess(idofn) = param(1+nmcod+idofn)
                    end do
                    tncod(1) % l(ncode) % kfl_funno = kfl_funno_tmp
                    tncod(1) % l(ncode) % kfl_fixrs = kfl_fixrs_tmp
                 end if
              end if

              if( itask == 100 ) then
                 do jcode = 1,nmcod
                    ioutp(jcode) = tncod(1) % l(ncode) % lcode(jcode)
                 end do

              end if

           end if
           if( itask == 100 ) call outfor(43_ip,momod(modul)%lun_outpu,' ')
           !
           ! Check errors
           !
           if( param(1) < 0 .and. ( .not. READ_AND_RUN() ) ) then
              ivcod = tncod(1) % l(ncode) % kfl_value
              if( ivcod > mfiel .or. ivcod < 0 ) then
                 call outfor(48_ip,0_ip,'WRONG VALUE CODE')
              end if
              if( kfl_field(1,ivcod) == 0 ) then
                 call outfor(48_ip,0_ip,'VALUE CODE WAS NOT DEFINED')
              end if
           end if

        end if

        call ecoute('reacod')
     end do

     tncod(1) % ncode = ncode

  else if( itask == 2 .or. itask == 200 ) then

     !-------------------------------------------------------------------
     !
     ! Master reads boundary codes
     !
     !-------------------------------------------------------------------

     !if( kfl_icodb == 0 ) call runend('REACOD: CODE CONDITIONS ARE ASSIGNED BUT THERE IS NO BOUNDARY CODE')
     ncode = 0
     ndofn = tbcod(1) % ndofn

     !if( ncodb == 0 ) call runend('REACOD: BOUNDARY CODES WERE NOT DEFINED IN DOMAIN')
     call ecoute('reacod')

     if( itask == 200 ) call outfor(85_ip,momod(modul)%lun_outpu,' ')

     do while( words(1) /= 'ENDCO' )

        ncode = ncode + 1
        !
        ! Code tag
        !
        if( exists('TAG  ') ) then
           wtag = getcha('TAG  ','     ','#Function name')
        else
           wtag = ''
        end if
        !
        ! KFL_FUNNO_TMP > 0: Time function
        ! KFL_FUNNO_TMP < 0: Space/time function
        !
        ! Time function
        !
        kfl_funbo_tmp = 0
        wfname = '     '
        if( exists('TIMEF') ) then
           wfname  = getcha('TIMEF','     ','#Time Function name')
           do ifunc = 1,number_time_function
              if( trim(wfname) == trim(time_function(ifunc) % name) ) then
                 kfl_funbo_tmp = ifunc
              end if
           end do
           if( kfl_funbo_tmp == 0 ) call runend('REACOD: TIME FUNCTION '//trim(wfname)//' IS NOT DEFINED')
        else
           if( param(1) < 0 ) then
              kfl_funbo_tmp = int( param(4) )
           else
              kfl_funbo_tmp = int( param(3+ndofn) )
           end if
        end if
        !
        ! Time/Space function
        !
        if( exists('FUNCT') ) then
           wfname  = getcha('FUNCT','     ','#Space/time Function name')
        else if( exists('SPACE') ) then
           wfname  = getcha('SPACE','     ','#Space/time Function name')
        end if
        if( exists('FUNCT') .or. exists('SPACE') ) then
           if( kfl_funbo_tmp /= 0 ) call runend('REACOD: CANNOT PRESCRIBE BOTH TIME AND SPACE/TIME FUNCTIONS')
           do ifunc = 1,number_space_time_function
              if( trim(wfname) == trim(space_time_function(ifunc) % name) ) then
                 kfl_funbo_tmp = -ifunc
              end if
           end do
           if( kfl_funbo_tmp == 0 ) call runend('REACOD: SPACE/TIME FUNCTION '//trim(wfname)//' IS NOT DEFINED')
!print *, "kfl_funbo_tmp", kfl_funbo_tmp
        end if
        !
        ! Values on nodes function
        !
        if( exists('VALUE') ) then
           param(1) = - param(1)
           param(3) = getrea('VALUE',1.0_rp,'#VALUES ON NODE FUNCTION')
        end if

        if( param(1) < 0.0_rp ) then

           tbcod(1) % l(ncode) % lcode     = int( param(1) )
           tbcod(1) % l(ncode) % kfl_fixbo = int( param(2) )
           tbcod(1) % l(ncode) % kfl_value = int( param(3) )
           tbcod(1) % l(ncode) % kfl_funbo = kfl_funbo_tmp
           tbcod(1) % l(ncode) % tag       = wtag

        else

           if( exists('FIELD') ) then
              tbcod(1) % l(ncode) % lcode     = int( param(1) )
              tbcod(1) % l(ncode) % kfl_fixbo = int( param(2) )
              tbcod(1) % l(ncode) % kfl_value = getint('FIELD',1_ip,'#FIELD TO USE AS A BOUNDARY CONDITION')
              tbcod(1) % l(ncode) % kfl_funbo = kfl_funbo_tmp
              tbcod(1) % l(ncode) % tag       = wtag
           else
              tbcod(1) % l(ncode) % lcode     = int( param(1) )
              tbcod(1) % l(ncode) % kfl_fixbo = int( param(2) )
              do idofn = 1,ndofn
                 tbcod(1) % l(ncode) % bvnat(idofn) = param(2+idofn)
              end do
              tbcod(1) % l(ncode) % kfl_funbo = kfl_funbo_tmp
              tbcod(1) % l(ncode) % tag       = wtag
           end if

        end if

        ioutp(1:3) = -100
        if( itask == 200 ) then
           ioutp(1) =  tbcod(1) % l(ncode) % lcode
        end if

        call ecoute('reacod')

     end do

     tbcod(1) % ncode = ncode

  else if ( itask == 4 ) then

     !-------------------------------------------------------------------
     !
     ! Master reads geometrical node codes
     !
     !-------------------------------------------------------------------

     call ecoute('reacod')
     ncode = 0
     ndofn = tgcod(1) % ndofn

     do while(words(1)/='ENDCO')

        if(words(2)/='OFF  ') then

           ncode = ncode + 1

           if( param(1) < 0 ) then
              tgcod(1) % l(ncode) % lcode(1)  = int(abs(param(1)) )
              tgcod(1) % l(ncode) % kfl_value = int( param(2) )
           else
              if( exists('FIELD') ) then
                 tgcod(1) % l(ncode) % lcode(1)  = int(abs(param(1)) )
                 tgcod(1) % l(ncode) % kfl_value = getint('FIELD',1_ip,'#FIELD TO USE AS A BOUNDARY CONDITION')
              else
                 tgcod(1) % l(ncode) % lcode(1)  = int( param(1) )
                 tgcod(1) % l(ncode) % kfl_value = 0
                 do idofn = 1,ndofn
                    tgcod(1) % l(ncode) % bvess(idofn) = param(1+idofn)
                 end do
              end if
           end if

        end if

        call ecoute('reacod')
     end do

     tgcod(1) % ncode = ncode

  else if ( itask == IMPOSE_NODE_CODES .or. itask == IMPOSE_EDGE_CODES ) then

     !-------------------------------------------------------------------
     !
     ! Impose node or edge codes (done only for untagged nodes)
     !
     !-------------------------------------------------------------------

     if( itask == IMPOSE_NODE_CODES ) then
        ntotn         =  npoin
        kfl_codes_tmp => kfl_codno
        mcono_tmp     =  mcono
     else if( itask == IMPOSE_EDGE_CODES ) then
        ntotn         =  meshe(ndivi) % nedge
        kfl_codes_tmp => kfl_coded
        mcono_tmp     =  2
     end if

     ndofn = memory_size(kfl_fixno,1_ip)

     if( itask == IMPOSE_NODE_CODES ) then
        if( ifbop == 1 .and. memory_size(kfl_fixno,2_ip) /= nbopo ) &
             call runend('REACOD: WRONG DIMENSIONS FOR FIXITY ARRAY')
     end if

!!!!!!!!!!!!!!!!!!!!
     ! ESTA LINEA ESTABA ASI PERO ES RARISIMA Y NO SE PARA QUE
     ! ME JODE EN UN PROBLEMA CUANDO HAY UN SOLO ELEMENTO
!!!!     if( ndofn == npoin ) ndofn = 1

     ! Y ME SIGUE JODIENDO ESTA MIERDA.. CHEQUEA LUEGO SI IBVES ES = A NPOIN
     if( ifbes == 1 ) ibves = memory_size(bvess,1_ip)
     ! ASI QUE LO CAGO Y QUE SE CAGUE
     ibves = ntotn + 1
!!!!!!!!!!!!!!!!!!!!

     do ncode = 1,tncod(1) % ncode

        !if (tncod(1) % l(ncode) % cotag  == '') then     !!! esta cosa hay que revisarla, no se puede descomentar!!!

        !
        ! Do it only when the code is untagged
        !
        if( itask == IMPOSE_NODE_CODES ) then
           if(      tncod(1) % l(ncode) % lcode(1) == mcodb+1 ) then
              nmcod = 0
           else if( tncod(1) % l(ncode) % lcode(2) == mcodb+1 ) then
              nmcod = 1
           else if( tncod(1) % l(ncode) % lcode(3) == mcodb+1 ) then
              nmcod = 2
           else
              nmcod = 3
           end if
        else
           if(      tncod(1) % l(ncode) % lcode(1) == mcodb+1 ) then
              nmcod = 0
           else if( tncod(1) % l(ncode) % lcode(2) == mcodb+1 ) then
              nmcod = 1
           else
              nmcod = 2
           end if
        end if
        !
        ! Check if value function exist
        !
        if( itask == IMPOSE_NODE_CODES ) then
           if( tncod(1) % l(ncode) % kfl_value > 0 ) then
              ivcod = tncod(1) % l(ncode) % kfl_value
              !call check_type(bvcod,ivcod,ndofn,npoin)
              call check_type(xfiel,ivcod,ndofn,npoin)
           end if
        end if
        !
        ! Order codes
        !
        do jcode = 1,mcono_tmp
           kcode(jcode) = tncod(1) % l(ncode) % lcode(jcode)
        end do
        call heapsorti1(2_ip,mcono_tmp,kcode)
        icode = tncod(1) % l(ncode) % lcode(1)

        do ipoin = 1,ntotn
           icono = 0
           do jcode = 1,mcono_tmp
              if( kfl_codes_tmp(jcode,ipoin) == abs(kcode(jcode)) ) icono = icono + 1
           end do

           if( icono == mcono_tmp ) then

              kpoin = ipoin
              if( itask == IMPOSE_NODE_CODES ) then
                 ibopo = lpoty(ipoin)
                 if( ifbop == 1 ) kpoin = ibopo
              end if

              if( kpoin == 0 ) then

                 messa = intost(ipoin)
                 ierro = ierro + 1
                 call outfor(2_ip,lun_outpu,&
                      'BOUNDARY CONDITION CANNOT BE IMPOSED ON INTERIOR NODE: '//trim(messa))
              else

                 kfl_fixno(1,kpoin) = tncod(1) % l(ncode) % kfl_fixno

                 call codfix(ndofn,kfl_fixno(1,kpoin))

                 if( ifbes == 1 ) then
                    if( ibves == ntotn ) then
                       if( tncod(1) % l(ncode) % kfl_value == 0 ) then
                          bvess(kpoin,1) = tncod(1) % l(ncode) % bvess(1)
                       else
                          ivcod = tncod(1) % l(ncode) % kfl_value
                          bvess(kpoin,1) = xfiel(ivcod) % a(1,ipoin,1)
                       end if
                    else
                       if( tncod(1) % l(ncode) % kfl_value == 0 ) then
                          do idofn = 1,ndofn
                             bvess(idofn,kpoin) = tncod(1) % l(ncode) % bvess(idofn)
                          end do
                       else
                          ivcod = tncod(1) % l(ncode) % kfl_value
                          do idofn = 1,ndofn
                             bvess(idofn,kpoin) = xfiel(ivcod) % a(idofn,ipoin,1)
                          end do
                       end if
                    end if
                 end if

                 if( iffun /= 0 ) then
                    kfl_funno(kpoin) = tncod(1) % l(ncode) % kfl_funno
                    if( ifloc == 1 ) then
                       ibopo = lpoty(kpoin)
                       if( ibopo /= 0 ) then
                          kfl_fixrs(ibopo) = tncod(1) % l(ncode) % kfl_fixrs
                          !if( kfl_fixrs(ibopo) == -2 .and. nskew > 0 ) call geofix(kpoin,ibopo)
                       end if
                    end if
                 else
                    if( ifloc == 1 ) then
                       ibopo = lpoty(kpoin)
                       if( ibopo /= 0 ) then
                          kfl_fixrs(ibopo) = tncod(1) % l(ncode) % kfl_fixrs
                          !if(kfl_fixrs(ibopo)==-2.and.nskew>0) call geofix(kpoin,ibopo)
                       end if
                    end if
                 end if
              end if

           end if

        end do

        !end if

     end do

  else if ( itask == 20 ) then

     !-------------------------------------------------------------------
     !
     ! Boundary codes
     !
     !-------------------------------------------------------------------

     do ncode = 1,tbcod(1) % ncode
        ivcob = tbcod(1) % l(ncode) % kfl_value
        if( ivcob == 0 ) then
           do iboun = 1,nboun
              if( kfl_codbo(iboun) == tbcod(1) % l(ncode) % lcode ) then
                 kfl_fixbo(iboun) = tbcod(1) % l(ncode) % kfl_fixbo
                 do iparb = 1,tbcod(1) % ndofn
                    bvnat(iparb,iboun) = tbcod(1) % l(ncode) % bvnat(iparb)
                 end do
                 if( iffun == 1 ) then
                    kfl_funbo(iboun) = tbcod(1) % l(ncode) % kfl_funbo
                 end if
              end if
           end do
        else
           do iboun = 1,nboun
              if( kfl_codbo(iboun) == abs(tbcod(1) % l(ncode) % lcode) ) then
                 kfl_fixbo(iboun) = tbcod(1) % l(ncode) % kfl_fixbo
                 do iparb = 1,tbcod(1) % ndofn
                    bvnat(iparb,iboun) = xfiel(ivcob) % a(iparb,iboun,1)
                 end do
                 if( iffun == 1 ) then
                    kfl_funbo(iboun) = tbcod(1) % l(ncode) % kfl_funbo
                 end if
              end if
           end do
        end if
     end do

  else if( itask == 30 ) then

     !-------------------------------------------------------------------
     !
     ! Initial conditions on nodes, only for the codes tagged as INITI
     !
     !-------------------------------------------------------------------


     ! OJO QUE NO ESTA LISTO AUN

     ndofn = memory_size(kfl_fixno,1_ip)
     if( memory_size(kfl_fixno,2_ip) /= nbopo ) &
          call runend('REACOD: WRONG DIMENSIONS FOR INITIAL CONDITIONS FIXITY ARRAY')
     if( ndofn == npoin ) ndofn = 1
     if( ifbes == 1 )     ibves = memory_size(bvess,1_ip)

     do ncode = 1,tncod(1) % ncode

        if(      tncod(1) % l(ncode) % lcode(1) == mcodb+1 ) then
           nmcod = 0
        else if( tncod(1) % l(ncode) % lcode(2) == mcodb+1 ) then
           nmcod = 1
        else if( tncod(1) % l(ncode) % lcode(3) == mcodb+1 ) then
           nmcod = 2
        else
           nmcod = 3
        end if
        !
        ! Order codes
        !
        do jcode = 1,mcono
           kcode(jcode) = tncod(1) % l(ncode) % lcode(jcode)
        end do
        call heapsorti1(2_ip,mcono,kcode)
        icode = tncod(1) % l(ncode) % lcode(1)

        do ipoin = 1,npoin
           icono = 0
           do jcode = 1,mcono
              if( kfl_codes_tmp(jcode,ipoin) == abs(kcode(jcode)) ) icono = icono + 1
           end do

           if( icono == mcono ) then

              ibopo = lpoty(ipoin)
              kpoin = ipoin
              if( ifbop == 1 ) kpoin = ibopo

              if( kpoin == 0 ) then

                 messa = intost(ipoin)
                 ierro = ierro + 1
                 call outfor(2_ip,lun_outpu,&
                      'BOUNDARY CONDITION CANNOT BE IMPOSED ON INTERIOR NODE: '//trim(messa))
              else

                 kfl_fixno(1,kpoin) = tncod(1) % l(ncode) % kfl_fixno
                 call codfix(ndofn,kfl_fixno(1,kpoin))

                 if( ifbes == 1 ) then
                    if( ibves == npoin ) then
                       if( tncod(1) % l(ncode) % kfl_value == 0 ) then
                          bvess(kpoin,1) = tncod(1) % l(ncode) % bvess(1)
                       else
                          ivcod = tncod(1) % l(ncode) % kfl_value
                          bvess(kpoin,1) = xfiel(ivcod) % a(1,ipoin,1)
                       end if
                    else
                       if( tncod(1) % l(ncode) % kfl_value == 0 ) then
                          do idofn = 1,ndofn
                             bvess(idofn,kpoin) = tncod(1) % l(ncode) % bvess(idofn)
                          end do
                       else
                          ivcod = tncod(1) % l(ncode) % kfl_value
                          do idofn = 1,ndofn
                             bvess(idofn,kpoin) = xfiel(ivcod) % a(idofn,ipoin,1)
                          end do
                       end if
                    end if
                 end if

                 if( iffun /= 0 ) then
                    kfl_funno(kpoin) = tncod(1) % l(ncode) % kfl_funno
                    if( ifloc == 1 ) then
                       ibopo = lpoty(kpoin)
                       if( ibopo /= 0 ) then
                          kfl_fixrs(ibopo) = tncod(1) % l(ncode) % kfl_fixrs
                          !if( kfl_fixrs(ibopo) == -2 .and. nskew > 0 ) call geofix(kpoin,ibopo)
                       end if
                    end if
                 else
                    if( ifloc == 1 ) then
                       ibopo = lpoty(kpoin)
                       if( ibopo /= 0 ) then
                          kfl_fixrs(ibopo) = tncod(1) % l(ncode) % kfl_fixrs
                          !if(kfl_fixrs(ibopo)==-2.and.nskew>0) call geofix(kpoin,ibopo)
                       end if
                    end if
                 end if
              end if

           end if
        end do

     end do

  end if

  messa=intost(ierro)
  if(ierro==1) then
     call runend('REACOD: '//trim(messa)//' ERROR HAS BEEN FOUND')
  else if(ierro>=2) then
     call runend('REACOD: '//trim(messa)//' ERRORS HAVE BEEN FOUND')
  end if
  !
  ! Recover original values
  !
  iffun = 0
  ifloc = 0
  ifbop = 0
  ifbes = 1

end subroutine reacod

subroutine geofix(kpoin,ibopo)
  use def_kintyp, only    :  ip
  use def_domain, only    :  kfl_fixno,kfl_fixrs,ndime,lpoin
  implicit none
  integer(ip), intent(in) :: kpoin,ibopo
  integer(ip)             :: idime,kfl_fixn2(2)

  kfl_fixn2(1)=kfl_fixno(    1,kpoin)
  kfl_fixn2(2)=kfl_fixno(ndime,kpoin)

  if(lpoin(ibopo)==0) then
     !
     ! 0 patch
     !
     kfl_fixno(1,kpoin)=kfl_fixn2(1)
     do idime=2,ndime
        kfl_fixno(idime,kpoin)=kfl_fixn2(2)
     end do

  else if(lpoin(ibopo)==1) then
     !
     ! 1 patch
     !
     kfl_fixno(1,kpoin)=kfl_fixn2(1)
     do idime=2,ndime
        kfl_fixno(idime,kpoin)=kfl_fixn2(2)
     end do

  else if(lpoin(ibopo)==2) then
     !
     ! 2 patches
     !
     do idime=1,min(2_ip,ndime)
        kfl_fixno(idime,kpoin)=kfl_fixn2(1)
     end do
     if(ndime==3) kfl_fixno(ndime,kpoin)=kfl_fixn2(2)

  else if(lpoin(ibopo)==3) then
     !
     ! 3 patches: corner of step type
     !
     do idime=1,ndime
        kfl_fixno(idime,kpoin)=kfl_fixn2(1)
     end do

  else if(lpoin(ibopo)==-3) then
     !
     ! 3 patches: corner of bottom type
     !
     do idime=1,ndime
        kfl_fixno(idime,kpoin)=kfl_fixn2(1)
     end do

  end if

end subroutine geofix
