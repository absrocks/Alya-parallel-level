subroutine posdef(itask,ivara) 
  !-----------------------------------------------------------------------
  !****f* master/posdef
  ! NAME
  !    posdef
  ! DESCRIPTION
  !    This subroutine deos the following:
  !    ITASK = 0 ... Initialize the solver type
  !    ITASK = 1 ... Bridge between modules and parall service
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master 
  use def_inpout
  use def_postpr
  use mod_memchk
  use mod_communications, only : PAR_SUM
  use mod_ecoute, only :  ecoute
  implicit none  
  integer(ip), intent(in)    :: itask
  integer(ip), intent(inout) :: ivara
  integer(ip)                :: ivari,imodu,ivarp,ivarw,ivart,ivars,itime
  integer(ip)                :: ipost,nvabs,nvaes,nvans
  integer(4)                 :: istat
  character(5)               :: CNULL='NULL '
  real(rp),    pointer       :: per_vbset(:,:)
  real(rp),    pointer       :: per_veset(:,:)
  
  if( itask == 0_ip ) then

     !-------------------------------------------------------------------
     !
     ! Initialization
     !
     !-------------------------------------------------------------------

     do imodu = 0,mmodu
        allocate( momod(imodu) % postp(1) )
     end do

     do imodu = 0,mmodu
        !
        ! Integer(ip)
        !
        postp => momod(imodu)%postp
        postp(1) % npp_inits = 0                           ! Postprocess initial step
        do ivarp = 1,nvarp
           postp(1) % kfl_oonce(ivarp) = 0                 ! Do not postprocess only once
           postp(1) % npp_stepi(ivarp) = 0                 ! Postprocess step interval
           postp(1) % pos_alrea(ivarp) = 0                 ! Already postprocessed
           postp(1) % vox_stepi(ivarp) = 0                 ! Postprocess step interval
           postp(1) % vox_alrea(ivarp) = 0                 ! Already postprocessed
           postp(1) % nfilt(ivarp)     = 0                 ! Filter 
        end do
        postp(1) %  npp_iniso = 0                          ! Postprocess initial condition
        do ivars = 1,nvars
           postp(1) % npp_setse(ivars) = 0                 ! Postprocess element sets calculation
           postp(1) % npp_setsb(ivars) = 0                 ! Postprocess boundary sets calculation
           postp(1) % npp_setsn(ivars) = 0                 ! Postprocess node sets calculation
           postp(1) % per_setse(ivars) = 0                 ! Postprocess permutation element sets calculation
           postp(1) % per_setsb(ivars) = 0                 ! Postprocess permutation boundary sets calculation
           postp(1) % per_setsn(ivars) = 0                 ! Postprocess permutation node sets calculation
        end do
        postp(1) % npp_stepw = 0                           ! Postprocess witness points at every step
        do ivarw = 1,nvarw
           postp(1) % npp_witne(ivarw) = 0                 ! Postprocess witness points
        end do
        postp(1) % nvaes     = 0                           ! Element  set variables
        postp(1) % nvabs     = 0                           ! Boundary set variables
        postp(1) % nvans     = 0                           ! Node     set variables
        postp(1) % per_nvaes = 0                           ! Element  set variables
        postp(1) % per_nvabs = 0                           ! Boundary set variables
        postp(1) % per_nvans = 0                           ! Node     set variables
        postp(1) % nvawi     = 0                           ! Node     set variables
        postp(1) % lun_setse = 0                           ! Element set unit imodu*10+6
        postp(1) % lun_setsb = 0                           ! Boundary set unit imodu*10+7
        postp(1) % lun_setsn = 0                           ! Node set unit imodu*10+8
        postp(1) % lun_witne = 0                           ! Node set unit imodu*10+8
        postp(1) % ipass     = 0                           ! Set memory allocated and header
        !
        ! Real(rp)
        !
        postp(1) % pos_tinit = 0.0_rp                      ! Postprocess initial time
        do ivarp = 1,nvarp
           do ivart = 1,nvart
              postp(1) % pos_times(ivart,ivarp) = 0.0_rp   ! Postprocess times 
              postp(1) % vox_times(ivart,ivarp) = 0.0_rp   ! Voxel times 
           end do
           postp(1) % pos_perio(ivarp) = 0.0_rp            ! Postprocess time period
           postp(1) % vox_perio(ivarp) = 0.0_rp            ! Voxel time period
        end do
        do ivars = 1,nvars
           do ivari = 1,5
              postp(1) % paese(ivari,ivars) = 0.0_rp       ! Element set parameters  
              postp(1) % pabse(ivari,ivars) = 0.0_rp       ! Boundary sets parameters  
              postp(1) % panse(ivari,ivars) = 0.0_rp       ! Node set parameters  
           end do
        end do
        !
        ! Character(5)
        !
        do ivars = 1,nvars
           postp(1) % woese(ivars) = cnull                 ! Name of the element set variables
           postp(1) % wobse(ivars) = cnull                 ! Name of the boundary set variables
           postp(1) % wonse(ivars) = cnull                 ! Name of the node set variables          
        end do
        do ivarw = 1,nvarw
           postp(1) % wowit(ivarw) = cnull                 ! Name of the witness 
        end do
        do ivarp = 1,nvarp
           postp(1) % wopos(1,ivarp) = cnull               ! Name of postprocess variable
           postp(1) % wopos(2,ivarp) = 'SCALA'             ! Character of postprocess variable
           postp(1) % wopos(3,ivarp) = 'NPOIN'             ! NPOIN type by default
        end do
        !
        ! Pointers
        !
        nullify( postp(1) % veset )
        nullify( postp(1) % vbset )
        nullify( postp(1) % vnset )
        nullify( postp(1) % witne )
     end do

  else if( itask == 1_ip ) then

     !-------------------------------------------------------------------
     !
     ! Used for Parall service
     !
     !-------------------------------------------------------------------

     !
     ! Integer(ip)
     !
     postp => momod(modul)%postp
     call iexcha(postp(1) % npp_inits)                           ! Postprocess initial step
     do ivarp = 1,nvarp
        call iexcha(postp(1) % kfl_oonce(ivarp))                 ! Postprocess only once
        call iexcha(postp(1) % npp_stepi(ivarp))                 ! Postprocess step interval
        call iexcha(postp(1) % pos_alrea(ivarp))                 ! Already postprocessed
        call iexcha(postp(1) % vox_stepi(ivarp))                 ! Postprocess step interval
        call iexcha(postp(1) % vox_alrea(ivarp))                 ! Already postprocessed
        call iexcha(postp(1) % nfilt(ivarp))                     ! Filter
     end do
     call iexcha(postp(1) %  npp_iniso)                          ! Postprocess initial condition
     do ivars = 1,nvars
        call iexcha(postp(1) % npp_setse(ivars))                 ! Postprocess element sets calculation
        call iexcha(postp(1) % npp_setsb(ivars))                 ! Postprocess boundary sets calculation
        call iexcha(postp(1) % npp_setsn(ivars))                 ! Postprocess node sets calculation
        call iexcha(postp(1) % per_setse(ivars))                 ! Postprocess element sets calculation
        call iexcha(postp(1) % per_setsb(ivars))                 ! Postprocess boundary sets calculation
        call iexcha(postp(1) % per_setsn(ivars))                 ! Postprocess node sets calculation
     end do
     call iexcha(postp(1) % npp_stepw)                           ! Postprocess witness points interval
     do ivarw = 1,nvarw
        call iexcha(postp(1) % npp_witne(ivarw))                 ! Postprocess witness points
     end do
     call iexcha(postp(1) % nvaes)                               ! Element  set variables  
     call iexcha(postp(1) % nvabs)                               ! Boundary set variables
     call iexcha(postp(1) % nvans)                               ! Node set variables
     call iexcha(postp(1) % per_nvaes)                           ! Element  set variables  
     call iexcha(postp(1) % per_nvabs)                           ! Boundary set variables
     call iexcha(postp(1) % per_nvans)                           ! Node set variables
     call iexcha(postp(1) % nvawi)                               ! Witness variables
     call iexcha(postp(1) % lun_setse)                           ! Element set unit imodu*10+6
     call iexcha(postp(1) % lun_setsb)                           ! Boundary set unit imodu*10+7
     call iexcha(postp(1) % lun_setsn)                           ! Node set unit imodu*10+8
     call iexcha(postp(1) % lun_witne)                           ! Node set unit imodu*10+8
     call iexcha(postp(1) % ipass)                               ! Set memory allocated and header
     !
     ! Real(rp)
     !
     call rexcha(postp(1) % pos_tinit)                           ! Postprocess initial time
     do ivarp = 1,nvarp
        do ivart = 1,nvart
           call rexcha(postp(1) % pos_times(ivart,ivarp))        ! Postprocess times 
           call rexcha(postp(1) % vox_times(ivart,ivarp))        ! Postprocess times 
        end do
        call rexcha(postp(1) % pos_perio(ivarp))                 ! Postprocess time period 
        call rexcha(postp(1) % vox_perio(ivarp))                 ! Postprocess time period 
     end do
     do ivars = 1,nvars
        do ivari = 1,5
           call rexcha(postp(1) % paese(ivari,ivars))            ! Element set parameters  
           call rexcha(postp(1) % pabse(ivari,ivars))            ! Boundary sets parameters  
           call rexcha(postp(1) % panse(ivari,ivars))            ! Node set parameters  
        end do
     end do
     !
     ! Character(5)
     !
     do ivars = 1,nvars
        if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)  = postp(1) % woese(ivars)
        if(parii==2.and.ISLAVE)  postp(1) % woese(ivars) = parch(nparc+1:nparc+5)
        nparc = nparc+5
        if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)  = postp(1) % wobse(ivars)
        if(parii==2.and.ISLAVE)  postp(1) % wobse(ivars) = parch(nparc+1:nparc+5)
        nparc = nparc+5
        if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)  = postp(1) % wonse(ivars)
        if(parii==2.and.ISLAVE)  postp(1) % wonse(ivars) = parch(nparc+1:nparc+5)
        nparc = nparc+5
     end do
     do ivarw = 1,nvarw
        if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)  = postp(1) % wowit(ivarw)
        if(parii==2.and.ISLAVE)  postp(1) % wowit(ivarw) = parch(nparc+1:nparc+5)
        nparc = nparc+5           
     end do
     do ivarp = 1,nvarp
        if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)    = postp(1) % wopos(1,ivarp) 
        if(parii==2.and.ISLAVE)  postp(1) % wopos(1,ivarp) = parch(nparc+1:nparc+5)
        nparc = nparc+5           
        if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)    = postp(1) % wopos(2,ivarp) 
        if(parii==2.and.ISLAVE)  postp(1) % wopos(2,ivarp) = parch(nparc+1:nparc+5)
        nparc = nparc+5           
        if(parii==2.and.IMASTER) parch(nparc+1:nparc+5)    = postp(1) % wopos(3,ivarp) 
        if(parii==2.and.ISLAVE)  postp(1) % wopos(3,ivarp) = parch(nparc+1:nparc+5)
        nparc = nparc+5           
     end do

     if( nparc > len(parch) ) call runend('POSDEF: TOO MANY CHARACTERS')

  else if( itask == 2_ip ) then

     !-------------------------------------------------------------------
     !
     ! Read postprocess
     !
     !-------------------------------------------------------------------

     if( words(1) == 'POSTP' ) then
        if( exists('ONVOX') ) then
           ipost = 2
        else
           ipost = 1
        end if
        if( words(2) == 'INITI' ) then
           postp(1) % npp_iniso = 1
        else
           do ivarp = 1,nvarp
              if( words(2) == postp(1) % wopos(1,ivarp) ) then
                 if( words(3) == 'STEPS' ) then
                    if( ipost == 1 ) then
                       postp(1) % npp_stepi(ivarp) = &
                            getint('STEPS',1_ip,&
                            '#Postprocess step interval for '// postp(1) % wopos(1,ivarp))
                       if (npp_stepo>-1) postp(1) % npp_stepi(ivarp) = npp_stepo
                    else
                       postp(1) % vox_stepi(ivarp) = &
                            getint('STEPS',1_ip,&
                            '#Postprocess step interval for '// postp(1) % wopos(1,ivarp))                       
                       if (npp_stepo>-1) postp(1) % vox_stepi(ivarp) = npp_stepo
                    end if
                    if(exists('ATTIM')) then
                       if( ipost == 1 ) then
                          do ivart = 1,nvart
                             postp(1) % pos_times(ivart,ivarp) = param(ivart+3)
                          end do
                       else
                          do ivart = 1,nvart
                             postp(1) % vox_times(ivart,ivarp) = param(ivart+3)
                          end do                          
                       end if
                    end if
                 else if(words(3)=='PERIO') then
                    if( ipost == 1 ) then
                       postp(1) % pos_perio(ivarp)= getrea('PERIO',0.0_rp,'#Postprocess time period')
                    else
                       postp(1) % vox_perio(ivarp)= getrea('PERIO',0.0_rp,'#Postprocess time period')
                    end if
                 else
                    if( ipost == 1 ) then
                       postp(1) % npp_stepi(ivarp) = 1
                       if (npp_stepo>-1) postp(1) % npp_stepi(ivarp) = npp_stepo
                    else
                       postp(1) % vox_stepi(ivarp) = 1
                       if (npp_stepo>-1) postp(1) % vox_stepi(ivarp) = npp_stepo
                    end if
                 end if
                 if( exists('FILTE') ) then
                    postp(1) % nfilt(ivarp) = getint('FILTE',1_ip,'#Filter number')        
                 end if
                 if( exists('ONLYO') ) then
                    postp(1) % kfl_oonce(ivarp) = 1
                 end if
              end if
           end do
        end if

     else if( words(1) == 'ELEME' ) then
        !
        ! Element sets
        !
        if( neset == 0 ) then
           do while( words(1) /= 'ENDEL' )
              call ecoute('posdef')
           end do
        else
           ivars = 0
           postp(1) % nvaes = 0
           do while( ivars < nvars )
              ivars = ivars + 1
              if( postp(1) % woese(ivars) == 'NULL ' ) then
                 postp(1) % nvaes = ivars - 1
                 ivars = nvars
              end if
           end do
           nvaes = 0
           do while( words(1) /= 'ENDEL' )
              call ecoute('posdef')
              do ivars = 1,nvars
                 if( words(1) == postp(1) % woese(ivars) ) then
                    nvaes                       = nvaes + 1 
                    postp(1) % per_setse(nvaes) = ivars     
                    postp(1) % npp_setse(ivars) = 1
                    postp(1) % paese(1:5,ivars) = param(2:6)
                 end if
              end do
           end do
           postp(1) % per_nvaes = nvaes
        end if

     else if( words(1) == 'BOUND' ) then
        !
        ! Boundary sets
        !
        if( nbset == 0 ) then
           do while( words(1) /= 'ENDBO' )
              call ecoute('posdef')
           end do
        else
           ivars = 0
           postp(1) % nvabs = 0
           do while( ivars < nvars )
              ivars = ivars + 1
              if( postp(1) % wobse(ivars) == 'NULL ' ) then
                 postp(1) % nvabs = ivars - 1
                 ivars = nvars
              end if
           end do
           nvabs = 0                                     
           do while( words(1) /= 'ENDBO' )
              call ecoute('posdef')
              do ivars = 1,nvars
                 if( words(1) == postp(1) % wobse(ivars) ) then
                    nvabs                       = nvabs + 1 
                    postp(1) % per_setsb(nvabs) = ivars     
                    postp(1) % npp_setsb(ivars) = 1
                    postp(1) % pabse(1:5,ivars) = param(2:6)
                 end if 
              end do
           end do 
           postp(1) % per_nvabs = nvabs
        end if

     else if (words(1) == 'NODES' ) then
        !
        ! Node sets
        !
        if( nnset == 0 ) then
           do while( words(1) /= 'ENDNO' )
              call ecoute('posdef')
           end do
        else
           ivars = 0
           postp(1) % nvans = 0
           do while( ivars < nvars )
              ivars = ivars + 1
              if( postp(1) % wonse(ivars) == 'NULL ' ) then
                 postp(1) % nvans = ivars - 1
                 ivars = nvars
              end if
           end do
           do while( words(1) /= 'ENDNO' )
              call ecoute('posdef')
              do ivars = 1,nvars
                 if( words(1) == postp(1) % wonse(ivars) ) then
                    postp(1) % npp_setsn(ivars) = 1
                    postp(1) % panse(1:5,ivars) = param(2:6)
                 end if
              end do
           end do
        end if

     else if (words(1) == 'WITNE' .and. .not. exists('NUMBE') ) then
        postp(1) % npp_stepw= 1  ! Default: postprocess at all time step
        postp(1) % npp_stepw= getint('STEPS',1_ip,'#When to postprocess witness points')
        !
        ! Witness nodes
        !
        if( nwitn == 0 ) then
           do while( words(1) /= 'ENDWI' )
              call ecoute('posdef')
           end do
        else
           ivarw = 0
           postp(1) % nvawi = 0
           do while( ivarw < nvarw )
              ivarw = ivarw + 1
              if( postp(1) % wowit(ivarw) == 'NULL ' ) then
                 postp(1) % nvawi = ivarw - 1
                 ivarw = nvarw
              end if
           end do
           do while( words(1) /= 'ENDWI' )
              call ecoute('posdef')
              do ivarw = 1,nvarw
                 if( words(1) == postp(1) % wowit(ivarw) ) then
                    postp(1) % npp_witne(ivarw) = 1
                 end if
              end do
           end do
        end if
     else if( words(1) == 'START' ) then
        !
        ! Starting time and step of post-process
        !
        if( exists('STEP ') ) then
           postp(1) % npp_inits = getint('STEP ',0_ip,'#Initial step to start postprocess')  
           if( postp(1) % npp_inits == 0 ) postp(1) % npp_iniso = 1
        end if
        if( exists('TIME ') ) then
           postp(1) % pos_tinit = getrea('TIME ',0.0_rp,'#Initial step to start postprocess')
        end if

     end if

  else if( itask == 3_ip ) then

     !-------------------------------------------------------------------
     !
     ! Allocate memory for sets and witness points
     !
     !-------------------------------------------------------------------

     if( postp(1) % ipass == 0 ) then

        postp(1) % ipass = 1

        if(maxval(postp(1) % npp_setse)>0) then
           allocate(postp(1) % veset(postp(1) % nvaes+1,neset),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'VESET','posdef',postp(1) % veset)
        end if
        if(maxval(postp(1) % npp_setsb)>0) then
           allocate(postp(1) % vbset(postp(1) % nvabs+1,nbset),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'VBSET','posdef',postp(1) % vbset)
           !if(nboib>0) then
           !   allocate(postp(1) % viset(postp(1) % nvabs+1,niset),stat=istat)
           !   call memchk(zero,istat,mem_modul(1:2,modul),'VISET_NSI','nsi_memose',postp(1) % viset)        
           !end if
        end if
        if(maxval(postp(1) % npp_setsn)>0) then
           allocate(postp(1) % vnset(postp(1) % nvans,nnset),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'VNSET','posdef',postp(1) % vnset)
        end if
        if(maxval(postp(1) % npp_witne)>0) then
           allocate(postp(1) % witne(postp(1) % nvawi,nwitn),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'WITNE','posdef',postp(1) % witne)
        end if

        if(maxval(postp(1) % npp_setse)>0)&
             call outset(&
             -1_ip,postp(1) % lun_setse,postp(1) % nvaes,  postp(1) % npp_setse,postp(1) % woese,postp(1) % veset)
        if(maxval(postp(1) % npp_setsb)>0)&
             call outset(&
             -2_ip,postp(1) % lun_setsb,postp(1) % nvabs,  postp(1) % npp_setsb,postp(1) % wobse,postp(1) % vbset)
        if(maxval(postp(1) % npp_setsn)>0)&
             call outset(&
             -3_ip,postp(1) % lun_setsn,postp(1) % nvans-1_ip,postp(1) % npp_setsn,postp(1) % wonse,postp(1) % vnset)
        if(maxval(postp(1) % npp_witne)>0)&
             call outset(&
             -5_ip,postp(1) % lun_witne,postp(1) % nvawi-1_ip,postp(1) % npp_witne,postp(1) % wowit,postp(1) % witne)

     end if

  else if( itask == 4_ip ) then

     !-------------------------------------------------------------------
     !
     ! Write set and witness results
     !
     !-------------------------------------------------------------------

     if( maxval(postp(1) % npp_setse) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
        call outset(&
             1_ip, postp(1) % lun_setse,postp(1) % nvaes,postp(1) % npp_setse,postp(1) % woese,postp(1) % veset)     
        call outset(&
             10_ip,postp(1) % lun_setse,postp(1) % nvaes,postp(1) % npp_setse,postp(1) % woese,postp(1) % veset)
     end if

     if( maxval(postp(1) % npp_setsb) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
        call outset(&
             2_ip, postp(1) % lun_setsb,postp(1) % nvabs,postp(1) % npp_setsb,postp(1) % wobse,postp(1) % vbset)
        call outset(&
             20_ip,postp(1) % lun_setsb,postp(1) % nvabs,postp(1) % npp_setsb,postp(1) % wobse,postp(1) % vbset)
     end if

     if( maxval(postp(1) % npp_setsn) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
        call outset(&
             3_ip, postp(1) % lun_setsn,postp(1) % nvans-1_ip,postp(1) % npp_setsn,postp(1) % wonse,postp(1) % vnset)           
        call outset(&
             30_ip,postp(1) % lun_setsn,postp(1) % nvans-1_ip,postp(1) % npp_setsn,postp(1) % wonse,postp(1) % vnset)
     end if

     if( maxval(postp(1) % npp_witne) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
        if(( mod(ittim, postp(1) % npp_stepw) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
           
           call outset(&
                5_ip, postp(1) % lun_witne,postp(1) % nvawi-1_ip,postp(1) % npp_witne,postp(1) % wowit,postp(1) % witne)           
           call outset(&
                50_ip,postp(1) % lun_witne,postp(1) % nvawi-1_ip,postp(1) % npp_witne,postp(1) % wowit,postp(1) % witne)
        end if
     end if

  else if( itask == 5_ip ) then

     !-------------------------------------------------------------------
     !
     ! Check: write witness or not?
     !
     !-------------------------------------------------------------------
     ivara = 0
     if( maxval(postp(1) % npp_witne) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
        if(( mod(ittim, postp(1) % npp_stepw) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
           ivara = 1
        end if
     end if

  else if( itask == 11 .or. itask == 111 ) then

     !-------------------------------------------------------------------
     !
     ! Initial solution, end of time step, end of the run
     !
     !-------------------------------------------------------------------

     ivari = ivara
     ivara = 0
     kfl_ivari = 0

     if( postp(1) % pos_alrea(ivari) /= 2 ) then
        if( ittyp == ITASK_INITIA .and. kfl_rstar < 2 ) then  ! outputs if INITI or no restart
           !
           ! Initial solution
           !
           if( postp(1) % npp_iniso == 1 .and. postp(1) % npp_stepi(ivari) > 0 ) then  
              ivara = ivari
              if( itask == 11 ) postp(1) % pos_alrea(ivari) = 1 ! To avoid repostprocess when no time step is required
           end if
 
        else if( ittyp == ITASK_ENDTIM ) then
           !
           ! End of a time step
           !
           if( cutim >= postp(1) % pos_tinit ) then
              if( itask == 11 ) postp(1) % pos_alrea(ivari) = 0
              !
              ! At a given time step
              !           
              if( &
                   ittim >= postp(1) % npp_inits .and. &
                   postp(1) % pos_alrea(ivari) == 0 .and. &
                   postp(1) % npp_stepi(ivari) > 0 ) then     
                 if( mod(ittim, postp(1) % npp_stepi(ivari)) == 0 ) then
                    ivara = ivari
                    if( itask == 11 ) postp(1) % pos_alrea(ivari) = 1
                 end if
              end if
              !
              ! At a given time
              !
              if( postp(1) % pos_alrea(ivari) == 0 ) then  
                 do itime = 1,nvart
                    if(   abs( postp(1) % pos_times(itime,ivari)-cutim) < (0.5_rp*dtime) .and. &
                         &     postp(1) % pos_times(itime,ivari)        > 0.0_rp) then
                       ivara = ivari
                       if( itask == 11 ) postp(1) % pos_alrea(ivari) = 1
                    end if
                 end do
              end if
              !
              ! At a given time period
              !
              if( ittim == 1 ) postp(1) % pos_times(1,ivari) = postp(1) % pos_perio(ivari)
              if( cutim >= postp(1) % pos_tinit .and. postp(1) % pos_alrea(ivari) == 0 ) then  
                 if(    abs(postp(1) % pos_times(1,ivari)-cutim) < (0.6_rp*dtime).and.&
                      &     postp(1) % pos_perio(ivari)          > 0.0_rp) then
                    ivara = ivari
                    if( itask == 11 ) then
                       postp(1) % pos_alrea(ivari)   = 1 
                       postp(1) % pos_times(1,ivari) = postp(1) % pos_times(1,ivari) + postp(1) % pos_perio(ivari) 
                    end if
                 end if
              end if

           end if

        else if( ittyp == ITASK_ENDRUN ) then
           !
           ! End of the run
           !        
           if( postp(1) % npp_stepi(ivari) /= 0 .and. postp(1) % pos_alrea(ivari) == 0 .and. &
                ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
              ivara = ivari   
           end if

        end if
     end if

     if( itask == 11 ) then
        if( ivara == ivari .and. postp(1) % kfl_oonce(ivari) == 1 ) then
           postp(1) % pos_alrea(ivari) = 2
        end if
     end if

     kfl_ivari(1) = ivara

     !-------------------------------------------------------------------
     !
     ! Voxel: Initial solution, end of time step, end of the run
     !
     !-------------------------------------------------------------------

     ivara = 0

     if( postp(1) % vox_alrea(ivari) /= 2 ) then
        if( ittyp == ITASK_INITIA .and. kfl_rstar == 0 ) then
           !
           ! Initial solution
           !
           if( postp(1) % npp_iniso == 1 .and. postp(1) % vox_stepi(ivari) > 0 ) then  
              ivara = ivari
              if( itask == 11 ) postp(1) % vox_alrea(ivari) = 1 ! To avoid repostprocess when no time step is required
           end if

        else if( ittyp == ITASK_ENDTIM ) then
           !
           ! End of a time step
           !
           if( cutim >= postp(1) % pos_tinit ) then
              if( itask == 11 ) postp(1) % vox_alrea(ivari) = 0
              !
              ! At a given time step
              !           
              if( &
                   ittim >= postp(1) % npp_inits .and. &
                   postp(1) % vox_alrea(ivari) == 0 .and. &
                   postp(1) % vox_stepi(ivari) > 0 ) then     
                 if( mod(ittim, postp(1) % vox_stepi(ivari)) == 0 ) then
                    ivara = ivari
                    if( itask == 11 ) postp(1) % vox_alrea(ivari) = 1
                 end if
              end if
              !
              ! At a given time
              !
              if( postp(1) % vox_alrea(ivari) == 0 ) then  
                 do itime = 1,nvart
                    if(   abs( postp(1) % vox_times(itime,ivari)-cutim) < (0.5_rp*dtime) .and. &
                         &     postp(1) % vox_times(itime,ivari)        > 0.0_rp) then
                       ivara = ivari
                       if( itask == 11 ) postp(1) % vox_alrea(ivari) = 1
                    end if
                 end do
              end if
              !
              ! At a given time period
              !
              if( ittim == 1 ) postp(1) % vox_times(1,ivari) = postp(1) % vox_perio(ivari)
              if( cutim >= postp(1) % pos_tinit .and. postp(1) % vox_alrea(ivari) == 0 ) then  
                 if(    abs(postp(1) % vox_times(1,ivari)-cutim) < (0.6_rp*dtime).and.&
                      &     postp(1) % vox_perio(ivari)          > 0.0_rp) then
                    ivara = ivari
                    if( itask == 11 ) then
                       postp(1) % vox_alrea(ivari)   = 1 
                       postp(1) % vox_times(1,ivari) = postp(1) % vox_times(1,ivari) + postp(1) % vox_perio(ivari) 
                    end if
                 end if
              end if

           end if

        else if( ittyp == ITASK_ENDRUN ) then
           !
           ! End of the run
           !        
           if( postp(1) % vox_stepi(ivari) /= 0 .and. postp(1) % vox_alrea(ivari) == 0 .and. &
                ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
              ivara = ivari   
           end if

        end if
     end if

     if( itask == 11 ) then
        if( ivara == ivari .and. postp(1) % kfl_oonce(ivari) == 1 ) then
           postp(1) % vox_alrea(ivari) = 2
        end if
     end if

     kfl_ivari(2) = ivara

     ivara = maxval(kfl_ivari)

  else if( itask == 21_ip ) then

     !-------------------------------------------------------------------
     !
     ! Element set: Parallel 
     !
     !-------------------------------------------------------------------

     if( IPARALL .and. maxval(postp(1) % npp_setse) > 0 ) then

        nullify(per_veset)
        allocate(per_veset(postp(1) % per_nvaes+1,neset))
        do nvaes = 1,postp(1) % per_nvaes
           ivars = postp(1) % per_setse(nvaes)
           per_veset(nvaes,:) = postp(1) % veset(ivars,:) 
        end do
        per_veset(postp(1) % per_nvaes+1,:) = postp(1) % veset(postp(1) % nvaes+1,:)
        
        call PAR_SUM(postp(1) % per_nvaes+1_ip,neset,per_veset,'IN MY CODE') 
        
        do nvaes = 1,postp(1) % per_nvaes
           ivars = postp(1) % per_setse(nvaes)
           postp(1) % veset(ivars,:) = per_veset(nvaes,:) 
        end do
        postp(1) % veset(postp(1) % nvaes+1,:) = per_veset(postp(1) % per_nvaes+1,:) 
        deallocate(per_veset)
        
     end if

  else if( itask == 22_ip ) then

     !-------------------------------------------------------------------
     !
     ! Boundary set: Parallel 
     !
     !-------------------------------------------------------------------

     if( IPARALL .and. maxval(postp(1) % npp_setsb) > 0 ) then

        nullify(per_vbset)
        allocate(per_vbset(postp(1) % per_nvabs+1,nbset))
        do nvabs = 1,postp(1) % per_nvabs
           ivars = postp(1) % per_setsb(nvabs)
           per_vbset(nvabs,:) = postp(1) % vbset(ivars,:) 
        end do
        per_vbset(postp(1) % per_nvabs+1,:) = postp(1) % vbset(postp(1) % nvabs+1,:)
        
        call PAR_SUM(postp(1) % per_nvabs+1_ip,nbset,per_vbset,'IN MY CODE')
        
        do nvabs = 1,postp(1) % per_nvabs
           ivars = postp(1) % per_setsb(nvabs)
           postp(1) % vbset(ivars,:) = per_vbset(nvabs,:) 
        end do
        postp(1) % vbset(postp(1) % nvabs+1,:) = per_vbset(postp(1) % per_nvabs+1,:) 
        deallocate(per_vbset)
        
     end if

  else if( itask == 23_ip ) then

     !-------------------------------------------------------------------
     !
     ! Node set: Parallel 
     !
     !-------------------------------------------------------------------

     if( IPARALL .and. maxval(postp(1) % npp_setsn) > 0 ) then
        nparr =  postp(1) % nvans*nnset
        vnset => postp(1) % vnset
        call Parall(15_ip)
     end if

  else if( itask == 24_ip ) then

     !-------------------------------------------------------------------
     !
     ! Witness: Parallel 
     !
     !-------------------------------------------------------------------

     if( IPARALL .and. nwitn > 0 .and. maxval(postp(1) % npp_witne) > 0 ) then
        call Parall(39_ip)
     end if

  else if( itask == 25_ip ) then

     !-------------------------------------------------------------------
     !
     ! Check if variable IVARA is postprocessed
     !
     !-------------------------------------------------------------------

     if(  postp(1) % npp_stepi(ivara)           /= 0      .or. &
          maxval(postp(1) % pos_times(:,ivara)) >  0.0_rp .or. &
          postp(1) % pos_perio(ivara)           /= 0.0_rp ) then
        continue
     else
        ivara = 0
     end if

  else if( itask == 26_ip ) then

     !-------------------------------------------------------------------
     !
     ! Cancel postprocess of variable IVARA
     !
     !-------------------------------------------------------------------

     postp(1) % npp_stepi(ivara)   = 0
     postp(1) % pos_times(:,ivara) = 0.0_rp
     postp(1) % pos_perio(ivara)   = 0.0_rp

  else if( itask == 27_ip ) then

     !-------------------------------------------------------------------
     !
     ! Redefine sets to be postprocessed... just in case modules have
     ! redefined NPP_SETSB outside of this subroutine
     !
     !-------------------------------------------------------------------

     nvaes = 0                
     do ivars = 1,nvars
        if( postp(1) % npp_setse(ivars) == 1 ) then
           nvaes                       = nvaes + 1
           postp(1) % per_setse(nvaes) = ivars    
        end if
     end do
     postp(1) % per_nvaes = nvaes
     
     nvabs = 0                
     do ivars = 1,nvars
        if( postp(1) % npp_setsb(ivars) == 1 ) then
           nvabs                       = nvabs + 1 
           postp(1) % per_setsb(nvabs) = ivars     
        end if
     end do
     postp(1) % per_nvabs = nvabs

  end if

end subroutine posdef
