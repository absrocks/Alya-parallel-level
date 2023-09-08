!-----------------------------------------------------------------------
!> @addtogroup Postprocess
!> @{
!> @file    mod_postprocess.f90
!> @author  houzeaux
!> @date    2018-11-08
!> @brief   Postprocess sruff
!> @details Output and postprocess stuffs
!-----------------------------------------------------------------------

module mod_output_postprocess

  use def_parame
  use def_domain
  use def_master 
  use def_inpout
  use def_postpr
  use mod_memory,         only : memory_alloca
  use mod_communications, only : PAR_SUM
  use mod_ecoute,         only : ecoute

  implicit none

  character(5) :: CNULL='NULL '

  private

  public :: output_postprocess_initialization                 !  0
  public :: output_postprocess_parall                         !  1
  public :: output_postprocess_read                           !  2
  public :: output_postprocess_allocate_witness               !  3
  public :: output_postprocess_output_witness                 !  4
  public :: output_postprocess_element_sets_parall            ! 21
  public :: output_postprocess_boundary_sets_parall           ! 22
  public :: output_postprocess_node_sets_parall               ! 23
  public :: output_postprocess_witness_parall                 ! 24
  public :: output_postprocess_cancel_variable_postprocess    ! 26
  public :: output_postprocess_sets_ordering                  ! 27
  public :: output_postprocess_check_witness_output           !  5
  public :: output_postprocess_check_variable_postprocess_now ! 11 and 111
  public :: output_postprocess_check_variable_postprocess     ! 25

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Initialization
  !> @details Allocation and inititalization of POSTP structure
  !>          POSDEF(0)
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_initialization()

    integer(ip) :: imodu,ivarp,ivars,ivarw,ivart,ivari

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

  end subroutine output_postprocess_initialization

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Parallelization
  !> @details Exchange POSTP structure 
  !>          POSDEF(1)
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_parall()

    integer(ip) :: ivarp,ivars,ivarw,ivart,ivari
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

  end subroutine output_postprocess_parall

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Read
  !> @details Read from data file the postprocess and output options
  !>          POSDEF(2)
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_read()

    integer(ip)                :: ivari,imodu,ivarp,ivarw,ivart,ivars,itime
    integer(ip)                :: ipost,nvabs,nvaes,nvans
    integer(4)                 :: istat

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

  end subroutine output_postprocess_read

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Allocate memory for witness
  !> @details Allocate memory for witness points
  !>          POSDEF(3)
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_allocate_witness()

    if( postp(1) % ipass == 0 ) then

       postp(1) % ipass = 1

       if(maxval(postp(1) % npp_setse)>0) then
          call memory_alloca(mem_modul(1:2,modul),'VESET','mod_output_postprocess',postp(1) % veset,postp(1) % nvaes+1,neset)
       end if
       if(maxval(postp(1) % npp_setsb)>0) then
          call memory_alloca(mem_modul(1:2,modul),'VBSET','mod_output_postprocess',postp(1) % vbset,postp(1) % nvabs+1,nbset)
       end if
       if(maxval(postp(1) % npp_setsn)>0) then
          call memory_alloca(mem_modul(1:2,modul),'VNSET','mod_output_postprocess',postp(1) % vnset,postp(1) % nvans,nnset)
       end if
       if(maxval(postp(1) % npp_witne)>0) then
          call memory_alloca(mem_modul(1:2,modul),'WITNE','mod_output_postprocess',postp(1) % witne,postp(1) % nvawi,nwitn)
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

  end subroutine output_postprocess_allocate_witness

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Output witness
  !> @details Output witness points values
  !>          POSDEF(4)
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_output_witness()

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

  end subroutine output_postprocess_output_witness
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Element set parallelization
  !> @details Element set parallelization
  !>          POSDEF(21)
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_element_sets_parall()

    integer(ip)         :: nvaes,ivars 
    real(rp),   pointer :: per_veset(:,:)
    
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

  end subroutine output_postprocess_element_sets_parall

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Element set parallelization
  !> @details Element set parallelization
  !>          POSDEF(22)
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_boundary_sets_parall()

   integer(ip)         :: nvabs,ivars 
   real(rp),   pointer :: per_vbset(:,:)

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

  end subroutine output_postprocess_boundary_sets_parall

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Element set parallelization
  !> @details Element set parallelization
  !>          POSDEF(23)
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_node_sets_parall()

    if( IPARALL .and. maxval(postp(1) % npp_setsn) > 0 ) then
       nparr =  postp(1) % nvans*nnset
       vnset => postp(1) % vnset
       call Parall(15_ip)
    end if

  end subroutine output_postprocess_node_sets_parall


  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Witness point parallelization
  !> @details Witness point parallelization
  !>          POSDEF(24)
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_witness_parall()

    if( IPARALL .and. nwitn > 0 .and. maxval(postp(1) % npp_witne) > 0 ) then
       call Parall(39_ip)
    end if
    
  end subroutine output_postprocess_witness_parall
    
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Check if a variable is postprocessed
  !> @details Check if a variable is postprocessed at one moment
  !>          POSDEF(26)
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_cancel_variable_postprocess(ivara)

    integer(ip), intent(in) :: ivara
    
    postp(1) % npp_stepi(ivara)   = 0
    postp(1) % pos_times(:,ivara) = 0.0_rp
    postp(1) % pos_perio(ivara)   = 0.0_rp

  end subroutine output_postprocess_cancel_variable_postprocess

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Sets reordering
  !> @details Redefine sets to be postprocessed... just in case modules have
  !>          redefined NPP_SETSB outside of this subroutine
  !>          POSDEF(27)
  !> 
  !-----------------------------------------------------------------------

  subroutine output_postprocess_sets_ordering()

    integer(ip) :: ivars,nvabs,nvaes

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

  end subroutine output_postprocess_sets_ordering

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Check if a witness point is output
  !> @details Check if a witness point is output
  !>          POSDEF(5)
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function output_postprocess_check_witness_output()
    
    output_postprocess_check_witness_output = .false.
    if( maxval(postp(1) % npp_witne) > 0 .and. ( kfl_rstar /= 2 .or. itti2 /= 0 ) ) then
       if(( mod(ittim, postp(1) % npp_stepw) == 0 ).or. (ittim .eq. mitim ) .or. (cutim >= timef)) then
          output_postprocess_check_witness_output = .true.
       end if
    end if
    
  end function output_postprocess_check_witness_output
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Check if a variable should be postprocessed
  !> @details Check if a variable should be postprocessed
  !>          POSDEF(11), POSDEF(111)
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function output_postprocess_check_variable_postprocess_now(FILTER)

    logical(lg), intent(in), optional :: FILTER
    integer(ip)                       :: ivari,imodu,ivarp,ivarw,ivart,ivars,itime
    integer(ip)                       :: ipost,nvabs,nvaes,nvans,ivara
    logical(lg)                       :: if_filter

    if( present(FILTER) ) then
       if_filter = FILTER
    else
       if_filter = .false.
    end if

    ivari     = ivara
    ivara     = 0
    kfl_ivari = 0

    if( postp(1) % pos_alrea(ivari) /= 2 ) then
       
       if( ittyp == ITASK_INITIA .and. kfl_rstar < 2 ) then  ! outputs if INITI or no restart
          !
          ! Initial solution
          !
          if( postp(1) % npp_iniso == 1 .and. postp(1) % npp_stepi(ivari) > 0 ) then  
             ivara = ivari
             if( .not. if_filter ) postp(1) % pos_alrea(ivari) = 1 ! To avoid repostprocess when no time step is required
          end if

       else if( ittyp == ITASK_ENDTIM ) then
          !
          ! End of a time step
          !
          if( cutim >= postp(1) % pos_tinit ) then
             if( .not. if_filter ) postp(1) % pos_alrea(ivari) = 0
             !
             ! At a given time step
             !           
             if( &
                  ittim >= postp(1) % npp_inits .and. &
                  postp(1) % pos_alrea(ivari) == 0 .and. &
                  postp(1) % npp_stepi(ivari) > 0 ) then     
                if( mod(ittim, postp(1) % npp_stepi(ivari)) == 0 ) then
                   ivara = ivari
                   if( .not. if_filter ) postp(1) % pos_alrea(ivari) = 1
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
                      if( .not. if_filter ) postp(1) % pos_alrea(ivari) = 1
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
                   if( .not. if_filter ) then
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

    if( .not. if_filter ) then
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
             if( .not. if_filter ) postp(1) % vox_alrea(ivari) = 1 ! To avoid repostprocess when no time step is required
          end if

       else if( ittyp == ITASK_ENDTIM ) then
          !
          ! End of a time step
          !
          if( cutim >= postp(1) % pos_tinit ) then
             if( .not. if_filter ) postp(1) % vox_alrea(ivari) = 0
             !
             ! At a given time step
             !           
             if( &
                  ittim >= postp(1) % npp_inits .and. &
                  postp(1) % vox_alrea(ivari) == 0 .and. &
                  postp(1) % vox_stepi(ivari) > 0 ) then     
                if( mod(ittim, postp(1) % vox_stepi(ivari)) == 0 ) then
                   ivara = ivari
                   if( .not. if_filter ) postp(1) % vox_alrea(ivari) = 1
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
                      if( .not. if_filter ) postp(1) % vox_alrea(ivari) = 1
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
                   if( .not. if_filter ) then
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

    if( .not. if_filter ) then
       if( ivara == ivari .and. postp(1) % kfl_oonce(ivari) == 1 ) then
          postp(1) % vox_alrea(ivari) = 2
       end if
    end if

    kfl_ivari(2) = ivara

    ivara = maxval(kfl_ivari)

    if( ivara == 0 ) then
       output_postprocess_check_variable_postprocess_now = .false.
    else
       output_postprocess_check_variable_postprocess_now = .true.
    end if
    
  end function output_postprocess_check_variable_postprocess_now
 
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-11-08
  !> @brief   Check if a variable is postprocessed
  !> @details Check if a variable is postprocessed at one moment
  !>          POSDEF(25)
  !> 
  !-----------------------------------------------------------------------

  logical(lg) function output_postprocess_check_variable_postprocess(ivara)

    integer(ip), intent(in) :: ivara
    
    output_postprocess_check_variable_postprocess = .true.
    
    if(  postp(1) % npp_stepi(ivara)           /= 0      .or. &
         maxval(postp(1) % pos_times(:,ivara)) >  0.0_rp .or. &
         postp(1) % pos_perio(ivara)           /= 0.0_rp ) then
       continue
    else
       output_postprocess_check_variable_postprocess = .false.
    end if

  end function output_postprocess_check_variable_postprocess
  
end module mod_output_postprocess
!> @}
