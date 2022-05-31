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
  use mod_communications,     only : PAR_SUM
  use mod_ecoute,             only : ecoute
  use mod_memory,             only : memory_alloca
  use mod_witness,            only : WITNESS_INSTANTANEAOUS
  use mod_witness,            only : WITNESS_AVERAGE
  use mod_witness,            only : WITNESS_ACCUMULATE
  use mod_witness,            only : witness_geometry_averaging_ini
  use mod_witness,            only : witness_geometry_averaging_end
  use mod_output_postprocess, only : output_postprocess_read
  use mod_output_postprocess, only : output_postprocess_parall
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess_now
  implicit none  
  integer(ip), intent(in)    :: itask
  integer(ip), intent(inout) :: ivara
  integer(ip)                :: ivari,imodu,ivarp,ivarw,ivart,ivars,itime
  integer(ip)                :: ipost,nvabs,nvaes,nvans,iwitg,iposi,ii,imesh
  integer(4)                 :: istat
  character(5)               :: CNULL='NULL '
  real(rp),    pointer       :: per_vbset(:,:)
  real(rp),    pointer       :: per_veset(:,:)
  
  if( itask == 0_ip ) then

     call runend('PODEF: DEPRECATED')
     
  else if( itask == 1_ip ) then

     !-------------------------------------------------------------------
     !
     ! Used for Parall service
     !
     !-------------------------------------------------------------------

     call output_postprocess_parall()

  else if( itask == 2_ip ) then

     !-------------------------------------------------------------------
     !
     ! Read postprocess
     !
     !-------------------------------------------------------------------

     call output_postprocess_read()

  else if( itask == 3_ip ) then

     !-------------------------------------------------------------------
     !
     ! Allocate memory for sets and witness points
     !
     !-------------------------------------------------------------------

     call runend('POSDEF: DEPRECATED')

  else if( itask == 4_ip .and. INOTSLAVE ) then

     call runend('POSDEF: DEPRECATED')

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

     if( .not. output_postprocess_check_variable_postprocess_now(ivara) ) ivara = 0
     
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
        call par_comset(1_ip)   
     end if

  else if( itask == 24_ip ) then

     call runend('POSDEF: NO MORE...')

  else if( itask == 25_ip ) then

     !-------------------------------------------------------------------
     !
     ! Check if variable IVARA is postprocessed
     !
     !-------------------------------------------------------------------

     do imesh = 0,nwith
        if(  postp(1) % npp_stepi(ivara,imesh)           /= 0      .or. &
             maxval(postp(1) % pos_times(:,ivara,imesh)) >  0.0_rp .or. &
             postp(1) % pos_perio(ivara,imesh)           /= 0.0_rp ) then
           continue
        else
           ivara = 0
           return
        end if
     end do
     
  else if( itask == 26_ip ) then

     !-------------------------------------------------------------------
     !
     ! Cancel postprocess of variable IVARA
     !
     !-------------------------------------------------------------------

     postp(1) % npp_stepi(ivara,:)   = 0
     postp(1) % pos_times(:,ivara,:) = 0.0_rp
     postp(1) % pos_perio(ivara,:)   = 0.0_rp

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
