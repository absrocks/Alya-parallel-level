!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    mod_nsi_arrays.f90
!> @author  houzeaux
!> @date    2019-11-16
!> @brief   Nastin arrays
!> @details Nastin arrays
!-----------------------------------------------------------------------

module mod_nsi_arrays

  use def_master
  use def_domain 
  use def_nastin
  use def_kermod
  use mod_output_postprocess,  only : output_postprocess_check_variable_postprocess
  use mod_arrays,              only : arrays
  use mod_memory,              only : memory_alloca

  
  implicit none

  private

  public :: nsi_arrays
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-11-16
  !> @brief   Nastin arrays
  !> @details Do what you have to do with nastin arrays
  !> 
  !-----------------------------------------------------------------------
  
  subroutine nsi_arrays(wtask)

    character(len=*), intent(in) :: wtask
    integer(ip)                  :: ncomp_loc,ielem,igaus
    integer(ip)                  :: pgaus,ncsgs,pelty,ipoin

    if (NSI_FRACTIONAL_STEP.and.kfl_stabi_nsi == NSI_ASGS) then
       ncomp_loc = ncomp_nsi+1
    else
       ncomp_loc = ncomp_nsi
    end if
    !
    ! Unknowns
    !
    call arrays(1_ip,trim(wtask),veloc,ndime,npoin,ncomp_loc)
    call arrays(2_ip,trim(wtask),press,npoin,ncomp_nsi)
    !
    ! DENSI: Compressible regime
    !
    if( kfl_regim_nsi == 1 ) then
       call arrays(7_ip,trim(wtask),densi,npoin,1_ip)
    else if( kfl_regim_nsi == 2 ) then
       call arrays(7_ip,trim(wtask),densi,npoin,ncomp_nsi)
    end if
    !
    ! VEOLD_NSI: Allocate memory for old (last iteration) velocity
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='RESID') ) then
       kfl_resid_nsi=1
       call arrays(111_ip,trim(wtask),veold_nsi,ndime,npoin)        
    else
       kfl_resid_nsi=0
    end if
    !
    ! VEPRO_NSI, PRPRO_NSI, GRPRO_NSI: Projections for orthogonal SGS
    ! 
    if( kfl_stabi_nsi > 0 ) then
       call arrays(23_ip,trim(wtask),vepro_nsi,ndime,npoin)        
       call arrays(24_ip,trim(wtask),prpro_nsi,npoin)  
       if( kfl_stabi_nsi == 2 ) then
          call arrays(27_ip,trim(wtask),grpro_nsi,ndime,npoin)
       end if
    end if
    !
    ! Immersed boundary method
    !
    if(      kfl_immer_nsi == 1 ) then
       call arrays(85_ip,trim(wtask),lagra_nsi,ndime,npoin,1_ip)
    else if( kfl_immer_nsi == 2 ) then
       call arrays(112_ip,trim(wtask),tauib_nsi,ndime,ndime,npoin)
    end if
    !
    ! VESGS: Subgrid scale velocity
    !
    if( kfl_sgsco_nsi == 1 .or. kfl_sgsti_nsi == 1 ) then
       call arrays(5_ip,trim(wtask),vesgs,nelem)
       if( trim(wtask) == 'ALLOCATE' ) then
          ncsgs = min(2_ip,2_ip*kfl_sgsti_nsi+kfl_sgsco_nsi)
          do ielem = 1,nelem
             pelty = abs(ltype(ielem))
             pgaus = ngaus(pelty)
             call memory_alloca(mem_modul(1:2,modul),'VESGS % A','nsi_memall',vesgs(ielem)%a,ndime,pgaus,ncsgs)
          end do
       end if
    end if
    !
    ! DUNKN_NSI, DUNKP_NSI: Delta velocity and pressure for Aitken relaxation strategy
    !
    if(kfl_relax_nsi==2) then
       call arrays(113_ip,trim(wtask),dunkn_nsi,ndime,npoin)
    end if
    if(kfl_relap_nsi==2) then
       call arrays(114_ip,trim(wtask),dunkp_nsi,npoin)
    end if
    !
    ! MASS_RHO_NSI
    !
    if( kfl_predi_nsi == 7 .or. kfl_predi_nsi == 8 .or. kfl_predi_nsi == 9 ) then
       call arrays(108_ip,trim(wtask),mass_rho_nsi,npoin,ncomp_nsi)
    end if
    !
    ! arrays for no slip wall law - later they can be used for other cases too
    ! I have intialized tehm to 0 here but ther better option might be possible
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='VAFOR') .or. kfl_noslw_ker /= 0 )  then
       call arrays(101_ip,trim(wtask),vafor_nsi,ndime,npoin)       
    end if
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVAF') .or. kfl_noslw_ker /= 0 )  then
       call arrays(102_ip,trim(wtask),avvaf_nsi,ndime,npoin)       
    end if
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVNTR') .or. kfl_noslw_ker /= 0 )  then
       call arrays(104_ip,trim(wtask),avntr_nsi,ndime,npoin)       
    end if
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVGTR') .or. kfl_noslw_ker /= 0 )  then
       call arrays(105_ip,trim(wtask),avgtr_nsi,ndime,npoin)       
    end if
    !
    ! Coupling with SOLIDZ
    !
    if( coupling('SOLIDZ','NASTIN') >= 1 .or. coupling('NASTIN','IMMBOU') >= 1 ) then
       call arrays(49_ip,trim(wtask),forcf,ndime,npoin)       
    end if
    !
    ! Traction on boundary nodes for alternative AVTAN postprocess
    !
    if(    kfl_noslw_ker /= 0                                                   .or. &
         & output_postprocess_check_variable_postprocess(VARIABLE_NAME='NOTRA') .or. &
         & output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTAN') .or. &
         & output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENTAN') .or. &
         & output_postprocess_check_variable_postprocess(VARIABLE_NAME='TANGE') .or. &
         & output_postprocess_check_variable_postprocess(VARIABLE_NAME='FTANG')      )  then
       call arrays(103_ip,trim(wtask),notra_nsi,ndime,npoin)       
    end if
    !
    ! AVVEL_NSI: average velocity (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVEL') ) then
       call arrays(21_ip,trim(wtask),avvel_nsi,ndime,npoin)       
    end if
    !
    ! AVPRE_NSI: average pressure (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVPRE') ) then
       call arrays(28_ip,trim(wtask),avpre_nsi,npoin)       
    end if
    !
    ! AVVE2_NSI: average velocity**2 (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVE2') ) then
       call arrays(57_ip,trim(wtask),avve2_nsi,ndime,npoin)       
    end if
    !
    ! AVVXY_NSI: average vx*vy (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVVXY') ) then
       call arrays(58_ip,trim(wtask),avvxy_nsi,ndime,npoin)       
    end if
    !
    ! AVPR2_NSI: average pressure (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVPR2') ) then
       call arrays(59_ip,trim(wtask),avpr2_nsi,npoin)       
    end if
    !
    ! AVTAN_NSI: average TANGE (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVTAN') ) then
       call arrays(46_ip,trim(wtask),avtan_nsi,ndime,npoin)       
    end if
    !
    ! AVMUT_NSI: Average turbulent viscosity (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMUT') ) then
       call arrays(75_ip,trim(wtask),avmut_nsi,npoin)       
    end if
    !
    ! AVSTX_NSI: Average stress mu_t * grad(u)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVSTX') ) then
       call arrays(76_ip,trim(wtask),avstx_nsi,ndime,npoin)       
    end if
    !
    ! AVSTY_NSI: Average stress mu_t * grad(v)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVSTY') ) then
       call arrays(77_ip,trim(wtask),avsty_nsi,ndime,npoin)       
    end if
    !
    ! AVSTZ_NSI: Average stress mu_t * grad(w)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVSTZ') ) then
       call arrays(78_ip,trim(wtask),avstz_nsi,ndime,npoin)       
    end if
    !
    ! AVMOS_NSI: average momentum source from spray
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='AVMOS') ) then
       call arrays(122_ip,trim(wtask),avmos_nsi,ndime,npoin)       
    end if
    !
    ! ENVEL_NSI: ensemble velocity (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENVEL') ) then
       call arrays(86_ip,trim(wtask),envel_nsi,ndime,npoin)       
    end if
    !
    ! ENPRE_NSI: ensemble pressure (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENPRE') ) then
       call arrays(87_ip,trim(wtask),enpre_nsi,npoin)       
    end if
    !
    ! ENVE2_NSI: ensemble velocity**2 (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENVE2') ) then
       call arrays(88_ip,trim(wtask),enve2_nsi,ndime,npoin)       
    end if
    !
    ! ENVXY_NSI: ensemble vx*vy (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENVXY') ) then
       call arrays(89_ip,trim(wtask),envxy_nsi,ndime,npoin)       
    end if
    !
    ! ENPRE_NSI: ensemble pressure (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENPR2') ) then
       call arrays(90_ip,trim(wtask),enpr2_nsi,npoin)       
    end if
    !
    ! ENTAN_NSI: ensemble TANGE (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENTAN') ) then
       call arrays(91_ip,trim(wtask),entan_nsi,ndime,npoin)       
    end if
    !
    ! ENMUT_NSI: Ensemble turbulent viscosity (postprocess)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENMUT') ) then
       call arrays(92_ip,trim(wtask),enmut_nsi,npoin)       
    end if
    !
    ! ENSTX_NSI: Ensemble stress mu_t * grad(u)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENSTX') ) then
       call arrays(93_ip,trim(wtask),enstx_nsi,ndime,npoin)       
    end if
    !
    ! ENSTY_NSI: Ensemble stress mu_t * grad(v)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENSTY') ) then
       call arrays(94_ip,trim(wtask),ensty_nsi,ndime,npoin)       
    end if
    !
    ! ENSTZ_NSI: Ensemble stress mu_t * grad(w)
    !
    if( output_postprocess_check_variable_postprocess(VARIABLE_NAME='ENSTZ') ) then
       call arrays(95_ip,trim(wtask),enstz_nsi,ndime,npoin)       
    end if
    !
    ! Bubble
    !
    if( kfl_bubbl_nsi /= 0 ) then
       call arrays(115_ip,trim(wtask),bubble_nsi    ,nelem)            
       call arrays(116_ip,trim(wtask),bubble_aqq_nsi,nelem)            
       call arrays(117_ip,trim(wtask),bubble_aqu_nsi,mnode*ndime,nelem)
       call arrays(118_ip,trim(wtask),bubble_aqp_nsi,mnode,nelem)      
       call arrays(119_ip,trim(wtask),bubble_bq_nsi ,nelem)                    
    end if
    !
    ! RANS/LES two-layer model
    !
    if( kfl_twola_ker > 0 ) then
       call arrays(97_ip,trim(wtask),btrac_nsi,ndime,npoin)                    
       call arrays(98_ip,trim(wtask),tracr_nsi,ndime,npoin)  
       call arrays(99_ip,trim(wtask),tluav_nsi,ndime,npoin) 
    end if

  end subroutine nsi_arrays
   
end module mod_nsi_arrays
!> @}
