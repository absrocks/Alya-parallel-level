!----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_plugin.f90
!> @author  J.C. Cajas
!> @date    17/04/2014
!> @brief   Receive force from Nastin and send displacement to Alefor
!> @details Receive force from Nastin and send displacement to Alefor
!> @        usind the coupling structures and functions.
!> @} 
!----------------------------------------------------------------------
subroutine sld_plugin(icoup)
  !
  ! Obligatory variables 
  !
  use def_coupli,         only :  coupling_type
  use def_domain,         only :  npoin,ndime, nmate
  use def_master,         only :  solve_sol, mem_modul,modul
  use def_master,         only :  kfl_eccty, kfl_cellmod, kfl_exm_max_nmaterials
  use def_kintyp,         only :  ip,rp,lg
  use def_master,         only :  current_code,modul,mem_modul
  use def_master,         only :  INOTMASTER
  use mod_parall,         only :  PAR_MY_CODE_RANK
  use mod_couplings,      only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_couplings,      only :  COU_SET_FIXITY_ON_TARGET
  use mod_communications, only :  PAR_SUM, PAR_BARRIER
  use mod_memory,         only :  memory_deallo
  use mod_memory,         only :  memory_alloca,memory_alloca_min
  use mod_matrix,         only :  matrix_initialize
  use mod_exm_cou,        only :  mod_exm_cou_initexchange
  use mod_exm_cou,        only :  mod_exm_cou_physics_initialised
  !
  ! Possible variables 
  !
  use def_master,        only :  displ
  use def_master,        only :  ID_NASTIN
  use def_master,        only :  gdepo
  use def_solidz,        only :  bvess_sld
  use def_solidz,        only :  kfl_gdepo
  use def_solidz,        only :  kfl_fixno_sld
  use def_solidz,        only :  kfl_minco_sld
  use mod_parall,        only :  PAR_GLOBAL_TO_LOCAL_NODE
  use mod_exm_sld_eccoupling, only: EXMSLD_EMEC_LAND, EXMSLD_EMEC_LAND_BIDIR
  use mod_exm_sld_eccoupling, only: calcium_ecc, state_ecc, troponin_ecc
  use mod_exm_sld_eccoupling, only: has_land

  implicit none 
  real(rp),    pointer    :: xvalu(:,:)
  real(rp),    pointer    :: svalu(:,:)
  real(rp),    pointer    :: svalu_no_dealloc(:,:)
  real(rp)                :: foref  ! for the coupling with nastin

  integer(ip)             :: i
  integer(ip)             :: idime,jdime, imate
  integer(ip)             :: ipoin,kpoin
  integer(ip), intent(in) :: icoup
  character(5)            :: variable
  integer(ip), save       :: ipass=0

  nullify(xvalu) 
  nullify(svalu)

  variable = coupling_type(icoup) % variable

  if( variable == 'DISPL' ) then   
     !
     ! Displacement
     !
     if( ipass == 0 ) then
        ipass = 1
        call COU_SET_FIXITY_ON_TARGET('DISPL',modul,kfl_fixno_sld)
     end if     
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,bvess_sld,displ,kfl_fixno_sld)

  else if( variable == 'CALCI' ) then   
     !
     ! Calcium concentration (vconc(1,ipoin,1))
     !
     if(ipass.eq.0_ip) then
        ipass = 1
        call mod_exm_cou_physics_initialised()
     endif

     if (INOTMASTER) then
        call memory_alloca(mem_modul(1:2,modul),'XVALU','sld_plugin',xvalu,1_ip,npoin)
        call memory_alloca(mem_modul(1:2,modul),'SVALU','sld_plugin',svalu,1_ip,npoin)          
     else
        call memory_alloca(mem_modul(1:2,modul),'XVALU','sld_plugin',xvalu,1_ip,1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SVALU','sld_plugin',svalu,1_ip,1_ip)          
     end if

     ! RECEIVE CALCIUM
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,calcium_ecc,svalu)

     if(has_land) then
      ! COMPUTE TROPONIN_PREV BEFORE RECEIVING NEW VALUE
      troponin_ecc(:,2) = troponin_ecc(:,1)

      ! RECEIVE TROPONIN
      call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,xvalu,svalu)
      troponin_ecc(:,1)=xvalu(1,:)


      ! RECEIVE STATELAND
      call runend('SLD_PLUGIN: LAND COUPLING NOT POSSIBLE IN VOLUMETRIC COUPLING. READ CODE')
      ! Sending stateland from exmedi to solidz is stupid for 2 reasons:
      !     1) It's element-wise (ielem vector)
      !     2) It's done nothing with it in exmedi, just a gather.
     endif


  else if( variable == 'RESID' .or. variable == 'MOMEN' .or. variable == 'FWALL' ) then
     !
     ! Residual
     !

     if(kfl_minco_sld.eq.0_ip) kfl_minco_sld=icoup
     if (icoup .eq. kfl_minco_sld) call matrix_initialize(solve_sol(1) % bvnat)

     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,solve_sol(1) % bvnat,solve_sol(1) % reaction)
     if( coupling_type(icoup) % module_source == ID_NASTIN )then
        if( kfl_gdepo == 0_ip ) call runend('ROTATION ON, MUST BE ADDED IN NUMERICAL TREATMENT .SLD.DAT')
        !
        ! Push forward for the force coming from nastin
        !
        if( INOTMASTER ) then

           do ipoin = 1,npoin

              do idime = 1,ndime

                 foref = 0.0_rp

                 if( kfl_fixno_sld(idime,ipoin) /= 1 ) then

                    do jdime = 1,ndime

                       foref = foref + gdepo(idime,jdime,ipoin) * solve_sol(1) % bvnat(jdime,ipoin) 

                    end do

                    solve_sol(1) % bvnat(idime,ipoin) = foref

                 end if

              end do

           end do

        end if

     end if

  else if( variable == 'ALEFO' ) then
     !
     ! Coupling with alefor
     !
     if ( INOTMASTER ) then

        allocate(svalu(ndime, npoin))

        svalu = 0.0_rp
        !
        ! If there are subcycles the displacement is calculated with the saved values
        !
        if( coupling_type(icoup) % frequ_send > 1_ip .or. coupling_type(icoup) % frequ_recv > 1_ip )then

           do kpoin = 1_ip, coupling_type(icoup) % geome % npoin_source

              ipoin = coupling_type(icoup) % geome % lpoin_source(kpoin)

              do idime = 1_ip, ndime
                 svalu(idime,ipoin) = displ(idime,ipoin,1_ip) - coupling_type(icoup) % values_frequ(idime,kpoin,1_ip)
                 coupling_type(icoup) % values_frequ(idime,kpoin,2_ip) = displ(idime,ipoin,1_ip)
              end do

           end do

        else
           !
           ! If the exchanges are on each time step, the displacement in calculated with displ only
           !
           do kpoin = 1_ip, coupling_type(icoup) % geome % npoin_source

              ipoin = coupling_type(icoup) % geome % lpoin_source(kpoin)
              do idime = 1_ip, ndime
                 svalu(idime,ipoin) = displ(idime,ipoin,1_ip) - displ(idime,ipoin,3_ip)
              end do

           end do

        end if

     else

        allocate(svalu(1_ip,1_ip))

     end if

     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,bvess_sld,svalu)

  end if

  if( associated(xvalu) ) call memory_deallo(mem_modul(1:2,modul),'XVALU','sld_plugin',xvalu)
  if( associated(svalu) ) call memory_deallo(mem_modul(1:2,modul),'SVALU','sld_plugin',svalu)

end subroutine sld_plugin
