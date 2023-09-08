!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_iniunk.f90
!> @author  Mariano Vazquez
!> @author  Eva Casoni
!> @todo    To finish X-FEM
!> @author  Gerard Guillamet
!> @date    July, 2017-Added velocities
!> @todo    <GGU> Revise initial conditions for displacement and stress
!> @brief   Initial conditions for solidz module
!> @details Set up the initial conditions for displacement and velocity.
!>          If this is a restart, initial conditions are read from files.
!>          The different options without a restart are:
!>          \verbatim
!>          Displacement
!>          --------
!>          to fill....
!>
!>          Velocity
!>          --------
!>          KFL_INVEL_SLD = 0 ... Deactivated
!>                        = 1 ... Values from CODES fields
!>
!>          X-FEM
!>          --------
!>          KFL_XFEME_SLD = 0 ... Deactivated
!>                        = 1 ... Activated
!>
!>          Cohesive elements (Pandolfi)
!>          --------
!>          KFL_COH       = 0 ... Deactivated
!>                        = 1 ... ini nodes cohe activated
!>
!>          \endverbatim
!>          Values are stored in position:
!>          \verbatim
!>          VELOC_SLD(1:NDIME,1:NPOIN,NCOMP_SLD)
!>          \endverbatim
!> @}
!-----------------------------------------------------------------------

subroutine sld_iniunk()

  use def_kintyp,  only : ip,rp
  use def_master,  only : TIME_N,ITER_K,INOTMASTER,ITASK_INIUNK
  use def_master,  only : mem_modul,modul
  use def_master,  only : kfl_rstar,nfacg,displ
  use def_domain,  only : ndime,npoin,xfiel
  use def_domain,  only : kfl_elcoh
  use mod_cutele,  only : cutele
  use mod_memory,  only : memory_alloca
  use def_solidz,  only : veloc_sld,accel_sld
  use def_solidz,  only : kfl_invel_sld, kfl_indis_sld
  use def_solidz,  only : cockf_sld,lcrkf_sld,crtip_sld,cranx_sld,crapx_sld
  use def_solidz,  only : bvess_sld,nprev_sld
  use def_solidz,  only : kfl_xfeme_sld
  use def_solidz,  only : kfl_rigid_sld
  use mod_sld_rbo, only : sld_rbo_iniunk

  implicit none

  integer(ip)         :: ipoin,idime

  if ( kfl_rstar == 0 ) then

     !-------------------------------------------------------------------
     !
     ! Normal run
     !
     !-------------------------------------------------------------------

     if ( kfl_rigid_sld == 0_ip ) then

        !----------------------------------------------------------------
        !
        ! Deformable body
        !
        !----------------------------------------------------------------

        if ( INOTMASTER ) then

           if ( kfl_indis_sld(1) < 0 ) then
              !
              ! Displacement for Pre-stress
              !
              do ipoin = 1,npoin
                 do idime=1,ndime
                    If (kfl_indis_sld(1) < 0) then ! an initial displacement is given as a field
                       displ(idime,ipoin,nprev_sld)= xfiel(-kfl_indis_sld(1)) % a(idime,ipoin,1)
                    else if (kfl_indis_sld(1) == 0) then
                       ! <GGU> Revisar esto
                       displ(idime,ipoin,nprev_sld)=bvess_sld(idime,ipoin,1)
                    end if
                    veloc_sld(idime,ipoin,nprev_sld) = 0.0_rp
                    accel_sld(idime,ipoin,nprev_sld) = 0.0_rp
                 end do
              end do
           end if

           if ( kfl_invel_sld == 1_ip ) then
              !
              ! Velocity
              !
              veloc_sld(1:ndime, 1:npoin, TIME_N) = xfiel(1) % a(1:ndime, 1:npoin, ITER_K)
           end if

           if ( kfl_xfeme_sld > 0 ) then
              !
              ! X-FEM
              !
              ! Allocate memory
              call memory_alloca(mem_modul(1:2,modul),'LCRKF_SLD','sld_iniunk',lcrkf_sld,nfacg)
              call memory_alloca(mem_modul(1:2,modul),'COCKF_SLD','sld_iniunk',cockf_sld,ndime,nfacg)
              call memory_alloca(mem_modul(1:2,modul),'CRTIP_SLD','sld_iniunk',crtip_sld,ndime,nfacg)
              !
              ! XFEM: initialize crapx: element center and cranx: preexisting crack
              !
              call sld_craini(1_ip)
              call sld_enrich(1_ip)
              !
              ! If the element is CUT
              !
              call cutele(1_ip,cranx_sld,crapx_sld,lcrkf_sld)
              !
              ! Initialize internal variables of cohesive laws and contact/friction
              !
              call sld_updcoh(0_ip)
           end if

           if ( kfl_elcoh > 0 ) then
              !
              ! Initialize vector of nodes where the cohesive law is activated
              !
              call sld_updcoh(0_ip)
           end if

        end if

     else if ( kfl_rigid_sld == 1_ip ) then

        !----------------------------------------------------------------
        !
        ! Rigid body
        !
        !----------------------------------------------------------------

        call sld_rbo_iniunk()

     end if

  else

     !-------------------------------------------------------------------
     !
     ! Read restart file and update unknowns
     !
     !-------------------------------------------------------------------

     call sld_restar(1_ip)
     call sld_updunk(6_ip)

  end if
  !
  ! Coupling initializations
  !
  call sld_coupli(ITASK_INIUNK)

end subroutine sld_iniunk
