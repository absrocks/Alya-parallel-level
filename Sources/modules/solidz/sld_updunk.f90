!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_updunk.f90
!> @author  Guillaume Houzeaux
!> @date    August, 2006
!>          - Subroutine creation
!> @author  Gerard Guillamet and Mariano Vazquez
!> @date    July, 2018
!>          - Subroutine refactoring
!>          - OpenMP parallelization
!> @brief   This routine performs several types of updates for solid's
!>          unknowns.
!> @details
!>          \verbatim
!>          TIME STEP AND ITERATION LABELS (Defined in def_master)
!>            ITER_K       = 1 ..... Current iteration at time step N
!>            ITER_K_STATE = 1 ..... Current iteration at time step N for state variables
!>            ITER_AUX     = 2 ..... Used for coupling iterations
!>            TIME_N       = 3 ..... Time step (Converged)
!>            TIME_N_STATE = 2 ..... Time step (Converged) for state variables
!>
!>          PRE-PROCESS
!>            sld_iniunk (itask=6) .......... if Restart
!>
!>          TIME LOOP
!>            do time
!>              sld_begste (itask=1) ........ (:,ITER_AUX)     <= (:,TIME_N)
!>              do outer
!>                sld_begite (itask=2) ...... (:,ITER_K)       <= (:,ITER_AUX)
!>                do inner
!>                  sld_endite (itask=3) .... (:,ITER_K)       <= UNKNO
!>                end do
!>                sld_endite (itask=4p) ..... (:,ITER_AUX)     <= (:,ITER_K)
!>              end do
!>              sld_endste (itask=5) ........ (:,TIME_N)       <= (:,ITER_K)
!>              sld_endste (itask=7) ........ (:,TIME_N_STATE) <= (:,ITER_K_STATE)
!>            end do
!>          \endverbatim
!> @}
!-----------------------------------------------------------------------

subroutine sld_updunk(itask)

  use def_kintyp,             only : ip,rp
  use def_master,             only : INOTMASTER, displ, unkno, coupling
  use def_master,             only : ITER_AUX, ITER_K, TIME_N
  use def_master,             only : ITER_K_STATE, TIME_N_STATE
  use def_master,             only : ITASK_BEGSTE, ITASK_ENDSTE
  use def_master,             only : ITASK_BEGITE, ITASK_ENDINN, ITASK_ENDITE
  use def_master,             only : ITASK_INIUNK
  use def_domain,             only : npoin, nelem, ndime
  use def_domain,             only : ltype, ngaus
  use mod_memory,             only : memory_size
  use def_solidz,             only : ndofn_sld, nsvar_sld
  use def_solidz,             only : veloc_sld, accel_sld
  use def_solidz,             only : dxfem_sld, vxfem_sld, axfem_sld
  use def_solidz,             only : kfl_xfeme_sld
  use def_solidz,             only : kfl_sdvar_sld, svegm_sld
  use def_solidz,             only : SLD_DYNAMIC_PROBLEM
  use mod_exm_sld_eccoupling, only : exm_sld_ecc_manage_state_variables
  use mod_exm_sld_eccoupling, only : has_exmsld_coupling

  implicit none

  integer(ip), intent(in)       :: itask  !> Update variables at selected case
  integer(ip)                   :: ipoin,itotv,idime,ielem
  integer(ip)                   :: pelty,pgaus,icomp

  if ( INOTMASTER ) then

     select case (itask)

     case ( ITASK_INIUNK )

        !-------------------------------------------------------------------
        !
        ! (:,*) <= (:,1) Initial solution
        !
        !-------------------------------------------------------------------

        do icomp = 1,memory_size(displ,3_ip)
           do ipoin = 1,npoin
              displ(:,ipoin,icomp)     = displ(:,ipoin,ITER_K)
              veloc_sld(:,ipoin,icomp) = veloc_sld(:,ipoin,ITER_K)
              accel_sld(:,ipoin,icomp) = accel_sld(:,ipoin,ITER_K)
           end do
        end do

     case( ITASK_BEGSTE )

        !------------------------------------------------------------------
        !
        !  Initial guess for outer iterations (begin time step)
        !
        !------------------------------------------------------------------

        do ipoin = 1,npoin
           displ(    1:ndime,ipoin,ITER_AUX) = displ(    1:ndime,ipoin,ITER_K)
           veloc_sld(1:ndime,ipoin,ITER_AUX) = veloc_sld(1:ndime,ipoin,ITER_K)
           accel_sld(1:ndime,ipoin,ITER_AUX) = accel_sld(1:ndime,ipoin,ITER_K)
        end do

        ! X-FEM
        if (kfl_xfeme_sld == 1) then
           do ipoin = 1,npoin
              dxfem_sld(1:ndime,ipoin,ITER_AUX) = dxfem_sld(1:ndime,ipoin,ITER_K)
              vxfem_sld(1:ndime,ipoin,ITER_AUX) = vxfem_sld(1:ndime,ipoin,ITER_K)
              axfem_sld(1:ndime,ipoin,ITER_AUX) = axfem_sld(1:ndime,ipoin,ITER_K)
           end do
        end if

     case( ITASK_BEGITE )

        !------------------------------------------------------------------
        !
        !  Initial guess for inner iterations (begin iterations)
        !
        !------------------------------------------------------------------
        !
        ! Assign to UNKNO
        !
        if (kfl_xfeme_sld == 1) then
           !
           ! X-FEM
           !   UNKNO = [uX1 uY1 uZ1 xX1 xY1 xZ1 uX2 ...]
           do ipoin = 1,npoin
              itotv = (ipoin - 1)*ndofn_sld
              do idime =1,ndime
                 itotv = itotv + 1
                 unkno(itotv) = displ(idime,ipoin,ITER_AUX)
              end do
              do idime = 1,ndime
                 itotv = itotv + 1
                 unkno(itotv) = dxfem_sld(idime,ipoin,ITER_AUX)
              end do
           end do

        else
           !
           ! Displacement field
           !   UNKNO = [uX1 uY1 uZ1 uX2 uY2 uZ2 ...
           do ipoin = 1,npoin
              itotv = (ipoin - 1)*ndofn_sld
              do idime=1,ndime
                 itotv = itotv + 1
                 unkno(itotv) = displ(idime,ipoin,ITER_AUX)
              end do
           end do
        end if

     case( ITASK_ENDINN )

        !------------------------------------------------------------------
        !
        !  Update of the displacement
        !
        !------------------------------------------------------------------

        ! X-FEM
        if (kfl_xfeme_sld == 1) then
           !
           ! X-FEM
           !   UNKNO = [uX1 uY1 uZ1 xX1 xY1 xZ1 uX2 ...]
           do ipoin = 1,npoin
              itotv = (ipoin - 1)*ndofn_sld
              do idime =1,ndime
                 itotv = itotv + 1
                 displ(idime,ipoin,ITER_K) = unkno(itotv)
              end do
              do idime = 1,ndime
                 itotv = itotv + 1
                 dxfem_sld(idime,ipoin,ITER_K) = unkno(itotv)
              end do
           end do

        else
           !
           ! Displacement field
           !   UNKNO = [uX1 uY1 uZ1 uX2 uY2 uZ2 ...]
           do ipoin = 1,npoin
              itotv = (ipoin - 1)*ndofn_sld
              do idime=1,ndime
                 itotv = itotv + 1
                 displ(idime,ipoin,1) = unkno(itotv)
              end do
           end do
        end if

     case( ITASK_ENDITE )

        !------------------------------------------------------------------
        !
        !  Updates (end iterations converged)
        !
        !------------------------------------------------------------------

        do ipoin = 1,npoin
           displ(    1:ndime,ipoin,ITER_AUX) = displ(    1:ndime,ipoin,ITER_K)
           veloc_sld(1:ndime,ipoin,ITER_AUX) = veloc_sld(1:ndime,ipoin,ITER_K)
           accel_sld(1:ndime,ipoin,ITER_AUX) = accel_sld(1:ndime,ipoin,ITER_K)
        end do

        ! X-FEM
        if (kfl_xfeme_sld == 1) then
           do ipoin = 1,npoin
              dxfem_sld(:,ipoin,ITER_AUX) = dxfem_sld(:,ipoin,ITER_K)
              vxfem_sld(:,ipoin,ITER_AUX) = vxfem_sld(:,ipoin,ITER_K)
              axfem_sld(:,ipoin,ITER_AUX) = axfem_sld(:,ipoin,ITER_K)
           end do
        end if

     case( ITASK_ENDSTE )

        !------------------------------------------------------------------
        !
        !  End Time Step
        !
        !------------------------------------------------------------------

        do ipoin = 1,npoin
           displ(    1:ndime,ipoin,TIME_N) = displ(    1:ndime,ipoin,ITER_K)
           veloc_sld(1:ndime,ipoin,TIME_N) = veloc_sld(1:ndime,ipoin,ITER_K)
           accel_sld(1:ndime,ipoin,TIME_N) = accel_sld(1:ndime,ipoin,ITER_K)
        end do

        ! X-FEM
        if (kfl_xfeme_sld == 1) then
           do ipoin = 1,npoin
              dxfem_sld(:,ipoin,TIME_N) = dxfem_sld(:,ipoin,ITER_K)
              vxfem_sld(:,ipoin,TIME_N) = vxfem_sld(:,ipoin,ITER_K)
              axfem_sld(:,ipoin,TIME_N) = axfem_sld(:,ipoin,ITER_K)
           end do
        end if
        !
        !  Update State Dependent Variables (SDVs)
        !
        if ( kfl_sdvar_sld == 1_ip ) then

           !--------------------------------------------------
           !$OMP PARALLEL DO SCHEDULE (STATIC)               &
           !$OMP DEFAULT  ( NONE )                           &
           !$OMP PRIVATE  ( ielem, pelty, pgaus )            &
           !$OMP SHARED   ( ltype, ngaus, nelem,             &
           !$OMP            nsvar_sld, svegm_sld )
           !--------------------------------------------------
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              svegm_sld(ielem)%a(1:nsvar_sld,1:pgaus,TIME_N_STATE) = svegm_sld(ielem)%a(1:nsvar_sld,1:pgaus,ITER_K_STATE)
           end do
           !$OMP END PARALLEL DO
           !--------------------------------------------------

        end if
        !
        ! Update coupling state variables
        !
        if( has_exmsld_coupling() )then
            call exm_sld_ecc_manage_state_variables(ITASK_ENDSTE)
        endif

     end select

  end if

end subroutine sld_updunk

