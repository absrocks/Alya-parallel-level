!-----------------------------------------------------------------------
!> @addtogroup Begste
!> @{
!> @file    Begste.f90
!> @author  Guillaume Houzeaux
!> @brief   Begin a time step
!> @details Begin a time step:
!>          - Modules update boundary conditions, unknowns, etc.
!> @} 
!-----------------------------------------------------------------------
subroutine Begste()
  use def_kintyp,    only : ip
  use def_master
  use def_parame
  use mod_ker_proper
  use def_coupli,    only : mcoup
  use def_coupli,    only : coupling_driver_iteration
  use mod_couplings, only : COU_TEMPORAL_PREDICTOR
  use mod_messages, only : livinf
  implicit none 

  integer(ip)            :: icoup
  !
  ! Turn back reset flag to previous step
  !
  if (kfl_reset == 1) then
     !
     ! Initializations
     !
     call iniste(2_ip)
     cutim  = cutim - dtime
     call setgts(2_ip)
     call livinf(201_ip, ' ',1_ip)
  endif
  !
  ! Temporal predction for zonal coupling (only call it after two time steps and before the modules)
  !
  if( itti2 > 2 )then 
     do icoup = 1_ip, mcoup
        call COU_TEMPORAL_PREDICTOR(icoup)
     end do
  end if
  !
  ! Begin a time step for each module
  !
  call moduls(ITASK_BEGSTE)
  !
  ! Turn back properties if reset
  !
  if (kfl_reset == 1) then
     call ker_updpro()
     ! Deactivate reset request
     kfl_reset = 0
  end if
  !
  ! Coupling: Put counters to zero
  ! 
  if( mcoup > 0 ) then
     coupling_driver_iteration(1:nblok) = 0
  end if

!!call coupling_set_dt()
  !
  ! Transient fields
  !
  call calc_kx_tran_fiel()

      !  if( trim(title)=='falling') call chimihole()

  contains
!-------------------------------------------------------------------------||---!
!                                                                              !
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
    subroutine coupling_set_dt()
      use mod_communications, only : par_min, par_max  
      use mod_couplings,      only : THERE_EXISTS_A_ZONE_COUPLING
      use def_coupli,         only : mcoup, coudt
      implicit none 
      real(rp) :: dt_inv(1_ip)
      !
      !
      !print *, coudt, "<-----"
      !
      if( (mcoup > 0).and.(coudt==1_ip).and.THERE_EXISTS_A_ZONE_COUPLING() ) then 
         dt_inv(1_ip) = dtinv  
         call PAR_MAX(1_ip, dt_inv, 'IN CURRENT COUPLING') !< 2016Mar23 
         !call PAR_MIN(1_ip, dt_inv, 'IN CURRENT COUPLING') 
         !
         !    if(imaster) print *, 1.0/dtinv, 1.0/dt_inv(1_ip), "<----"
         !
         dtinv = dt_inv(1_ip)      
         !
         cutim  = cutim - dtime          ! \
         !call iniste(2_ip)               ! |__ 2014Dic10, <-- Begste, if(kfl_reset == 1) 
         call setgts(2_ip)               ! |
         call livinf(201_ip, ' ',1_ip)   ! /
         !
      endif
      ! 

    end subroutine coupling_set_dt
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!

  end subroutine Begste


subroutine chimihole()
  use def_master
  use def_domain
  use def_elmtyp
  use mod_memory
  use mod_communications
  implicit none
  integer(ip)          :: ipoin,inode,knode,ielem,pnode,kpoin
  real(rp)             :: xc(3),dista
  logical(lg), pointer :: lmark(:)
  integer(ip), pointer :: lnsub(:)

  nullify(lmark)
  nullify(lnsub)

  !xc(1) = 0.5_rp
  !xc(2) = 0.5_rp - 0.5_rp * cutim
  !xc(3) = 0.0_rp

  xc(1) = 0.0_rp
  xc(2) = 0.0_rp
  xc(3) = 8.0_rp - 1.0_rp * cutim

  call memory_alloca(memor_dom,'LMARK','memgeo',lmark,npoin)
  call memory_alloca(memor_dom,'LNSUB','memgeo',lnsub,npoin)

  do ielem = 1,nelem
     do inode = 1,lnnod(ielem)
        ipoin = lnods(inode,ielem)
        lnsub(ipoin) = lesub(ielem)
        lnoch(ipoin) = NOFEM
     end do
  end do
 
  kpoin = 0
  do ipoin = 1,npoin
     if( lnsub(ipoin) == 1 ) then
        Dista = sqrt( dot_product(coord(1:ndime,ipoin)-xc(1:ndime),coord(1:ndime,ipoin)-xc(1:ndime)) ) 
        if( dista <= 0.4_rp ) then
           lmark(ipoin) = .true.
           kpoin = kpoin + 1
        end if
     end if
  end do

  call PAR_SUM(kpoin)
  call PAR_INTERFACE_NODE_EXCHANGE(lmark,'OR','IN MY CODE')

  do ielem = 1,nelem
     ltype(ielem) = abs(ltype(ielem))
     lelch(ielem) = ELFEM
     knode = 0
     pnode = lnnod(ielem)
     do inode = 1,pnode
        ipoin = lnods(inode,ielem)
        if( lmark(ipoin) ) knode = knode + 1
     end do
     if( knode == pnode ) then
        lelch(ielem) =  ELHOL
        ltype(ielem) = -ltype(ielem)
        do inode = 1,pnode
           ipoin = lnods(inode,ielem)
           lnoch(ipoin) = NOHOL
        end do
     end if     
  end do

  call memory_deallo(memor_dom,'LNSUB','memgeo',lnsub)
  call memory_deallo(memor_dom,'LMARK','memgeo',lmark)
  
  call par_element_loop()

end subroutine chimihole
