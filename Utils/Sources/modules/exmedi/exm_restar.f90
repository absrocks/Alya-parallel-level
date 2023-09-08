!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_restar.f90
!> @date    16/11/1966
!> @author  Mariano Vazquez
!> @brief   Restart subroutine
!> @details Restart subroutine
!> @}
!------------------------------------------------------------------------
subroutine exm_restar(itask)
  use def_parame
  use def_master
  use def_domain
  use def_exmedi
  use mod_postpr
  use mod_memchk
  use mod_messages, only : livinf
  implicit none
  integer(ip), intent(in) :: itask !> read or write
  integer(ip)             :: icomp,icomp_cell,jcomp,kfl_gores,iconc,iauxi,iicel
  character(5)            :: wopoe(3)
  character(300)          :: messa
  !
  ! Check if restart file should be read or written
  !
  call respre(itask,kfl_gores)
  if( kfl_gores == 0 ) return

  if( itask == READ_RESTART_FILE ) then
     icomp = ncomp_exm
     icomp_cell = 2
  else
     icomp = 1
     icomp_cell = 1
  end if

  !----------------------------------------------------------------------
  !
  ! Variables
  !
  !----------------------------------------------------------------------

!  vauxi_exm(:,:2)
!  vconc(:,:,2)
!  elmag(:,icomp)
  
  ! all exmedi restart variables are scalars and defined as node-wise 
  wopoe(2) = 'SCALA'
  wopoe(3) = 'NPOIN'

  wopoe(1) = 'INTRA'
  messa = '        VARIABLE: ELMAG -> INTRA - TYPE: '//wopoe(2)//' '//wopoe(3)
  if( INOTSLAVE  ) call livinf(0_ip,messa,one)
  if( INOTMASTER ) gesca => elmag(:,icomp)
  call postpr(gesca,wopoe,ittim,cutim)

  messa = '        VARIABLE: VCONC -> VCO01 ... VCO'//trim(intost(nconc_exm))//' - TYPE: '//wopoe(2)//' '//wopoe(3)
  if( INOTSLAVE  ) call livinf(0_ip,messa,one)

  do iconc= 1,nconc_exm
     if (iconc.le.9) then
        wopoe(1) = 'VCO0'//trim(intost(iconc))
     else
        wopoe(1) = 'VCO'//trim(intost(iconc))
     end if
     messa = '        VARIABLE: '//wopoe(1)//' - TYPE: '//wopoe(2)//' '//wopoe(3)
     if( INOTMASTER ) gesca => vconc(iconc,1:npoin,icomp_cell)
     call postpr(gesca,wopoe,ittim,cutim)
  end do

  messa = '        VARIABLE: VAUXI_EXM -> VAU01 ... VAU'//trim(intost(nauxi_exm))//' - TYPE: '//wopoe(2)//' '//wopoe(3)
  if( INOTSLAVE  ) call livinf(0_ip,messa,one)

  do iauxi= 1,nauxi_exm
     if (iauxi.le.9) then
        wopoe(1) = 'VAU0'//trim(intost(iauxi))
     else
        wopoe(1) = 'VAU'//trim(intost(iauxi))
     end if
     messa = '        VARIABLE: '//wopoe(1)//' - TYPE: '//wopoe(2)//' '//wopoe(3)
     if( INOTMASTER ) gesca => vauxi_exm(iauxi,1:npoin,icomp_cell)
     call postpr(gesca,wopoe,ittim,cutim)
  end do

  messa = '        VARIABLE: VICEL_EXM -> VIC01 ... VIC'//trim(intost(nicel_exm))//' - TYPE: '//wopoe(2)//' '//wopoe(3)
  if( INOTSLAVE  ) call livinf(0_ip,messa,one)

  do iicel= 1,nicel_exm
     if (iicel.le.9) then
        wopoe(1) = 'VIC0'//trim(intost(iicel))
     else
        wopoe(1) = 'VIC'//trim(intost(iicel))
     end if
     messa = '        VARIABLE: '//wopoe(1)//' - TYPE: '//wopoe(2)//' '//wopoe(3)
     if( INOTMASTER ) gesca => vicel_exm(iicel,1:npoin,icomp_cell)
     call postpr(gesca,wopoe,ittim,cutim)
  end do

  if( itask == READ_RESTART_FILE ) then
     if (INOTMASTER) then
        do jcomp=1,ncomp_exm-1
           elmag(1:npoin,ncomp_exm - jcomp) = elmag(1:npoin,icomp)
        end do
        do iconc= 1,nconc_exm
           vconc(iconc,1:npoin,1)= vconc(iconc,1:npoin,icomp_cell)
        end do
        do iauxi= 1,nauxi_exm
           vauxi_exm(iauxi,1:npoin,1)= vauxi_exm(iauxi,1:npoin,icomp_cell)
        end do
        do iicel= 1,nicel_exm
           vicel_exm(iicel,1:npoin,1)= vicel_exm(iicel,1:npoin,icomp_cell)
        end do
     end if
  end if
  
  
!  wopoe(1) = 'REFHN'
!  wopoe(2) = 'SCALA'
!  if( INOTMASTER ) gesca => refhn_exm(:,icomp)
!  call postpr(gesca,wopoe,ittim,cutim)     

  !----------------------------------------------------------------------
  !
  ! Finish
  !
  !----------------------------------------------------------------------

  call respre(3_ip,kfl_gores)

end subroutine exm_restar
 
