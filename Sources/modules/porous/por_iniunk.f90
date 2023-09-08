!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_iniunk.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Sets up the initial condition
!> @details Sets up the initial condition
!> @} 
!------------------------------------------------------------------------
subroutine por_iniunk()
  use def_parame
  use def_master
  use def_domain
  use def_porous
  use def_kermod
  use mod_ker_proper 
  implicit none
  integer(ip) :: ipoin,icomp,ielem,iwell,iheiw

  kprsa_por = 1_ip
  if( kfl_rstar == 0 ) then  

     !-------------------------------------------------------------------
     !
     ! Normal run
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then
        !
        ! Initial conditions for pressure and saturation - also obtain nodal porosity
        !
        icomp = min(3_ip,ncomp_por)
        !
        ! For the moment I use xfiel(1) for Satin,  xfiel(2) for poro0 ...
        do ielem = 1,nelem
           satin_por(ielem)         = xfiel(1) % a(1,ielem,1)
           poro0_por(ielem)         = xfiel(2) % a(1,ielem,1)
           perme_por(1:ndime,ielem) = xfiel(3) % a(1:ndime,ielem,1)
        end do
        do ipoin = 1,npoin
           iwell_por(ipoin) = max(nint(xfiel(4)%a(1,ipoin,1)),0_ip)
           !iheip_por(ipoin) = nint(xfiel(4)%a(2,ipoin))
        end do

        do ipoin = 1,npoin
           iwell = iwell_por(ipoin)
           if (iwell /= 0_ip) then
              ! ipwel_por(iheip_por(ipoin),iwell) = ipoin
           end if
        end do
        !
        ! Check that ipwel_por has been set for all nheiw_por points
        !
        do iwell=1,nwell_por
           do iheiw=1,nheiw_por(iwell)
              ! if ( ipwel_por(iheiw,iwell) == 0 ) then
              !    print*,'iheiw,iwell,ipwel_por(iheiw,iwell),nheiw_por(iwell)',iheiw,iwell,ipwel_por(iheiw,iwell),nheiw_por(iwell)
              !    call runend('POR_INIUNK:ipwel_por(iheiw,iwell) == 0 ')
              ! end if
           end do
        end do

        call smoot2(satin_por,wasat(1,1))
        call smoot2(poro0_por,nodpo_por)
        call smoot3(perme_por,nodpe_por,ndime)
        do ipoin = 1,npoin
           wasat(ipoin,icomp) = wasat(ipoin,1)
           press(ipoin,icomp) = prini_por
           press(ipoin,1) = prini_por
        end do
     end if

     if ( grnor_por /= 0.0_rp ) then
        call por_inipre()
        if( INOTMASTER ) then
           do ipoin = 1,npoin
              press(ipoin,icomp) = press(ipoin,1)
           end do
        end if
     end if
     if( INOTMASTER ) then
        do ipoin = 1,npoin
           if( kfl_fixno_por(1,ipoin,1) == 1 ) press(ipoin,icomp) = bvess_por(1,ipoin,1)
           if( kfl_fixno_por(1,ipoin,2) == 1 ) wasat(ipoin,icomp) = bvess_por(1,ipoin,2)
        end do
     end if

  else

     !-------------------------------------------------------------------
     !
     ! Read restart file
     !
     !-------------------------------------------------------------------

     !     call por_restar(1_ip)  ! for the moment restart not ready

  end if

  call por_wellin()

end subroutine por_iniunk
