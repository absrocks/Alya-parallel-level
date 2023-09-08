!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_outvar.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Bridge to the kernel's posptrocess
!> @details Bridge to the kernel's posptrocess
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_outvar(ivari)
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal
  use      def_postpr
  use      mod_memchk
  use      mod_postpr
  use      mod_iofile

use mod_lodi_nsa

  implicit none
  integer(ip), intent(in)  :: ivari
  integer(ip)              :: ipoin,kpoin,idime
  real(rp)                 :: rauxi,rmodu,xadim,xadi2,xdelo,xprlo,xtelo,auxi

  select case ( ivari )  

  case ( 0_ip )
     !
     ! Do not postprocess anything
     !
     return

  case ( 1_ip)
     !
     ! VELOC
     !    
     gevec => veloc(:,:,1)

  case ( 2_ip )
     !
     ! UMOME
     !
     gevec => umome(:,:,1)

  case ( 4_ip )
     !
     ! VORTICITY
     !
     if(ndime==2) then
        if( INOTMASTER ) then
           call memgen(zero,npoin,0_ip)
           do ipoin = 1,npoin
              gesca(ipoin)=sqrt(&
                   &  vorti(1,ipoin)*vorti(1,ipoin)&
                   & +vorti(2,ipoin)*vorti(2,ipoin))
           end do
        end if
     else
        if( INOTMASTER ) then
           gevec => vorti
        end if
     end if

  case ( 20_ip )
     !
     ! PRESS
     !
     gesca =>  press(:,1)

  case ( 21_ip )
     !
     ! CDENS or DENSI
     ! 
     if ( INOTMASTER ) then
        if(kfl_pertu_nsa == 1) then
           do ipoin = 1,npoin
              xdelo= rekee_nsa(ndime+1    ,ipoin)
              unkno(ipoin)= densi(ipoin,1) + xdelo
           end do
        else
           do ipoin = 1,npoin
              unkno(ipoin)= densi(ipoin,1)
           end do
        end if
     end if
     gesca => unkno

  case ( 22_ip )
     !
     ! TEMPE
     !
     if ( INOTMASTER ) then
        if(kfl_pertu_nsa == 1) then
           do ipoin = 1,npoin
              xtelo= rekee_nsa(ndime+2    ,ipoin)
              unkno(ipoin)= tempe(ipoin,1) + xtelo
           end do
        else
           do ipoin = 1,npoin
              unkno(ipoin)= tempe(ipoin,1)
           end do
        end if
     end if
     gesca => unkno
     !gesca => tempe(:,1)

  case ( 23_ip ) 
     !
     ! ENERG
     !
     gesca => energ(:,1)

  case ( 24_ip ) 
     !
     ! VISCO
     !
     gesca => visco(:,1)

  case ( 25_ip ) 
     !
     ! MACH
     !
     if( INOTMASTER ) call nsa_outmac(0,0,1,rauxi)
     gesca => unkno

  case ( 26_ip )
     !
     ! ENTHA
     !
     if (kfl_relat_nsa == 1) then   ! enthalpy only for relativistic flows
        gesca => entha(:,1)
     else
        return
     end if

  case ( 27_ip ) 
     !
     ! T in Celsius
     !
     if( INOTMASTER ) then        
        do ipoin = 1,npoin
           unkno(ipoin) = tempe(ipoin,1) - 273.0_rp         ! to output T in Celsius
        end do
     end if
     gesca => unkno

  case ( 28_ip ) 
     !
     !
     !
     if( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin) = press(ipoin,1) / densi(ipoin,1) ** adgam_nsa
        end do
     end if
     gesca => unkno

  case ( 29_ip ) 
     !
     ! VORTI: CHECK!
     !
     if( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin) = vorti(1,ipoin)*xadim
        end do
     end if
     gesca => unkno

  case ( 30_ip ) 
     !
     ! RMODU
     !
     if( INOTMASTER ) then
        if ( ndime == 1 ) then
           do ipoin = 1,npoin
              rmodu        = veloc(1,ipoin,1)*veloc(1,ipoin,1)
              unkno(ipoin) = sqrt(rmodu)
           end do
        else if ( ndime == 2 ) then
           do ipoin = 1,npoin
              rmodu        = veloc(1,ipoin,1)*veloc(1,ipoin,1) + &
                   &         veloc(2,ipoin,1)*veloc(2,ipoin,1)
              unkno(ipoin) = sqrt(rmodu)
           end do
        else if ( ndime == 3 ) then
           do ipoin = 1,npoin
              rmodu        = veloc(1,ipoin,1)*veloc(1,ipoin,1) + &
                   &         veloc(2,ipoin,1)*veloc(2,ipoin,1) + &
                   &         veloc(3,ipoin,1)*veloc(3,ipoin,1)
              unkno(ipoin) = sqrt(rmodu)
           end do
        end if
     end if
     gesca => unkno

  case ( 31_ip ) 
     !
     ! RMODU
     !
     if( INOTMASTER ) then
        do ipoin = 1,npoin
           if (ndime==1 .or. ndime==2) then
              rmodu= vorti(1,ipoin)*vorti(1,ipoin)
           else if (ndime==3) then
              rmodu= vorti(1,ipoin)*vorti(1,ipoin) + &
                   vorti(2,ipoin)*vorti(2,ipoin) + vorti(3,ipoin)*vorti(3,ipoin)
           end if
           unkno(ipoin)= sqrt(rmodu)*xadim
        end do
     end if
     gesca => unkno

  case ( 32_ip ) 
     !
     !
     !
     if( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin)= vorti(ndime+1,ipoin)*xadi2
        end do
     end if
     gesca => unkno

  case ( 33_ip )
     !
     ! VELOX
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin)= veloc(1,ipoin,1)
        end do
     end if
     gesca => unkno

  case ( 34_ip )              
     !
     ! VELOY
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin)= veloc(2,ipoin,1)
        end do
     end if
     gesca => unkno

  case ( 35_ip )              
     !
     ! VELOZ
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin)= veloc(3,ipoin,1)
        end do
     end if
     gesca => unkno

  case ( 40_ip ) 
     !
     ! PRESS - XPRLO
     !    
     if ( INOTMASTER ) then
        do ipoin = 1,npoin

           xprlo= rekee_nsa(ndime+3,ipoin)
           unkno(ipoin)= press(ipoin,1) - xprlo
        end do
     end if
     gesca => unkno

  case ( 41_ip ) 
     !
     ! DENSI - XDELO
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           xdelo= rekee_nsa(ndime+1,ipoin)
           unkno(ipoin)= densi(ipoin,1) - (1.0_rp - kfl_pertu_nsa)*xdelo      
        end do
     end if
     gesca => unkno

  case ( 42_ip ) 
     !
     ! THBAS
     !
     if ( INOTMASTER ) then

        do ipoin = 1,npoin
           xtelo= rekee_nsa(ndime+2    ,ipoin)
           unkno(ipoin)= tempe(ipoin,1) - (1.0_rp - kfl_pertu_nsa)*xtelo
        end do
     end if
     gesca => unkno

  case ( 43_ip ) 
     !
     ! SGSMX
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin) = umoss_nsa(1,ipoin,2)
        end do
     end if
     gesca => unkno

  case ( 44_ip ) 
     !
     ! SGSMY
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin) = umoss_nsa(ndime,ipoin,2)
        end do
     end if
     gesca => unkno

  case ( 45_ip ) 
     !
     ! SGSDE
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin) = denss_nsa(ipoin,2)
        end do
     end if
     gesca => unkno

  case ( 46_ip ) 
     !
     ! SGSEN (energy or theta sgs)
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin) = eness_nsa(ipoin,2)
        end do
     end if
     gesca => unkno

  case ( 47_ip ) 
     !
     ! SGS MOMENTUM (as vector)
     !
     gevec => umoss_nsa(:,:,2)

  case ( 48_ip )
     !
     ! XMOME
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin)= umome(1,ipoin,1)
        end do
     end if
     gesca => unkno

  case ( 49_ip )
     !
     ! YMOME
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin)= umome(2,ipoin,1)
        end do
     end if
     gesca => unkno

  case ( 50_ip )
     !
     ! ZMOME
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin)= umome(3,ipoin,1)
        end do
     end if
     gesca => unkno

  case ( 51_ip ) 
     !
     ! FREQU
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin) = frequ_nsa(ipoin)
        end do
     end if
     gesca => unkno

  case ( 54_ip ) 
     !
     ! NUTUR
     !

     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin) = turmu(ipoin)
        end do
     end if
     call rhsmod(1_ip,unkno)
     gesca => unkno
  case ( 55_ip ) 
     !
     !  AVVEL
     !
     if (cutim > avtim_nsa) then
        auxi = cutim - avtim_nsa
        if ( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin = 1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin)     = avvel_nsa(idime,ipoin) / auxi
                 avvel_nsa(idime,ipoin) = 0.0_rp
              enddo
           end do
        endif
     else
        return
     endif

  case ( 56_ip ) 
     !
     !  AVVE2
     !
     if (cutim > avtim_nsa) then
        auxi = cutim - avtim_nsa
        if ( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin = 1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin)     = avve2_nsa(idime,ipoin) / auxi
                 avve2_nsa(idime,ipoin) = 0.0_rp
              enddo
           end do
        endif
     else
        return
     endif

  case ( 57_ip ) 
     !
     !  AVVXY
     !
     if (cutim > avtim_nsa) then
        auxi = cutim - avtim_nsa
        if ( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin = 1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin)     = avvxy_nsa(idime,ipoin) / auxi
                 avvxy_nsa(idime,ipoin) = 0.0_rp
              enddo
           end do
        endif
     else
        return
     endif

  case ( 58_ip ) 
     !
     !  AVMOM
     !
     if (cutim > avtim_nsa) then
        auxi = cutim - avtim_nsa
        if ( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin = 1,npoin
              do idime=1,ndime
                 gevec(idime,ipoin)     = avmom_nsa(idime,ipoin) / auxi
                 avmom_nsa(idime,ipoin) = 0.0_rp
              enddo
           end do
        endif
     else
        return
     endif

  case ( 59_ip ) 
     !
     !  AVPRE
     !
     if (cutim > avtim_nsa) then
        auxi = cutim - avtim_nsa
        if ( INOTMASTER ) then
           do ipoin = 1,npoin
              unkno(ipoin)     = avpre_nsa(ipoin) / auxi
              avpre_nsa(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

  case ( 60_ip ) 
     !
     !  AVTEM
     !
     if (cutim > avtim_nsa) then
        auxi = cutim - avtim_nsa
        if ( INOTMASTER ) then
           do ipoin = 1,npoin
              unkno(ipoin)     = avtem_nsa(ipoin) / auxi
              avtem_nsa(ipoin) = 0.0_rp
           end do
           gesca => unkno
        endif
     else
        return
     endif

     !-----------------------------------------------------------------------!
     ! LODI
     !-----------------------------------------------------------------------! 
  case( 61_ip )
!     if(INOTMASTER) then
!        !
!        CHARAs%gamme = adgam_nsa
!        call lodi_nsa_allocate( CHARAs )
!        call lodi_nsa_get_characteristics( CHARAs )
!        !
!        call memgen(zero,ndime,npoin)
!        do idime = 1,ndime
!          do kpoin = 1,npoiz(izone_nsa)
!            ipoin = lpoiz(izone_nsa)%l(kpoin) 
!            gevec(idime,ipoin) = CHARAs%chrc(CHARAs%idofn,ipoin,idime) ! chrc_nsa(ndofn_nsa,npoin,ndime) 
!          enddo 
!        enddo
!        !
!        call lodi_nsa_deallocate( CHARAs )
!        !
!     endif
     call get_characteristic( 1_ip )


  case( 62_ip )
     call get_characteristic( 2_ip )
     !if(INOTMASTER) then 
     !   !do kpoin = 1,npoiz(izone_nsa)
     !   !  ipoin = lpoiz(izone_nsa) % l(kpoin)
     !   !  unkno(ipoin) = schlrn_nsa(ipoin)
     !   !enddo 
     !   gesca => unkno
     !endif
     !-----------------------------------------------------------------------!

  case( 63_ip )
     if(INOTMASTER) then 
        call memgen(zero,npoin,0_ip)
!        do kpoin = 1,npoiz(izone_nsa)
!           ipoin = lpoiz(izone_nsa) % l(kpoin)
!           gesca(ipoin) = crens_nsa(1,ipoin)
!        end do

!!!        call nsa_outvar_residual(1_ip,npoin,gesca)

     endif
     call get_characteristic( ndime )

  case( 64_ip )
     if(INOTMASTER) then 
        call memgen(zero,npoin,0_ip)
!        do kpoin = 1,npoiz(izone_nsa)
!           ipoin = lpoiz(izone_nsa) % l(kpoin)
!           gesca(ipoin) = crens_nsa(2,ipoin)
!        end do

!!!        call nsa_outvar_residual(2_ip,npoin,gesca)

     endif
     call get_characteristic( ndime+1_ip )

  case( 65_ip )
     if(INOTMASTER) then 
        call memgen(zero,npoin,0_ip)
!        do kpoin = 1,npoiz(izone_nsa)
!           ipoin = lpoiz(izone_nsa) % l(kpoin)
!           gesca(ipoin) = crens_nsa(3,ipoin)
!        end do

!!!!        call nsa_outvar_residual(3_ip,npoin,gesca)

     endif
     call get_characteristic( ndime+2_ip )

  case( 66_ip )
     if(INOTMASTER) then 
        call memgen(zero,npoin,0_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = dtieq_nsa(1,ipoin,1)
        end do
     endif

  case( 67_ip )
     if(INOTMASTER) then 
        call memgen(zero,npoin,0_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = dtieq_nsa(2,ipoin,1)
        end do
     endif

  case( 68_ip )
     if(INOTMASTER) then 
        call memgen(zero,npoin,0_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = dtieq_nsa(3,ipoin,1)
        end do
     endif

  case ( 69_ip )
     !
     ! TEMPE
     !
     if ( INOTMASTER ) then
        do ipoin = 1,npoin
           unkno(ipoin)= (press(ipoin,1) - press_nsa) / press_dynamic_nsa 
        end do
     end if
     gesca => unkno
     !gesca => tempe(:,1)

  case(70_ip)
     !
     ! GROUP: GROUPS
     !
     if( INOTMASTER ) then
        do ipoin=1,npoin
           rhsid(ipoin)=real(solve(1)%lgrou(ipoin))
        end do
        gesca => rhsid
     end if

  end select
  !
  ! Postprocess
  !
  call outvar(&
       ivari,&
       ittim,cutim,postp(1) % wopos(1,ivari))


contains

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine get_characteristic( idofn )
  implicit none
  integer(ip), intent(in) :: idofn
  !
  if(INOTMASTER) then
        !
        CHARAs%gamme = adgam_nsa
        call lodi_nsa_allocate( CHARAs )
        call lodi_nsa_get_characteristics( CHARAs )
        !
        call memgen(zero,ndime,npoin)
        do idime = 1,ndime
          do ipoin = 1,npoin
            gevec(idime,ipoin) = CHARAs%chrc(idofn,ipoin,idime) ! chrc_nsa(ndofn_nsa,npoin,ndime) 
          enddo
        enddo
        !
        call lodi_nsa_deallocate( CHARAs )
        !
  endif
  !
  end subroutine

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

!!$  subroutine nsa_outvar_residual(itask,ngisc,gisca)
!!$    use      def_parame
!!$    use      def_master
!!$    use      def_domain
!!$    use      def_nastal
!!$
!!$    implicit none
!!$    integer(ip)    :: itask !> 
!!$    integer(ip)    :: ngisc !> 
!!$    integer(ip)    :: kpoin,ipoin,idime
!!$    real(rp)       :: gisca(ngisc)  !> ngisc is npoin (so far)
!!$    real(rp)       :: va,vd,va_plus,vd_plus,resid(3,2),debug_resi(3)
!!$    
!!$    if (itask == 1) then
!!$
!!$       do kpoin = 1,npoiz(izone_nsa)
!!$          ipoin = lpoiz(izone_nsa) % l(kpoin)
!!$          va_plus = 0.0_rp
!!$          vd_plus = 0.0_rp
!!$          do idime=1,ndime
!!$             vd = umome(idime,ipoin,1) - umome(idime,ipoin,3)
!!$             va = umome(idime,ipoin,1)
!!$             va_plus = va_plus + va * va
!!$             vd_plus = vd_plus + vd * vd
!!$          end do
!!$          gisca(ipoin) = sqrt(vd_plus)
!!$       end do
!!$       
!!$    else if (itask == 2) then
!!$       
!!$       do kpoin = 1,npoiz(izone_nsa)
!!$          vd = densi(ipoin,1) - densi(ipoin,3)
!!$          va = densi(ipoin,1)
!!$          va_plus = va * va
!!$          vd_plus = vd * vd
!!$          gisca(ipoin) = sqrt(vd_plus)
!!$       end do
!!$       
!!$    else if (itask == 3) then
!!$       
!!$       do kpoin = 1,npoiz(izone_nsa)
!!$          vd = energ(ipoin,1) - energ(ipoin,3)
!!$          va = energ(ipoin,1)
!!$          va_plus = va * va
!!$          vd_plus = vd * vd
!!$          gisca(ipoin) = sqrt(vd_plus)
!!$       end do
!!$       
!!$    end if
!!$
!!$
!!$
!!$  end subroutine nsa_outvar_residual

end subroutine nsa_outvar
    
