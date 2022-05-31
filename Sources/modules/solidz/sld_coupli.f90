!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_coupli.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   This routine manages solidz coupling with IMMBOU, ALEFOR, NASTIN...
!> @details This routine manages solidz coupling with IMMBOU, ALEFOR, NASTIN...
!> @} 
!-----------------------------------------------------------------------
subroutine sld_coupli(itask)
  use def_master
  use def_elmtyp
  use def_domain
  use def_solidz
  implicit none
  integer(ip), intent(in) :: itask   !< who is calling sld_coupli
  integer(ip)             :: ipoin,kpoin,iimbo,idime,jdime,idofn,incnt
  integer(ip)             :: ipoin_fluid,ipoin_solid
  real(rp)                :: foref, aitkw, prev_coupling_dispm, xfact(2), deltu, rdiff

  select case ( itask )
     
  case (ITASK_INIUNK )

  case (ITASK_BEGITE )
     !
     ! Coupling with IMMBOU: prescribe displacement
     !
     if( INOTMASTER ) then
        if( coupling('SOLIDZ','IMMBOU') >= 1 ) then
           do iimbo = 1,nimbo
              if( imbou(iimbo) % kfl_typeb >= 1 ) then
                 do kpoin = 1,imbou(iimbo) % npoib
                    ipoin = imbou(iimbo) % lninv(kpoin)
                    do idime = 1,ndime    
                       bvess_sld(idime,ipoin,1) = imbou(iimbo) % cooib(idime,kpoin) - coord(idime,ipoin)
                    end do
                 end do
              end if
           end do
        end if
     end if

  case (ITASK_ENDITE )
     !
     ! Coupling with ALEFOR: prescribe displacement
     !


     !
     ! Aitken relaxation factor on the contacts (NOTE: CORRECTION DONE BY WALL, RELCO COMING FROM AITKEN SUBRU IS AITKW)
     !
     !     aitkw= (1.0_rp - relco_sld)
     !     if (micou(iblok) == 1) aitkw= 1.0_rp      ! there are no coupling iterations in this block
     aitkw = relco_sld
     if( micou(iblok) == 1 ) aitkw = 1.0_rp      ! there are no coupling iterations in this block
     !     write(6,*)
     !     write(6,*) 'totoooo', itcou, aitkw
     !     write(6,*)

  case ( ITASK_MATRIX )

     !-------------------------------------------------------------------
     !
     ! MATRIX: After assembly
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then

        if( coupling('SOLIDZ','IMMBOU') >= 1 ) then
           !
           ! Coupling with Immbou: impose force FORCF, interpolate b.c.
           !     
           call sld_immbou(amatr,rhsid,unkno)

        end if
        !
        ! <GGU> Obsolete 
        !
        !if( coupling('SOLIDZ','NASTIN') >= 1 ) then
        !   !
        !   ! Add forces in the reference frame foref
        !   !
        !   call memgen(1_ip,npoin,0_ip)
        !   do ipoin = 1,npoin
        !      if( lnoch(ipoin) == NODE_CONTACT_SOLID ) then
        !         if( ipoin <= npoi1 .or. ( ipoin >= npoi2 .and. ipoin <= npoi3 ) ) then  
        !            gisca(ipoin) = 1
        !         end if
        !      end if
        !   end do
        !   do kpoin = 1,npoi4
        !      ipoin = lpoi4(kpoin)
        !      if( lnoch(ipoin) == NODE_CONTACT_SOLID ) then
        !         gisca(ipoin) = 1
        !      end if
        !   end do
        !   do ipoin = 1,npoin
        !      if( gisca(ipoin) == 1 ) then
        !         idofn = (ipoin-1) * ndime
        !         do idime = 1,ndime
        !            if( kfl_fixno_sld(idime,ipoin) /= 1 ) then
        !               foref = 0.0_rp
        !               do jdime = 1,ndime
        !                  foref = foref + gdepo(idime,jdime,ipoin) * forcf(jdime,ipoin)
        !               end do
        !               idofn        = idofn + 1
        !               rhsid(idofn) = rhsid(idofn) + foref
        !            end if
        !         end do
        !      end if
        !   end do
        !   call memgen(3_ip,npoin,0_ip)

!!$           do incnt = 1,nncnt
!!$              ipoin_solid = lncnt(2,incnt)
!!$              idofn       = (ipoin_solid-1)*ndime
!!$              do idime = 1,ndime
!!$                 pipix(idime) = pipix(idime) + forcf(idime,ipoin_solid)
!!$                 foref = 0.0_rp
!!$                 do jdime = 1,ndime
!!$                    foref = foref + gdepo(idime,jdime,ipoin_solid) * forcf(jdime,ipoin_solid)
!!$                 end do
!!$                 idofn        = idofn + 1
!!$                 rhsid(idofn) = rhsid(idofn) + foref
!!$                 pipox(idime) = pipox(idime) + foref
!!$              end do
!!$           end do
!!$        end if

        !end if

     end if

  end select

!===============================================================| contains |===!
contains
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_sld_get_svalu( svalu, dummy )
  use def_solidz, only: forc_fsi                       !< (ndime,npoin) used FSI
  use def_master, only: displ                          !< (:,:,:)
  implicit none
  real(rp), intent(out) :: svalu(:,:)
  real(rp), intent(out) :: dummy(:,:)
  !
  svalu(1:ndime,1:npoin) = 0.0_rp
  !
  if(INOTMASTER) then
    svalu(1:ndime,1:npoin) = displ(1:ndime,1:npoin,1) - displ(1:ndime,1:npoin,3)
  endif
  !
  dummy(1:ndime,1:npoin) = 0.0_rp
  !
  !
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_sld_set_force( dummy, xvalu )
  use def_solidz, only: forc_fsi                       !< (ndime,npoin) used FSI
  use def_master, only: displ                          !< (:,:,:)
  implicit none
  real(rp), intent(out)   :: dummy(:,:)
  real(rp), intent(inout) :: xvalu(:,:)
  !
  if(INOTMASTER) then
    forc_fsi(1:ndime,1:npoin) = xvalu(1:ndime,1:npoin)
  endif
  !
  dummy(1:ndime,1:npoin) = 0.0_rp
  xvalu(1:ndime,1:npoin) = 0.0_rp
  !
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_sld_matrix( rhs_id )
  use def_solidz, only: forc_fsi  !< (ndime,npoin)       To be used in FSI problems
  use def_master, only: gdepo !< (ndime,ndime,npoin) Deformation gradient
  use def_solidz, only: kfl_fixno_sld
  implicit none
  real(rp), intent(out) :: rhs_id(:)
  real(rp)    :: foref(3)
  integer(ip) :: idofn, idime, ipoin
  !
  !< call sld_coupli(ITASK_MATRIX)
  if( INOTMASTER ) then
    !
    do ipoin = 1,npoin
      !
      do idime = 1,ndime
        foref(idime) = dot_product( gdepo(idime,1:ndime,ipoin), forc_fsi(1:ndime,ipoin) )
      enddo
      !
      idofn = (ipoin-1) * ndime
      do idime = 1,ndime
        idofn = idofn + 1
        if( kfl_fixno_sld(idime,ipoin) /= 1 ) rhs_id(idofn) = rhs_id(idofn) + foref(idime)
      end do
    end do
    !
  endif
  !
  !-----------------------------------------------------------------------||---!
  if(.not.INOTMASTER) then
    !print *, "commdom_alya_sld_matrix"
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

!===============================================================| contains |===!

end subroutine sld_coupli
