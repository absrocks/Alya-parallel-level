subroutine chm_accumu()
  !-----------------------------------------------------------------------
  !****f* chemic/chm_accumu
  ! NAME 
  !    chm_accumu
  ! DESCRIPTION
  !    This routine computes the ground load accumulation (mass per unit area)
  !    Ground (surface of accumulation) is defined by kfl_fixno_chm = 2
  ! USES
  !
  ! USED BY
  !    chm_endste
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_domain
  use def_master
  use def_chemic
  implicit none
!  
  integer(ip) :: ipoin,ibopo,iclas
  real(rp)    :: xfact
!

  if( INOTMASTER ) then 
       !
       !  Computes pts_accum
       !
       do ipoin=1,npoin
          if(abs(kfl_fixno_chm(1,ipoin))==2) then
              !
              ! This is a ground node (KFL_FIXNO_PTS=2). Condition
              ! does not depend on classes
              ! 
              ibopo = lpoty(ipoin) 
!
              do iclas=1,nclas_chm
!
!                Conce is dimensionless, so it is necessary to multipy
!                by the density of the media
!
                 xfact = dtime*densi_chm(ipoin)*exnor(ndime,1,ibopo)
!
!                Compute accumulation at the time mid-point
!
                 accum_chm(ipoin) = accum_chm(ipoin) &
                      -xfact*0.5_rp*(conce(ipoin,iclas,1)+conce(ipoin,iclas,ncomp_chm))*vterm_chm(ipoin,iclas)
!
              end do
          else
              !
              ! Other nodes
              !
              accum_chm(ipoin) = 0.0_rp
         end if
      end do
!
  end if
!
end subroutine chm_accumu
