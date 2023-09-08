subroutine chm_upwmea(itask) 
  !-----------------------------------------------------------------------
  ! Update
  ! NAME 
  !    chm_upwmea
  ! DESCRIPTION
  !    Update each of the species specific heat and viscosity
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,npoin,mnode,coord
  use def_chemic, only      :  nspec_chm,kfl_warni_chm,nclas_chm,kfl_model_chm,kfl_tiaor_chm, &
                               kfl_tisch_chm
  use def_master, only      :  conce,wmean,inotmaster,speci
  implicit none
  integer(ip),intent(in)    :: itask   ! 1 = init step, 2=local update, 3=final update
  integer(ip)               :: ipoin,ispec,iclas,itime

  if (INOTMASTER) then
     select case(itask)
     case (1_ip) 
        !
        ! Compute mean molecular weight for time = 0
        !
        do ipoin = 1,npoin
           wmean(ipoin,3) = 0.0_rp
           do ispec = 1,nspec_chm
              wmean(ipoin,3) = wmean(ipoin,3) +  conce(ipoin,ispec,3) / speci(ispec)%weigh 
           enddo
           if(wmean(ipoin,3) /= 0.0_rp) then 
              wmean(ipoin,3) = 1.0_rp / wmean(ipoin,3)
           else
              if (kfl_warni_chm==1) print *, ' WARNING: INITIAL point with zero mean weight',ipoin
           endif
           wmean(ipoin,2) = wmean(ipoin,3)
           wmean(ipoin,1) = wmean(ipoin,2)
        enddo

     case (2_ip) ! At Begste
        do ipoin = 1,npoin
           wmean(ipoin,2) = wmean(ipoin,3)
        enddo
     case (3_ip) ! At Begite
        do ipoin = 1,npoin
           wmean(ipoin,1) = wmean(ipoin,2)
        enddo
     case (4_ip) ! At Endite(2)
        do ipoin = 1,npoin
           wmean(ipoin,2) = wmean(ipoin,1)
        enddo
     case (5_ip) ! At Endste
        if(kfl_tisch_chm==2) then
           !
           ! BDF scheme
           !
           do ipoin=1,npoin
              do itime=2+kfl_tiaor_chm,4,-1
                 wmean(ipoin,itime) = wmean(ipoin,itime-1)
              end do
           end do
        end if
        do ipoin = 1,npoin
           wmean(ipoin,3) = wmean(ipoin,1)
        enddo
     case (6_ip) ! Update mean molecular weight, usually at endite(1)

        do ipoin = 1,npoin
           wmean(ipoin,1) = 0.0_rp
           do ispec = 1,nspec_chm
              wmean(ipoin,1) = wmean(ipoin,1) +  conce(ipoin,ispec,1) / speci(ispec)%weigh 
           enddo
           if (wmean(ipoin,1) /= 0.0_rp) then
              wmean(ipoin,1) = 1.0_rp / wmean(ipoin,1)
           else
              if (kfl_warni_chm==1) print *, ' WARNING: point with zero mean weight',ipoin
           endif
        enddo

     case (7_ip) ! For restart: wmean(:,3) = wmean(:,1)
        do ipoin = 1,npoin
           wmean(ipoin,3) = wmean(ipoin,1)
        enddo

     end select
  endif

end subroutine chm_upwmea
