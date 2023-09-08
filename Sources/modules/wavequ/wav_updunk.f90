subroutine wav_updunk(itask)
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_updunk
  ! NAME 
  !    wav_updunk
  ! DESCRIPTION
  !    This routine performs several types of updates for the wave
  !    amplitude
  ! USED BY
  !    wav_begste (itask=1)
  !    wav_begite (itask=2)
  !    wav_endite (itask=3, inner loop) 
  !    wav_endite (itask=4, outer loop) 
  !    wav_endste (itask=5)
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_wavequ
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: itime,icomp,ipoin,pelty,pgaus,ielem,igaus
  real(rp)                :: dtn,dtn2

  if(kfl_paral/=0) then

     select case (itask)

     case(1)
        !
        ! Assign u(n,0,*) <-- u(n-1,*,*), initial guess for outer iterations
        !
        do ipoin=1,npoin
           wavam(ipoin,2) = wavam(ipoin,min(3,ncomp_wav))
        end do

     case(2)
        !
        ! Assign u(n,i,0) <-- u(n,i-1,*), initial guess for inner iterations
        !
        do ipoin=1,npoin
           wavam(ipoin,1) = wavam(ipoin,2)
        end do
        if(kfl_timet_wav==2) then
           do ipoin=1,npoin
              unkno(ipoin) = wavam(ipoin,2)
           end do
        end if

     case(3)
        !
        ! Assign u(n,i,j-1) <-- u(n,i,j), update of the wave amplitude
        !
        do ipoin=1,npoin
           wavam(ipoin,1) = unkno(ipoin)
        end do

     case(4)
        !
        ! Assign u(n,i-1,*) <-- u(n,i,*)
        !        
        do ipoin=1,npoin
           wavam(ipoin,2) = wavam(ipoin,1)
        end do

     case(5)
        !
        ! Update previous time steps
        !         
        if(kfl_tisch_wav==2) then
           !
           ! Leap-frog scheme
           !
           do ipoin=1,npoin
              wavam(ipoin,4)=wavam(ipoin,3) ! u^{n-2}
              wavam(ipoin,3)=wavam(ipoin,1) ! u^{n-1}
           end do
           if(kfl_subgs_wav==1) then
              !
              ! Subgrid scale
              !
              do ielem=1,nelem
                 pelty=ltype(ielem)
                 pgaus=ngaus(pelty)
                 do igaus=1,pgaus
                    wasgs(ielem)%a(igaus,3)=wasgs(ielem)%a(igaus,2)
                    wasgs(ielem)%a(igaus,2)=wasgs(ielem)%a(igaus,1)
                 end do
              end do
           end if

        else if(kfl_tisch_wav==3) then
           !
           ! Newmark
           !
           dtn   = 1.0_rp/dtinv
           dtn2  = dtn*dtn
           do ipoin=1,npoin
              wavac_wav(ipoin,1) = (wavam(ipoin,1)-wavam(ipoin,3))/(nebet_wav*dtn2)&
                   &               -wavve_wav(ipoin,2)/(nebet_wav*dtn)&
                   &               -(0.5_rp/nebet_wav-1.0_rp)*wavac_wav(ipoin,2)
              wavve_wav(ipoin,1) = wavve_wav(ipoin,2)+((1.0_rp-negam_wav)*wavac_wav(ipoin,2)&
                   &               +negam_wav*wavac_wav(ipoin,1))*dtn
              wavam(ipoin,3)     = wavam(ipoin,1)
              wavac_wav(ipoin,2) = wavac_wav(ipoin,1)
              wavve_wav(ipoin,2) = wavve_wav(ipoin,1)
           end do

        end if

     case(6) 
        !
        ! Assign u(n,1,*) <-- u(n-1,*,*), when reading from restart file
        ! 
        icomp=min(3,ncomp_wav) 
        do ipoin=1,npoin
           wavam(ipoin,1) = wavam(ipoin,icomp)
        end do

     end select

  end if

end subroutine wav_updunk

