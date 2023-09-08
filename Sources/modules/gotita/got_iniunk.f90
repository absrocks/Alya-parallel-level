subroutine got_iniunk
!-----------------------------------------------------------------------
!****f* Gotita/got_iniunk
! NAME 
!    got_iniunk
! DESCRIPTION
!    This routine sets up the initial condition for the droplet velocity.
!    If this is a restart, initial condition are loaded somewhere else.
! USED BY
!    got_begste
!***
!-----------------------------------------------------------------------
  use      def_parame 
  use      def_master
  use      def_domain
  use      def_gotita
  use      mod_memchk
  implicit none
  integer(ip) :: ipoin,idime,icomp
  !
  ! Problem is transient
  !
  kfl_stead_got = 0
  !
  ! Total iterations counter then
  !
  if(kfl_rstar==0) ittot_got=0
  !
  ! Initial value of droplet velocity
  !
  if(kfl_rstar==0) then   

     if(kfl_paral/=0) then
        !
        ! Actualize initial values
        !        
        icomp=min(3,ncomp_got) 
        if(kfl_probl_got==1) then
           do ipoin=1,npoin
              do idime=1,ndime
                 vdrop(idime,ipoin,icomp)=bvess_got(idime,ipoin)
              end do
              cdrop(ipoin,icomp)=bvess_got(ndofn_got(3),ipoin)
           end do
        else if(kfl_probl_got==2) then
           do ipoin=1,npoin
              do idime=1,ndime
                 vdrop(idime,ipoin,icomp)=bvess_got(idime,ipoin)
              end do
           end do
        else if(kfl_probl_got==3) then
           do ipoin=1,npoin
              cdrop(ipoin,icomp)=bvess_got(1,ipoin)
           end do
        end if
        !
        ! Actualize current iteration VELOC(:,:,1) <= VELOC(:,:,icomp)
        !
        call got_updunk(11_ip)
     end if

  else
     !
     ! Restart
     !
     call got_restar(1_ip)

  end if

end subroutine got_iniunk
