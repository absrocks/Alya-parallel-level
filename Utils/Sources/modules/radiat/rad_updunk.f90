subroutine rad_updunk(itask)
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_updunk
  ! NAME 
  !    rad_updunk
  ! DESCRIPTION
  !    This routine performs several types of updates for the radiation
  ! USED BY
  !    rad_begste (itask=1)
  !    rad_begite (itask=2)
  !    rad_endite (itask=3, inner loop) 
  !    rad_endite (itask=4, outer loop) 
  !    rad_endste (itask=5)
  !    rad_restar (itask=6)
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  use def_kermod
  use mod_ker_proper 
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ielem,igaus,itime,pgaus,ipoin,icomp,ispec,dummi
  real(rp)                :: rela1_rad,absor_xx,densi_rad(npoin),dummr(1)

  if( INOTMASTER ) then

     select case (itask)

     case(1_ip)
        !
        ! Assign G(n,0,*) <-- G(n-1,*,*), initial guess for outer iterations
        !
        do ipoin=1,npoin
           radav_rad(ipoin,2) = radav_rad(ipoin,ncomp_rad)
        end do

     case(2_ip)
        !
        ! Assign G(n,i,0) <-- G(n,i-1,*), initial guess for inner iterations
        !
        do ipoin=1,npoin
           radav_rad(ipoin,1) = radav_rad(ipoin,2)
           unkno(ipoin)   = radav_rad(ipoin,2)
        end do

     case(3_ip)
        !
        ! Assign G(n,i,j-1) <-- G(n,i,j), update of the radiation
        !
        if(postp(1)%npp_stepi(8)>0) then
           do ipoin=1,npoin
              raold_rad(ipoin) = radav_rad(ipoin,1)
           end do           
        end if
        do ipoin=1,npoin
           radav_rad(ipoin,1) = unkno(ipoin)
        end do

     case(4_ip)
        !
        ! Assign G(n,i-1,*) <-- G(n,i,*)
        !        
        do ipoin=1,npoin
           radav_rad(ipoin,2) = radav_rad(ipoin,1)
        end do

     case(5_ip)
!!$        !
!!$        ! Obtain G(n,*,*) for the Crank-Nicolson method and assign
!!$        ! G(n-1,*,*) <-- G(n,*,*)
!!$        !        
!!$        if(kfl_tisch_rad==1.and.kfl_tiacc_rad==2) then
!!$           !
!!$           ! Crank-Nicolson method 
!!$           !        
!!$           do ipoin=1,npoin
!!$              radav_rad(ipoin,1) = 2.0_rp*radav_rad(ipoin,1)-radav_rad(ipoin,3)
!!$           end do
!!$
!!$        else if(kfl_tisch_rad==2) then
!!$           !
!!$           ! BDF scheme
!!$           !
!!$           do ipoin=1,npoin
!!$              do itime=2+kfl_tiacc_rad,4,-1
!!$                 radav_rad(ipoin,itime) = radav_rad(ipoin,itime-1)
!!$              end do
!!$           end do
!!$        end if
!!$        call rad_averag()
        do ipoin=1,npoin
           radav_rad(ipoin,3) = radav_rad(ipoin,1)
        end do

!!$        if(kfl_sgsti_rad==1) then
!!$           !
!!$           ! Time tracking of the subscales
!!$           !        
!!$           if(kfl_tiacc_rad==2) then 
!!$              do ielem=1,nelem
!!$                 pgaus=ngaus(ltype(ielem))
!!$                 do igaus=1,pgaus
!!$                    rasgs_rad(ielem)%a(igaus,1)=2.0_rp*rasgs_rad(ielem)%a(igaus,1)&
!!$                         -rasgs_rad(ielem)%a(igaus,2)               
!!$                 end do
!!$              end do
!!$           end if
!!$
!!$           do ielem=1,nelem
!!$              pgaus=ngaus(ltype(ielem))
!!$              do igaus=1,pgaus
!!$                 rasgs_rad(ielem)%a(igaus,2)=rasgs_rad(ielem)%a(igaus,1)
!!$              end do
!!$           end do
!!$
!!$        end if

     case(6_ip) 
        !
        ! Assign G(n,1,*) <-- G(n-1,*,*), when readin from restart file
        ! 
        icomp=min(3_ip,ncomp_rad) 
        do ipoin=1,npoin
           radav_rad(ipoin,1) = radav_rad(ipoin,icomp)
        end do

     case(7_ip)
        !
        ! Update heat radiation source
        !
        call ker_proper('DENSI','NPOIN',dummi,dummi,densi_rad)
        do icomp=1,ncomp_rad
           do ipoin=1,npoin
              absor_xx = 0.0_rp
              do ispec=1,nspec_rad
                 absor_xx = absor_xx + absor_rad(ispec)* densi_rad(ipoin)*conce_rad(ipoin,ispec)
              enddo
              radso(ipoin,icomp)= 4.0_rp * absor_xx* steph_rad * tempe_rad(ipoin)**4 &
                   - absor_xx * radav_rad(ipoin,icomp)
           end do
        enddo
        
     end select

  end if

end subroutine rad_updunk

