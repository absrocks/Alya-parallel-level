subroutine tem_matrix()
  !------------------------------------------------------------------------
  !
  ! Compute elemental matrix and RHS of the following problem:
  !   _                             _
  !  |                             |
  !  |  rho*cp/(theta*dt) T*v dw + |  rho*cp*[u.grad(T)] dw
  ! _|W                           _|W
  !
  !    _                         _               _
  !   |                         |               |
  ! + |  k*grad(T).grad(v) dw + |  ar*Tr*v ds = |  Q*v dw
  !  _|W                       _|S             _|
  !
  !   _                               _
  !  |                               |
  !  |  rho*cp/(theta*dt) T^n*v dw - |  (qr-ar*Tr) v ds
  ! _|W                             _|S
  ! 
  ! All terms are calculated at current time step, i.e. n+theta
  !
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_temper
  use def_domain
  use def_coupli
  use def_kermod,         only :  wallcoupling_waexl, kfl_waexl_ker, temel_ker
  use mod_interpolation
  use mod_couplings
  use mod_communications, only : PAR_MIN
  use mod_ker_timeline,   only : ker_timeline
  use mod_messages,       only : messages_live
  use mod_timings,        only : timings_assembly
  implicit none
  integer(ip) :: ipoin  
  real(rp)    :: time1,time2,time3,time_element,time_boundary
  
  call ker_timeline('INI_ASSEMBLY')
  !
  ! Initializations
  !
  call inisol()
  resgs_tem(1) = 0.0_rp
  resgs_tem(2) = 0.0_rp

  call cputim(time1) 

  if( kfl_waexl_ker == 1_ip ) &
       call COU_GET_INTERPOLATE_POINTS_VALUES(tempe,temel_ker,wallcoupling_waexl)
  
  if( INOTMASTER ) then 

     if( kfl_discr_tem == NODAL_SCHEME ) then
        !
        ! FE: Element assembly
        !        
        call tem_elmope_new(1_ip)

     else if( kfl_discr_tem == CELL_CENTERED_SCHEME ) then
        !
        ! FV: Assembly
        !
        call tem_finite_volume()
     end if

  end if
  call cputim(time2)

  if( INOTMASTER ) then 

     if( kfl_discr_tem == NODAL_SCHEME ) then
        !
        ! FE: Boundary assembly
        !
        call tem_bouope()
     end if
  end if
  call cputim(time3)

  if( INOTMASTER ) then 
     !
     ! Impose non-uniform heat flux (from FIELDS)
     !    
     if( kfl_flux_tem > 0_ip ) then
        do ipoin = 1,npoin
           if( ipoin <= npoin_own ) then
              rhsid(ipoin) = rhsid(ipoin) + heat_flux(1,ipoin)
           end if
        end do
     end if
  end if
  !
  ! Timings
  !
  time_element  = time2 - time1
  time_boundary = time3 - time2
  call timings_assembly(time_element,time_boundary)  

  if( kfl_reset /= -1 ) then
     !
     ! Look for minimum over subdomains
     ! 
     call PAR_MIN(dtcri_tem,'IN MY ZONE')
     !
     !
     if (dtcri_tem /= 0.0_rp) dtinv_tem = 1.0_rp/(dtcri_tem*safet_tem)
     !
     ! If changing dtinv is necessary, activate reset mode
     if (dtinv_tem > reset_factor * dtinv) then
        dtinv     = dtinv_tem
        kfl_reset = 1
        if( INOTSLAVE ) call messages_live('REQUESTED RESET OF TIME STEP')
     endif
  end if

  call ker_timeline('END_ASSEMBLY')
  
end subroutine tem_matrix

subroutine tem_test()
  use def_master
  use mod_parall
  use def_domain
  use def_coupli
  use mod_couplings
  implicit none

  integer(ip)         :: kcoup,icoup,ndofn,kpoin,ipoin
  real(rp),   pointer :: xx_send(:)
  real(rp),   pointer :: xx_recv(:)

  nullify(xx_send,xx_recv)
  
  if( INOTMASTER ) then
     allocate( xx_send(npoin) )
     allocate( xx_recv(npoin) )
     xx_send = 1.0_rp
     xx_recv = 0.0_rp
     ndofn    = 1
     do kcoup = 1,ncoup_implicit_n
        icoup = lcoup_implicit_n(kcoup)
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndofn,xx_recv,xx_send)
        color_target = coupling_type(icoup) % color_target
        if(I_AM_IN_COLOR(color_target)) then
           do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
              ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
              write(*,'(a,i1,1x,i3,100(1x,e12.6))')'N: ',kfl_paral,lninv_loc(ipoin),xx_recv(ipoin)
           end do
        end if
     end do    
     deallocate( xx_send )
     deallocate( xx_recv )
  end if
    
  if( INOTMASTER ) then
     allocate( xx_send(npoin) )
     allocate( xx_recv(npoin) )
     xx_send = 1.0_rp
     xx_recv = 0.0_rp
     ndofn    = 1
     do kcoup = 1,ncoup_implicit_d
        icoup = lcoup_implicit_d(kcoup)
        call COU_INTERPOLATE_NODAL_VALUES(icoup,ndofn,xx_recv,xx_send)
        color_target = coupling_type(icoup) % color_target
        if(I_AM_IN_COLOR(color_target)) then
           do kpoin = 1,coupling_type(icoup) % wet % npoin_wet
              ipoin = coupling_type(icoup) % wet % lpoin_wet(kpoin)
              write(*,'(a,i1,1x,i3,100(1x,e12.6))')'D: ',kfl_paral,lninv_loc(ipoin),xx_recv(ipoin)
           end do
        end if
     end do    
     deallocate( xx_send )
     deallocate( xx_recv )
  end if
  
  
end subroutine tem_test


!!$subroutine youssef()
!!$  
!!$  use def_parame
!!$  use def_master
!!$  use def_temper
!!$  use def_domain
!!$  use def_coupli
!!$  use mod_interpolation
!!$  use mod_couplings
!!$  use mod_gradie
!!$  use mod_postpr
!!$  use mod_maillage
!!$  
!!$  implicit none
!!$  
!!$  integer(ip)       :: ipoin,idime,ii,jj
!!$
!!$  real(rp), pointer :: hessi(:,:,:)
!!$  real(rp), pointer :: metric(:,:,:)
!!$  
!!$  real(rp)          :: eigen_value(ndime)
!!$  real(rp)          :: eigen_vector(ndime)
!!$  real(rp)          :: x,y,c(5),eps
!!$
!!$
!!$  do jj = 1,10
!!$  nullify(hessi)
!!$  nullify(metric)
!!$  
!!$  do ipoin = 1,npoin
!!$     tempe(ipoin,1) = 0.25_rp*( (coord(1,ipoin)-0.5_rp)**2 + (coord(2,ipoin)-0.5_rp)**2 )
!!$  end do
!!$
!!$  tempe = 0.0_rp
!!$  c = (/ 0.0_rp,-0.6_rp,0.6_rp,-1.2_rp,1.2_rp /)
!!$  eps = 0.01_rp
!!$  do ipoin = 1,npoin
!!$     x = coord(1,ipoin)
!!$     y = coord(2,ipoin)
!!$     do ii = 1,5
!!$        tempe(ipoin,1) = tempe(ipoin,1) + 1.0_rp / ( 1.0_rp+exp( (x+y-c(ii))/(2.0_rp*eps) ) ) + 1.0_rp / ( 1.0_rp+exp( (x-y-c(ii))/(2.0_rp*eps) ) ) 
!!$     end do
!!$  end do
!!$  
!!$  call projec_hessian(tempe(:,1),hessi)
!!$  !
!!$  !allocate(hessi(2,2,npoin)) 
!!$  ! Add small value to Hessian
!!$  !
!!$  do ipoin = 1,npoin
!!$     do idime = 1,ndime
!!$        if( hessi(idime,idime,ipoin) == 0.0_rp ) & 
!!$             hessi(idime,idime,ipoin) = hessi(idime,idime,ipoin) + zeror
!!$     end do
!!$  end do
!!$
!!$  call maillage_invocation(hessi)
!!$  deallocate(hessi)
!!$end do
!!$  !
!!$  ! Eigenvalues and eigenvectors
!!$  !
!!$  do ipoin = 1,npoin
!!$     call eigen_2by2(hessi(1,1,ipoin),eigen_vector,eigen_value)
!!$  end do
!!$  
!!$  stop
!!$  
!!$end subroutine youssef
!!$!
!!$! http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/
!!$!
!!$subroutine eigen_2by2(M,E,V)
!!$  use def_kintyp
!!$  implicit none
!!$  real(rp), intent(in)  :: M(2,2)
!!$  real(rp), intent(out) :: E(2,2)
!!$  real(rp), intent(out) :: V(2)
!!$  real(rp)              :: a,b,c,d,T,Dd,L1,L2
!!$
!!$  a  = M(1,1)
!!$  b  = M(1,2)
!!$  c  = M(2,1)
!!$  d  = M(2,2)
!!$  
!!$  T  = a+d
!!$  Dd = a*d-b*c
!!$  L1 = 0.5_rp*T + sqrt(0.25_rp*T**2-Dd)
!!$  L2 = 0.5_rp*T - sqrt(0.25_rp*T**2-Dd)
!!$ 
!!$  if( c /= 0.0_rp ) then
!!$     E(1,1) = L1-d
!!$     E(2,1) = c
!!$     E(1,2) = L2-d
!!$     E(2,2) = c
!!$  else if( b /= 0.0_rp ) then
!!$     E(1,1) = b
!!$     E(2,1) = L1-a
!!$     E(1,2) = b
!!$     E(2,2) = L2-a
!!$  else
!!$     E(1,1) = 1.0_rp
!!$     E(2,1) = 0.0_rp
!!$     E(1,2) = 0.0_rp
!!$     E(2,2) = 1.0_rp   
!!$  end if
!!$
!!$  V(1) = L1
!!$  V(2) = L2
!!$
!!$end subroutine eigen_2by2

