subroutine nsi_elmrhs(&
     pgaus,plapl,gpden,gpcod,gpvel,gppre,gptem,gpgve,&
     gpgv2,gpsgs,gpgrk,gpvis,gpgrt,gpgde,gpdiv,gprhs)
  !----------------------------------------------------------------------
  !****f* Nastin/nsi_elmrhs
  ! NAME 
  !    nsi_elmrhs
  ! DESCRIPTION
  !    This routine computes the right-hand side of the momentum 
  !    equations, including:
  !    The total force term:
  !    - Gravity ............. rho*g 
  !    - Frame accelration ... rho*[ - w x [w x (r-r0)] + dw/dt x r + a ]
  !    - Time derivative ..... rho/(theta*dt)*u^{n-1}
  !    - Boussinesq approx ... -rho*beta*|gb|*g*(T-Tr)
  !    - Turbulence .......... -2/3*rho*grad(K)
  !    Linearization term:
  !    - Newton-Raphson ...... rho*(u.grad)u
  !    Time derivative:
  !    - Trapezoidal ......... rho/(theta*dt)*u^{n-1}
  !    - BDF ................. sumk rho/dt*Bk*u^{k}
  ! INPUT
  !    PABDF_NSI ............. BDF Paremeters Bk
  !    BOUGR_NSI ............. Boussinesq gravity norm |gb|
  !    BOUBE_NSI ............. Volume expansion beta
  !    BOUTR_NSI ............. Boussinesq reference temperature Tr
  !    DTINV_NSI ............. Time step inverse 1/dt
  !    GPDEN ................. Density rho
  !    GPCOD ................. Coordinate r
  !    GPVEL ................. Velocity u
  !    GPTEM ................. Temperature T
  !    GPGVE ................. Velocity gradients grad(u)
  !    GRNOR_NSI ............. Gravity norm |g|
  !    GRAVI_NSI ............. Gravity vector g
  !    FACCA_NSI ............. Angular accelaration dw/dt
  !    FACCL_NSI ............. Linear acceleration a
  !    FVELA_NSI ............. Angular velocity w
  !    FROTC_NSI ............. Center of rotation r0
  !    GPSGS ................. Subgrid scale u'
  !    GPGRK ................. Turb. kinetic energy gradients grad(K)
  ! OUTPUT
  !    GPRHS ................. Right hand side
  ! USES
  ! USED BY
  !    nsi_elmope 
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  dpthe,prthe
  use def_domain, only     :  ndime,mgaus,ntens
  use def_nastin, only     :  kfl_timei_nsi,kfl_cotem_nsi,kfl_linea_nsi,&
       &                      kfl_sgsti_nsi,kfl_advec_nsi,kfl_exacs_nsi,&
       &                      kfl_tiacc_nsi,kfl_tisch_nsi,kfl_grtur_nsi,&
       &                      bougr_nsi,boube_nsi,boutr_nsi,dtinv_nsi,&
       &                      grnor_nsi,gravi_nsi,facca_nsi,faccl_nsi,&
       &                      fvela_nsi,frotc_nsi,pabdf_nsi,pcoef_nsi,&
       &                      kfl_regim_nsi,kfl_visco_nsi,kfl_grvis_nsi,&
       &                      dtsgs_nsi,nbdfp_nsi, gravb_nsi
  use def_kermod, only       :  gasco
  implicit none
  integer(ip), intent(in)  :: pgaus,plapl
  real(rp),    intent(in)  :: gpden(pgaus,2),gpgrk(ndime,pgaus)
  real(rp),    intent(in)  :: gpcod(ndime,pgaus),gpvel(ndime,pgaus,*)
  real(rp),    intent(in)  :: gppre(pgaus,*),gptem(pgaus,*)
  real(rp),    intent(in)  :: gpgve(ndime,ndime,pgaus)
  real(rp),    intent(in)  :: gpgv2(ntens,ndime,pgaus)
  real(rp),    intent(in)  :: gpvis(pgaus)
  real(rp),    intent(in)  :: gpsgs(ndime,pgaus)
  real(rp),    intent(in)  :: gpgrt(ndime,pgaus),gpgde(ndime,pgaus)
  real(rp),    intent(in)  :: gpdiv(pgaus)
  real(rp),    intent(out) :: gprhs(ndime+1,pgaus)
  real(rp)                 :: alpha(3),centf(3),gppos(3)
  real(rp)                 :: fact1,fact2,factt
  integer(ip)              :: igaus,idime,jdime,itime
  !
  ! Initialization
  ! 
  do igaus=1,pgaus
     do idime=1,ndime+1
        gprhs(idime,igaus)=0.0_rp
     end do
  end do
  !----------------------------------------------------------------------
  !
  ! MOMENTUM EQUATION
  !
  !----------------------------------------------------------------------
  !
  ! External forces: rho*[ g - w x (w x r) + dw/dt x r + a 
  !
  if(kfl_exacs_nsi==0) then
     do igaus=1,pgaus
        gppos(3)=0.0_rp
        do idime=1,ndime
           gppos(idime)=&
                gpcod(idime,igaus)-frotc_nsi(idime) ! Relative position r-r0
        end do
        call vecpro(facca_nsi,gppos,alpha,3)        ! Angular acceleration dw/dt x (r-r0)
        call vecpro(fvela_nsi,gppos,centf,3)        ! Centrifugal force w x [w x (r-r0)]
        call vecpro(fvela_nsi,centf,centf,3)
        do idime=1,ndime
           gprhs(idime,igaus)=gpden(igaus,1)*(&
                grnor_nsi*gravi_nsi(idime)&         ! g
                -centf(idime)&                      ! w x (w x r)
                -alpha(idime)&                      ! dw/dt x r
                -faccl_nsi(idime))                  ! a
        end do
     end do
  end if
  !
  ! Time derivative: rho/(theta*dt)*u^{n-1}
  !
  if(kfl_timei_nsi==1) then
     do igaus=1,pgaus
        fact1=gpden(igaus,1)*dtinv_nsi
        do itime=2,nbdfp_nsi
           factt=fact1*pabdf_nsi(itime)
           do idime=1,ndime
              gprhs(idime,igaus)=gprhs(idime,igaus)&
                   -factt*gpvel(idime,igaus,itime)
           end do
        end do
     end do
  end if
  !
  ! Boussinesq coupling: -rho*beta*g*(T-Tr)
  !
  if(kfl_cotem_nsi==1) then
     do igaus=1,pgaus
        fact2=gpden(igaus,1)*bougr_nsi*boube_nsi
        do idime=1,ndime
           gprhs(idime,igaus)=gprhs(idime,igaus)&
                -fact2*(gptem(igaus,1)-boutr_nsi)*gravb_nsi(idime)
        end do
     end do
  end if
  !
  ! Newton-Raphson
  !
  if(kfl_advec_nsi==1.and.kfl_linea_nsi==2) then
     do igaus=1,pgaus
        do idime=1,ndime
           do jdime=1,ndime
              gprhs(idime,igaus)=gprhs(idime,igaus)&
                   +gpden(igaus,1)*gpvel(jdime,igaus,1)&
                   *gpgve(jdime,idime,igaus)
           end do
        end do
     end do
  end if
  !
  ! Turbulence kinetic energy: -2/3*rho*K
  !
  if( kfl_grtur_nsi /= 0 ) then
     do igaus = 1,pgaus
        fact1 = gpden(igaus,1)*2.0_rp/3.0_rp
        do idime = 1,ndime
           gprhs(idime,igaus) = gprhs(idime,igaus)&
                - fact1 * gpgrk(idime,pgaus)
        end do
     end do
  end if
  !
  ! Subgrid scale time derivative: rho/(theta*dt)*u^{n-1}
  !
  if( kfl_timei_nsi == 1 .and. kfl_sgsti_nsi == 1 ) then
     do igaus = 1,pgaus
        fact1 = gpden(igaus,1)*dtsgs_nsi
        do idime = 1,ndime
           gprhs(idime,igaus) = gprhs(idime,igaus)&
                + fact1 * gpsgs(idime,igaus)
        end do
     end do
  end if

  !----------------------------------------------------------------------
  !
  ! CONTINUITY EQUATION
  !
  !----------------------------------------------------------------------
  !
  ! Time derivative: (1-kappa)/(p*theta*dt)*p^{n-1}
  !
  if(kfl_timei_nsi==1) then

     if(kfl_regim_nsi==111) then
        if(kfl_tisch_nsi==1) then
           !
           ! Trapezoidal rule: OJO GPDEN NO SE HA CALCULADO
           !
           do igaus=1,pgaus
              gprhs(ndime+1,igaus)=gprhs(ndime+1,igaus)&
                   -dtinv_nsi*pabdf_nsi(2)*gpden(igaus,2)
           end do
        else if(kfl_tisch_nsi==2) then
           !
           ! BDF scheme
           !
           do igaus=1,pgaus
              fact1=dtinv_nsi
              do itime=2,kfl_tiacc_nsi+1
                 factt=fact1*pabdf_nsi(itime)
                 gprhs(ndime+1,igaus)=gprhs(ndime+1,igaus)&
                      -factt*gppre(igaus,itime)
              end do
           end do
        end if
     end if
  end if

  if(kfl_regim_nsi==1) then
     !
     ! Compressible flow: (rho/T)*DT/Dt + (rho/p)*p^n/dt
     !
     do igaus=1,pgaus
        fact1=gpden(igaus,1)/gptem(igaus,1)
        do idime=1,ndime              
           gprhs(ndime+1,igaus)=gprhs(ndime+1,igaus)&
                +fact1*gpvel(idime,igaus,1)*gpgrt(idime,igaus)  ! (rho/T)*[u.grad(T)]
        end do
     end do
     if(kfl_timei_nsi==1) then
        do igaus=1,pgaus
           gprhs(ndime+1,igaus)=gprhs(ndime+1,igaus)&
                +dtinv_nsi*pabdf_nsi(1)/gptem(igaus,1)*(gpden(igaus,1)&
                *(gptem(igaus,1)-gptem(igaus,2))&               ! +(rho/T)*d(rho)/dt
                +gppre(igaus,2)/gasco)                      ! +p^n/(dt*RT)
        end do        
     end if

  else if(kfl_regim_nsi==2) then
     !
     ! Compressible flow: rho*div(u) + rho^n/dt
     !
     do igaus=1,pgaus
        gprhs(ndime+1,igaus)=gprhs(ndime+1,igaus)&
             +gpden(igaus,1)*gpdiv(igaus)
     end do
     if(kfl_timei_nsi==1) then
        do igaus=1,pgaus
           gprhs(ndime+1,igaus)=gprhs(ndime+1,igaus)&
                +dtinv_nsi*pabdf_nsi(1)*gpden(igaus,2)
        end do
     end if

  else if(kfl_regim_nsi==3) then
     !
     ! Low-Mach: (rho/T)*DT/Dt - (rho/p0)*Dp0/dt
     !
     do igaus=1,pgaus
        fact1=gpden(igaus,1)/gptem(igaus,1)
        gprhs(ndime+1,igaus)=gprhs(ndime+1,igaus)&     ! -(rho/p0)*Dp0/dt
             -gpden(igaus,1)/prthe(1)*dpthe
        do idime=1,ndime                               ! rho/T*[u.grad(T)]
           gprhs(ndime+1,igaus)=gprhs(ndime+1,igaus)&
                +fact1*gpvel(idime,igaus,1)*gpgrt(idime,igaus)
        end do
     end do
     if(kfl_timei_nsi==1) then                         ! rho/T*dT/dt
        do igaus=1,pgaus
           gprhs(ndime+1,igaus)=gprhs(ndime+1,igaus)&
                +dtinv_nsi*pabdf_nsi(1)*gpden(igaus,1)/gptem(igaus,1)&
                *(gptem(igaus,1)-gptem(igaus,2))
        end do        
     end if
  end if

end subroutine nsi_elmrhs
