subroutine elmcel(&
     pnode,pgaus,ptopo,plapl,pelty,ndofn,kfl_shock,kfl_taust,&
     staco,gpdif,gprea,gpden,gpvel,gpvol,gpgrd,gprhs,&
     gpsha,gpcar,gphes,elunk,chale, tau2, tau1,dtinv,&
     shock,gplap,elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* mathru/elmcel
  ! NAME
  !   elmcel
  ! DESCRIPTION
  !   General assembly for ADR equation. The units are:
  !   GPSTT ... [1]
  !   GPSTP ... L/(rho*U)
  !   
  ! OUTPUT
  ! USES
  ! USED BY
  !   *_elmope
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,ntens,mnode
  use def_elmtyp, only       :  TRI03
  use mod_tauadr, only       :  tauadr
  implicit none

  integer(ip), intent(in)    :: pnode,pgaus,ptopo,plapl,pelty,ndofn
  integer(ip), intent(in)    :: kfl_shock,kfl_taust
  real(rp),    intent(in)    :: staco(3)
  real(rp),    intent(in)    :: gpdif(pgaus)
  real(rp),    intent(in)    :: gprea(pgaus)
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpvel(ndime,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: gpgrd(ndime,pgaus)
  real(rp),    intent(in)    :: gprhs(ndofn,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(in)    :: elunk(pnode)
  real(rp),    intent(in)    :: chale(2)
  real(rp),    intent(inout) :: tau2(pgaus)
  real(rp),    intent(inout) :: tau1(pgaus)
  real(rp),    intent(in)    :: dtinv
  real(rp),    intent(in)    :: shock
  real(rp),    intent(out)   :: gplap(pnode,pgaus)
  real(rp),    intent(out)   :: elmat(pnode,pnode)
  real(rp),    intent(out)   :: elrhs(ndofn,pnode)

  integer(ip)                :: inode,jnode,igaus,idofn,idime
  real(rp)                   :: A1(pgaus) ,A2(pgaus), A3(pgaus), A4(pgaus)
  real(rp)                   :: A5(pgaus), A6(pgaus), A7(pgaus), A8(pgaus)
  real(rp)                   :: A9(pgaus), A10(pgaus)
  real(rp)                   :: A11(pgaus),A12(pgaus),A13(pgaus),A14(pgaus)
  real(rp)                   :: A15(pgaus),A16(pgaus)
  real(rp)                   :: B1(pgaus) ,B2(pgaus) ,B3(pgaus) ,B4(pgaus)
  real(rp)                   :: B5(pgaus) ,B6(pgaus), B7(pgaus)
  real(rp)                   :: B8(pgaus), B9(pgaus)
  real(rp)                   :: D1(pgaus), D2(pgaus), D3(pgaus), D4(pgaus)
  real(rp)                   :: D5(pgaus), CF(pgaus), CD(pgaus)
  real(rp)                   :: gpadv(ndime,pgaus)
  real(rp)                   :: c1 ,c2, c3, c4, c5, c6, c7, c8, c9, c10, c11
  real(rp)                   :: dkdx,dkdy,dkdz,rho
  real(rp)                   :: s,k,e1,e2,e3,ktau2,rtau2,t1st2
  real(rp)                   :: axCF,ayCF,azCF,rhodts,c10CD
  real(rp)                   :: a1x,a1y,a1z,a2x,a2y,a2z,gpnve,x1

  !-------------------------------------------------------------------
  !
  ! Initialization: ELMAT and ELRHS
  !
  !-------------------------------------------------------------------

  do jnode=1,pnode
     do idofn=1,ndofn
        elrhs(idofn,jnode)=0.0_rp
     end do
     do inode=1,pnode
        elmat(inode,jnode)=0.0_rp
     end do
  end do

  !-------------------------------------------------------------------
  !
  ! Redefine advection + compute stabilization
  !
  !-------------------------------------------------------------------

  if(kfl_taust==0_ip.and.kfl_shock/=0_ip) then
     do igaus=1,pgaus
        do idime=1,ndime
           gpadv(idime,igaus)=gpvel(idime,igaus)-gpgrd(idime,igaus)/gpden(igaus)
        end do
        call vecnor(gpadv(1,igaus),ndime,gpnve,2_ip) 
        gpnve=gpden(igaus)*gpnve
        call tauadr(&
             2_ip,staco,gpnve,gpdif(igaus),gprea(igaus),&
             chale(1),chale(2),tau2(igaus))
        tau1(igaus)=1.0_rp
     end do     

  else if(kfl_taust>0_ip) then
     do igaus=1,pgaus
        do idime=1,ndime
           gpadv(idime,igaus)=gpvel(idime,igaus)-gpgrd(idime,igaus)/gpden(igaus)
        end do
        call vecnor(gpadv(1,igaus),ndime,gpnve,2_ip) 
        gpnve=gpden(igaus)*gpnve
        call tauadr(&
             kfl_taust,staco,gpnve,gpdif(igaus),gprea(igaus),&
             chale(1),chale(2),tau2(igaus))
        tau1(igaus)=1.0_rp
     end do

  else
     if(kfl_taust==0_ip) then
        do igaus=1,pgaus
           tau1(igaus)=1.0_rp
           tau2(igaus)=0.0_rp
        end do
     end if
     do igaus=1,pgaus
        do idime=1,ndime
           gpadv(idime,igaus)=gpvel(idime,igaus)-gpgrd(idime,igaus)/gpden(igaus)
        end do
     end do
  end if

  !-------------------------------------------------------------------
  !
  ! Laplacian matrix: GPLAP
  !
  !-------------------------------------------------------------------

  if(plapl/=2) then
     if(plapl==1) then
        if(ndime==1) then
           do igaus=1,pgaus
              do inode=1,pnode
                 gplap(inode,igaus)=&
                      gphes(1,inode,igaus)
              end do
           end do           
        else if(ndime==2) then
           do igaus=1,pgaus
              do inode=1,pnode
                 gplap(inode,igaus)=&
                      gphes(1,inode,igaus)&
                      +gphes(2,inode,igaus)
              end do
           end do
        else
           do igaus=1,pgaus
              do inode=1,pnode
                 gplap(inode,igaus)=&
                       gphes(1,inode,igaus)&
                      +gphes(2,inode,igaus)&
                      +gphes(3,inode,igaus)
              end do
           end do
        end if
     else
        do igaus=1,pgaus
           do inode=1,pnode
              gplap(inode,igaus)=0.0_rp
           end do
        end do
     end if
  end if

  !-------------------------------------------------------------------
  !
  ! Shock capturing: CF and CD
  !
  !-------------------------------------------------------------------

  do igaus=1,pgaus
     CD(igaus) = 0.0_rp
     CF(igaus) = 0.0_rp
  end do
  if(kfl_shock/=0) then
     call shocel(&
          pnode,pgaus,ptopo,pelty,kfl_shock,shock,dtinv,gpden,&
          gpadv,gprhs, tau2,gpdif,gprea,elunk,gpsha,gpcar,&
          gplap,chale,CD,CF)
     if(kfl_taust==0) then
        do igaus=1,pgaus
           tau2(igaus)=0.0_rp
        end do
     end if
  end if

  !-------------------------------------------------------------------
  !
  ! Assembly in 2D for P1 elements and 1 Gauss point
  !
  !-------------------------------------------------------------------

  if(pgaus==1.and.pelty==TRI03) then

     k      = gpdif(1)                        ! k
     s      = gprea(1)                        ! s
     rho    = gpden(1)                        ! rho

     a1x    = gpadv(1,1)                      ! Residual convection
     a1y    = gpadv(2,1)
     a2x    = gpvel(1,1)+gpgrd(1,1)/rho       ! Adjoint convection
     a2y    = gpvel(2,1)+gpgrd(2,1)/rho 
     dkdx   = gpgrd(1,1)                      ! grad(k)
     dkdy   = gpgrd(2,1)

     ktau2  = tau2(1)*k
     rtau2  = tau2(1)*rho
     t1st2  = tau1(1) - tau2(1)*s
     rhodts = rho*dtinv + s

     c1     = t1st2*rhodts
     c2     = ktau2*rhodts
     c3     = rtau2*rhodts
     c4     = rho*t1st2
     c5     = rho*ktau2
     c6     = rtau2*rho
     c7     = k*(-t1st2 + 1.0_rp)
     c8     = -ktau2*k
     c9     = -c5
     c10    = k
     c11    = 1.0_rp

     e1     = a1x*c6 
     e2     = a1y*c6 
     axCF   = CF(1)*a1x
     ayCF   = CF(1)*a1y
     c10CD  = c10 + CD(1)

     A1(1)  = (c1)
     A2(1)  = (a2x*c3)
     A3(1)  = (a2y*c3)
     A5(1)  = (a1x*c4 + dkdx*c11)
     A6(1)  = (a1y*c4 + dkdy*c11)
     A8(1)  = a2x*e1 + c10CD  + axCF*a1x
     A9(1)  = a2x*e2          + axCF*a1y
     A11(1) = a2y*e1          + ayCF*a1x
     A12(1) = a2y*e2 + c10CD  + ayCF*a1y
     D1(1)  = t1st2
     D2(1)  = rtau2*a2x
     D3(1)  = rtau2*a2y
     D5(1)  = ktau2

     do inode=1,pnode
        x1 = gpvol(1)*(D1(1)*gpsha(inode,1)+D2(1)*gpcar(1,inode,1) &
             +D3(1)*gpcar(2,inode,1))
        do idofn=1,ndofn
           elrhs(idofn,inode)=elrhs(idofn,inode)+gprhs(idofn,1)*x1
        end do
     end do

     do inode=1,pnode
        do jnode=1,pnode
           elmat(inode,jnode)=elmat(inode,jnode)+gpvol(1)*(&
                & +gpsha(inode,1)*(    A1(1)*gpsha(jnode,1)      &  ! A1  Na Nb
                &                   +  A5(1)*gpcar(1,jnode,1)    &  ! A5  Na dNb/dx
                &                   +  A6(1)*gpcar(2,jnode,1)  ) &  ! A6  Na dNb/dy
                & +gpsha(jnode,1)*(    A2(1)*gpcar(1,inode,1)    &  ! A2  Nb dNa/dx
                &                   +  A3(1)*gpcar(2,inode,1)  ) &  ! A3  Nb dNa/dy
                & +gpcar(1,inode,1)*(  A8(1)*gpcar(1,jnode,1)    &  ! A8  dNa/dx dNb/dx
                &                   +  A9(1)*gpcar(2,jnode,1)  ) &  ! A9  dNa/dx dNb/dy
                & +gpcar(2,inode,1)*( A11(1)*gpcar(1,jnode,1)    &  ! A11 dNa/dy dNb/dx
                &                   + A12(1)*gpcar(2,jnode,1)  ) )  ! A12 dNa/dy dNb/dy
        end do
     end do

  !-------------------------------------------------------------------
  !
  ! Assembly in 1D
  !
  !-------------------------------------------------------------------

  else if(ndime==1) then

     do igaus=1,pgaus

        k          = gpdif(igaus)                        ! k
        s          = gprea(igaus)                        ! s
        rho        = gpden(igaus)                        ! rho

        a1x        = gpadv(1,igaus)                      ! Residual convection
        a2x        = gpvel(1,igaus)+gpgrd(1,igaus)/rho   ! Adjoint convection
        dkdx       = gpgrd(1,igaus)                      ! grad(k)

        ktau2      = tau2(igaus)*k
        rtau2      = tau2(igaus)*rho
        t1st2      = tau1(igaus) - tau2(igaus)*s
        rhodts     = rho*dtinv + s

        c1         = t1st2*rhodts
        c2         = ktau2*rhodts
        c3         = rtau2*rhodts
        c4         = rho*t1st2
        c5         = rho*ktau2
        c6         = rtau2*rho
        c7         = k*(-t1st2 + 1.0_rp)
        c8         = -ktau2*k
        c9         = -c5
        c10        = k
        c11        = 1.0_rp

        e1         = a1x*c6 
        axCF       = CF(igaus)*a1x
        c10CD      = c10 + CD(igaus)

        A1(igaus)  = (c1)
        A2(igaus)  = (a2x*c3)
        A5(igaus)  = (a1x*c4 + dkdx*c11)
        A8(igaus)  = a2x*e1 + c10CD  + axCF*a1x
        D1(igaus)  = t1st2
        D2(igaus)  = rtau2*a2x
        D5(igaus)  = ktau2

        B1(igaus)  = c2
        B2(igaus)  = c7
        B3(igaus)  = a1x*c5

        B6(igaus)  = c8
        B7(igaus)  = a2x*c9

     end do

     do igaus=1,pgaus
        do inode=1,pnode
           do idofn=1,ndofn
              elrhs(idofn,inode)=elrhs(idofn,inode)&
                   +gprhs(idofn,igaus)*gpvol(igaus)&
                   &            *(   D1(igaus)*gpsha(inode,igaus)    &
                   &                +D2(igaus)*gpcar(1,inode,igaus)  &
                   &                +D5(igaus)*gplap(inode,igaus)    )
           end do
        end do
     end do

     do igaus=1,pgaus
        do inode=1,pnode
           do jnode=1,pnode
              elmat(inode,jnode)=elmat(inode,jnode)+gpvol(igaus)*(&
                                ! first order term
                   & +gpsha(inode,igaus)*(    A1(igaus)*gpsha(jnode,igaus)      &  ! A1  Na Nb
                   &                       +  A5(igaus)*gpcar(1,jnode,igaus) )  &  ! A5  Na dNb/dx
                   & +gpsha(jnode,igaus)*(    A2(igaus)*gpcar(1,inode,igaus) )  &  ! A2  Nb dNa/dx
                   & +gpcar(1,inode,igaus)*(  A8(igaus)*gpcar(1,jnode,igaus) )  &  ! A8  dNa/dx dNb/dx
                                ! second order terms
                   & +gplap(inode,igaus)*(   B1(igaus)*gpsha(jnode,igaus)       &  ! B1  DNa Nb
                   &                       + B3(igaus)*gpcar(1,jnode,igaus)     &  ! B3  DNa dNb/dx
                   &                       + B6(igaus)*gplap(jnode,igaus)     ) &  ! B6  DNa DNb
                   & +gplap(jnode,igaus)*(   B2(igaus)*gpsha(inode,igaus)       &  ! B2  DNb Na
                   &                       + B7(igaus)*gpcar(1,inode,igaus)   ) )  ! B7  DNb dNa/dx
           end do
        end do
     end do

  !-------------------------------------------------------------------
  !
  ! Assembly in 2D
  !
  !-------------------------------------------------------------------

  else if(ndime==2) then

     do igaus=1,pgaus

        k          = gpdif(igaus)                        ! k
        s          = gprea(igaus)                        ! s
        rho        = gpden(igaus)                        ! rho

        a1x        = gpadv(1,igaus)                      ! Residual convection
        a1y        = gpadv(2,igaus)
        a2x        = gpvel(1,igaus)+gpgrd(1,igaus)/rho   ! Adjoint convection
        a2y        = gpvel(2,igaus)+gpgrd(2,igaus)/rho 
        dkdx       = gpgrd(1,igaus)                      ! grad(k)
        dkdy       = gpgrd(2,igaus)

        ktau2      = tau2(igaus)*k
        rtau2      = tau2(igaus)*rho
        t1st2      = tau1(igaus) - tau2(igaus)*s
        rhodts     = rho*dtinv + s

        c1         = t1st2*rhodts
        c2         = ktau2*rhodts
        c3         = rtau2*rhodts
        c4         = rho*t1st2
        c5         = rho*ktau2
        c6         = rtau2*rho
        c7         = k*(-t1st2 + 1.0_rp)
        c8         = -ktau2*k
        c9         = -c5
        c10        = k
        c11        = 1.0_rp

        e1         = a1x*c6 
        e2         = a1y*c6 
        axCF       = CF(igaus)*a1x
        ayCF       = CF(igaus)*a1y
        c10CD      = c10 + CD(igaus)

        A1(igaus)  = (c1)
        A2(igaus)  = (a2x*c3)
        A3(igaus)  = (a2y*c3)
        A5(igaus)  = (a1x*c4 + dkdx*c11)
        A6(igaus)  = (a1y*c4 + dkdy*c11)
        A8(igaus)  = a2x*e1 + c10CD  + axCF*a1x
        A9(igaus)  = a2x*e2          + axCF*a1y
        A11(igaus) = a2y*e1          + ayCF*a1x
        A12(igaus) = a2y*e2 + c10CD  + ayCF*a1y
        D1(igaus)  = t1st2
        D2(igaus)  = rtau2*a2x
        D3(igaus)  = rtau2*a2y
        D5(igaus)  = ktau2

        B1(igaus)  = c2
        B2(igaus)  = c7
        B3(igaus)  = a1x*c5
        B4(igaus)  = a1y*c5

        B6(igaus)  = c8
        B7(igaus)  = a2x*c9
        B8(igaus)  = a2y*c9

     end do

     do igaus=1,pgaus
        do inode=1,pnode
           do idofn=1,ndofn
              elrhs(idofn,inode)=elrhs(idofn,inode)&
                   +gprhs(idofn,igaus)*gpvol(igaus)&
                   &            *(   D1(igaus)*gpsha(inode,igaus)    &
                   &                +D2(igaus)*gpcar(1,inode,igaus)  &
                   &                +D3(igaus)*gpcar(2,inode,igaus)  &
                   &                +D5(igaus)*gplap(inode,igaus)    )
           end do
        end do
     end do

     do igaus=1,pgaus
        do inode=1,pnode
           do jnode=1,pnode
              elmat(inode,jnode)=elmat(inode,jnode)+gpvol(igaus)*(&
                                ! first order term
                   & +gpsha(inode,igaus)*(    A1(igaus)*gpsha(jnode,igaus)      &  ! A1  Na Nb
                   &                       +  A5(igaus)*gpcar(1,jnode,igaus)    &  ! A5  Na dNb/dx
                   &                       +  A6(igaus)*gpcar(2,jnode,igaus)  ) &  ! A6  Na dNb/dy
                   & +gpsha(jnode,igaus)*(    A2(igaus)*gpcar(1,inode,igaus)    &  ! A2  Nb dNa/dx
                   &                       +  A3(igaus)*gpcar(2,inode,igaus)  ) &  ! A3  Nb dNa/dy
                   & +gpcar(1,inode,igaus)*(  A8(igaus)*gpcar(1,jnode,igaus)    &  ! A8  dNa/dx dNb/dx
                   &                       +  A9(igaus)*gpcar(2,jnode,igaus)  ) &  ! A9  dNa/dx dNb/dy
                   & +gpcar(2,inode,igaus)*( A11(igaus)*gpcar(1,jnode,igaus)    &  ! A11 dNa/dy dNb/dx
                   &                       + A12(igaus)*gpcar(2,jnode,igaus)  ) &  ! A12 dNa/dy dNb/dy
                                ! second order terms
                   & +gplap(inode,igaus)*(   B1(igaus)*gpsha(jnode,igaus)       &  ! B1  DNa Nb
                   &                       + B3(igaus)*gpcar(1,jnode,igaus)     &  ! B3  DNa dNb/dx
                   &                       + B4(igaus)*gpcar(2,jnode,igaus)     &  ! B4  DNa dNb/dy
                   &                       + B6(igaus)*gplap(jnode,igaus)     ) &  ! B6  DNa DNb
                   & +gplap(jnode,igaus)*(   B2(igaus)*gpsha(inode,igaus)       &  ! B2  DNb Na
                   &                       + B7(igaus)*gpcar(1,inode,igaus)     &  ! B7  DNb dNa/dx
                   &                       + B8(igaus)*gpcar(2,inode,igaus)   ) )  ! B8  DNb dNa/dy
           end do
        end do
     end do

  else

     !-------------------------------------------------------------------
     !
     ! Assembly in 3D
     !
     !-------------------------------------------------------------------

     do igaus=1,pgaus

        k          = gpdif(igaus)
        s          = gprea(igaus)
        rho        = gpden(igaus)

        a1x        = gpadv(1,igaus)
        a1y        = gpadv(2,igaus)
        a1z        = gpadv(3,igaus)
        a2x        = gpvel(1,igaus)+gpgrd(1,igaus)/rho
        a2y        = gpvel(2,igaus)+gpgrd(2,igaus)/rho
        a2z        = gpvel(3,igaus)+gpgrd(3,igaus)/rho
        dkdx       = gpgrd(1,igaus) 
        dkdy       = gpgrd(2,igaus)
        dkdz       = gpgrd(3,igaus)

        ktau2      = tau2(igaus)*k
        rtau2      = tau2(igaus)*rho
        t1st2      = tau1(igaus) - tau2(igaus)*s
        rhodts     = rho*dtinv + s

        c1         = t1st2*rhodts
        c2         = ktau2*rhodts
        c3         = rtau2*rhodts
        c4         = rho*t1st2
        c5         = rho*ktau2
        c6         = rtau2*rho
        c7         = k*(-t1st2 + 1.0_rp)
        c8         = -ktau2*k
        c9         = -c5
        c10        = k
        c11        = 1.0_rp

        e1         = a1x*c6 
        e2         = a1y*c6 
        e3         = a1z*c6 
        axCF       = CF(igaus)*a1x
        ayCF       = CF(igaus)*a1y
        azCF       = CF(igaus)*a1z
        c10CD      = c10 + CD(igaus)

        A1(igaus)  = (c1)
        A2(igaus)  = (a2x*c3)
        A3(igaus)  = (a2y*c3)
        A4(igaus)  = (a2z*c3)
        A5(igaus)  = (a1x*c4 + dkdx*c11)
        A6(igaus)  = (a1y*c4 + dkdy*c11)
        A7(igaus)  = (a1z*c4 + dkdz*c11)
        A8(igaus)  = a2x*e1 + c10CD  + axCF*a1x
        A9(igaus)  = a2x*e2          + axCF*a1y
        A10(igaus) = a2x*e3          + axCF*a1z
        A11(igaus) = a2y*e1          + ayCF*a1x
        A12(igaus) = a2y*e2 + c10CD  + ayCF*a1y
        A13(igaus) = a2y*e3          + ayCF*a1z
        A14(igaus) = a2z*e1          + azCF*a1x
        A15(igaus) = a2z*e2          + azCF*a1y
        A16(igaus) = a2z*e3 + c10CD  + azCF*a1z
        D1(igaus)  = t1st2
        D2(igaus)  = rtau2*a2x
        D3(igaus)  = rtau2*a2y
        D4(igaus)  = rtau2*a2z
        D5(igaus)  = ktau2

        B1(igaus)  = c2
        B2(igaus)  = c7
        B3(igaus)  = a1x*c5
        B4(igaus)  = a1y*c5
        B5(igaus)  = a1z*c5

        B6(igaus)  = c8
        B7(igaus)  = a2x*c9
        B8(igaus)  = a2y*c9
        B9(igaus)  = a2z*c9

     end do

     do igaus=1,pgaus
        do inode=1,pnode
           do idofn=1,ndofn
              elrhs(idofn,inode)=elrhs(idofn,inode)&
                   +gprhs(idofn,igaus)*gpvol(igaus)&
                   &            *(   D1(igaus)*gpsha(inode,igaus)    &
                   &                +D2(igaus)*gpcar(1,inode,igaus)  &
                   &                +D3(igaus)*gpcar(2,inode,igaus)  &
                   &                +D4(igaus)*gpcar(3,inode,igaus)  &
                   &                +D5(igaus)*gplap(inode,igaus)    )
           end do
        end do
     end do

     do igaus=1,pgaus
        do inode=1,pnode
           do jnode=1,pnode
              elmat(inode,jnode)=elmat(inode,jnode)+gpvol(igaus)*(&
                                ! first order term
                   & +gpsha(inode,igaus)*(    A1(igaus)*gpsha(jnode,igaus)      & ! N_a N_b
                   &                       +  A5(igaus)*gpcar(1,jnode,igaus)    & ! N_a dN_b/dx
                   &                       +  A6(igaus)*gpcar(2,jnode,igaus)    & ! N_a dN_b/dy
                   &                       +  A7(igaus)*gpcar(3,jnode,igaus)  ) & ! N_a dN_b/dz
                   & +gpsha(jnode,igaus)*(    A2(igaus)*gpcar(1,inode,igaus)    & ! N_b dN_a/dx
                   &                       +  A3(igaus)*gpcar(2,inode,igaus)    &
                   &                       +  A4(igaus)*gpcar(3,inode,igaus)  ) &
                   & +gpcar(1,inode,igaus)*(  A8(igaus)*gpcar(1,jnode,igaus)    &
                   &                       +  A9(igaus)*gpcar(2,jnode,igaus)    &
                   &                       + A10(igaus)*gpcar(3,jnode,igaus)  ) &
                   & +gpcar(2,inode,igaus)*( A11(igaus)*gpcar(1,jnode,igaus)    &
                   &                       + A12(igaus)*gpcar(2,jnode,igaus)    &
                   &                       + A13(igaus)*gpcar(3,jnode,igaus)  ) &
                   & +gpcar(3,inode,igaus)*( A14(igaus)*gpcar(1,jnode,igaus)    &
                   &                       + A15(igaus)*gpcar(2,jnode,igaus)    &
                   &                       + A16(igaus)*gpcar(3,jnode,igaus)  ) &
                                ! second order terms
                   & +gplap(inode,igaus)*(   B1(igaus)*gpsha(jnode,igaus)       &
                   &                       + B3(igaus)*gpcar(1,jnode,igaus)     &
                   &                       + B4(igaus)*gpcar(2,jnode,igaus)     &
                   &                       + B5(igaus)*gpcar(3,jnode,igaus)     &
                   &                       + B6(igaus)*gplap(jnode,igaus)     ) &
                   & +gplap(jnode,igaus)*(   B2(igaus)*gpsha(inode,igaus)       &
                   &                       + B7(igaus)*gpcar(1,inode,igaus)     &
                   &                       + B8(igaus)*gpcar(2,inode,igaus)     &
                   &                       + B9(igaus)*gpcar(3,inode,igaus)   ) )
           end do
        end do
     end do

  end if

end subroutine elmcel

subroutine shocel(&
     pnode,pgaus,ptopo,pelty,kfl_shock,shock,dtinv,gpden,&
     gpadv,gprhs,gpstp,gpdif,gprea,elunk,gpsha,gpcar,&
     gplap,chale,CD,CF)
  !-----------------------------------------------------------------------
  !****f* mathru/shocel
  ! NAME
  !   shocel
  ! DESCRIPTION
  !   This routine computes the contribution from the shock capturing
  ! OUTPUT
  ! USES
  ! USED BY
  !   elmcel
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mgaus,mnode
  use def_elmtyp, only       :  TRI03,TET04
  implicit none  
  integer(ip), intent(in)    :: pnode,pgaus,ptopo,pelty,kfl_shock
  real(rp),    intent(in)    :: shock,dtinv
  real(rp),    intent(in)    :: gpadv(ndime,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gplap(pnode,pgaus)
  real(rp),    intent(in)    :: elunk(pnode)
  real(rp),    intent(in)    :: gpden(pgaus),gprhs(pgaus)
  real(rp),    intent(in)    :: gpdif(pgaus)
  real(rp),    intent(in)    :: gprea(pgaus)
  real(rp),    intent(in)    :: chale,gpstp(pgaus)
  real(rp),    intent(out)   :: CD(pgaus),CF(pgaus)                 
  real(rp)                   :: gpve2,grnor,vepar,facta
  real(rp)                   :: gpgru(3),factt,gpres
  real(rp)                   :: F1,F2,F3,F4,cdv2,umbra
  integer(ip)                :: inode,igaus
  !
  ! Factor according to element topology
  !
  if(ptopo==0.or.ptopo==-1) then
     factt=1.5_rp                ! Quadrilateral/Hexahedra or Bars
  else if(ptopo==1) then
     factt=0.7_rp                ! Triangles/Tetrahedra
  else
     factt=1.0_rp                ! Other elements
  end if
  factt=0.5_rp*factt
  umbra=1.0e-8_rp
  !
  ! Isotropic/anisotropic shock capturing
  !
  facta=real(kfl_shock-1_ip,rp)
  !
  ! Coarse grid residual |R(u)|
  ! |R(u)| = | f - rho*u/(theta*dt) - rho*a'.grad(u) + k*Lap(u) - s*u|
  ! a' = a - grad(k)/rho 
  ! f  = Q + rho*u^n/(theta*dt)
  !
  if(pgaus==1.and.pelty==TRI03) then

     gpres =  gprhs(1)                                    ! Residual
     F1    = -gpden(1)*dtinv - gprea(1)
     F2    = -gpden(1)*gpadv(1,1) 
     F3    = -gpden(1)*gpadv(2,1)     
     do inode=1,pnode
        gpres = gpres+elunk(inode)*(&
             &  +F1*gpsha(inode,1)&
             &  +F2*gpcar(1,inode,1)&
             &  +F3*gpcar(2,inode,1)&
             &  +gpdif(1)*gplap(inode,1))
     end do
     gpres    = abs(gpres)
     gpve2    =   gpadv(1,1)*gpadv(1,1)&              ! a^2
          &     + gpadv(2,1)*gpadv(2,1)
     gpgru(1) = 0.0_rp
     gpgru(2) = 0.0_rp
     do inode=1,pnode                                         ! grad(u)
        gpgru(1) = gpgru(1)+gpcar(1,inode,1)*elunk(inode)
        gpgru(2) = gpgru(2)+gpcar(2,inode,1)*elunk(inode)
     end do
     grnor=sqrt( gpgru(1)*gpgru(1) + gpgru(2)*gpgru(2) )      ! | grad(u) | 

     if( grnor>umbra ) then
        vepar = gpres/grnor
        CD(1) = factt*shock*chale*vepar-gpdif(1)
        if( CD(1)>0.0_rp .and. gpve2>1.0d-6 ) then
           cdv2  = CD(1)/gpve2
           CF(1) = max(0.0_rp,cdv2-gpden(1)*gpstp(1)*gpden(1)) - cdv2
           CF(1) = facta*CF(1)
        end if
     end if

  else if(ndime==1) then

     do igaus=1,pgaus
        gpres =  gprhs(igaus)                                    ! Residual
        F1    = -gpden(igaus)*dtinv - gprea(igaus)
        F2    = -gpden(igaus)*gpadv(1,igaus)    
        do inode=1,pnode
           gpres = gpres+elunk(inode)*(&
                &  +F1*gpsha(inode,igaus)&
                &  +F2*gpcar(1,inode,igaus)&
                &  +gpdif(igaus)*gplap(inode,igaus))
        end do
        gpres    = abs(gpres)
        gpve2    =   gpadv(1,igaus)*gpadv(1,igaus)               ! a^2
        gpgru(1) = 0.0_rp
        do inode=1,pnode                                         ! grad(u)
           gpgru(1) = gpgru(1)+gpcar(1,inode,igaus)*elunk(inode)
        end do
        grnor=abs(gpgru(1))                                      ! | grad(u) | 

        if( grnor>umbra ) then
           vepar     = gpres/grnor
           CD(igaus) = factt*shock*chale*vepar-gpdif(igaus)
           if( CD(igaus)>0.0_rp .and. gpve2>1.0d-6 ) then
              cdv2      = CD(igaus)/gpve2
              CF(igaus) = max(0.0_rp,cdv2-gpden(igaus)*gpstp(igaus)*gpden(igaus)) - cdv2
              CF(igaus) = facta*CF(igaus)
           end if
        end if

     end do

  else if(ndime==2) then

     do igaus=1,pgaus
        gpres =  gprhs(igaus)                                    ! Residual
        F1    = -gpden(igaus)*dtinv - gprea(igaus)
        F2    = -gpden(igaus)*gpadv(1,igaus) 
        F3    = -gpden(igaus)*gpadv(2,igaus)     
        do inode=1,pnode
           gpres = gpres+elunk(inode)*(&
                &  +F1*gpsha(inode,igaus)&
                &  +F2*gpcar(1,inode,igaus)&
                &  +F3*gpcar(2,inode,igaus)&
                &  +gpdif(igaus)*gplap(inode,igaus))
        end do
        gpres    = abs(gpres)
        gpve2    =   gpadv(1,igaus)*gpadv(1,igaus)&              ! a^2
             &     + gpadv(2,igaus)*gpadv(2,igaus)
        gpgru(1) = 0.0_rp
        gpgru(2) = 0.0_rp
        do inode=1,pnode                                         ! grad(u)
           gpgru(1) = gpgru(1)+gpcar(1,inode,igaus)*elunk(inode)
           gpgru(2) = gpgru(2)+gpcar(2,inode,igaus)*elunk(inode)
        end do
        grnor=sqrt( gpgru(1)*gpgru(1) + gpgru(2)*gpgru(2) )      ! | grad(u) | 

        if( grnor>umbra ) then
           vepar     = gpres/grnor
           CD(igaus) = factt*shock*chale*vepar-gpdif(igaus)
           if( CD(igaus)>0.0_rp .and. gpve2>1.0d-6 ) then
              cdv2      = CD(igaus)/gpve2
              CF(igaus) = max(0.0_rp,cdv2-gpden(igaus)*gpstp(igaus)*gpden(igaus)) - cdv2
              CF(igaus) = facta*CF(igaus)
           end if
        end if

     end do

  else

     do igaus=1,pgaus
        gpres =  gprhs(igaus)                                    ! Residual
        F1    = -gpden(igaus)*dtinv - gprea(igaus)
        F2    = -gpden(igaus)*gpadv(1,igaus) 
        F3    = -gpden(igaus)*gpadv(2,igaus)    
        F4    = -gpden(igaus)*gpadv(3,igaus)    
        do inode=1,pnode
           gpres = gpres+elunk(inode)*(&
                &  +F1*gpsha(inode,igaus)&
                &  +F2*gpcar(1,inode,igaus)&
                &  +F3*gpcar(2,inode,igaus)&
                &  +F4*gpcar(3,inode,igaus)&
                &  +gpdif(igaus)*gplap(inode,igaus))
        end do
        gpres    = abs(gpres)
        gpve2    =   gpadv(1,igaus)*gpadv(1,igaus)&              ! a^2
             &     + gpadv(2,igaus)*gpadv(2,igaus)& 
             &     + gpadv(3,igaus)*gpadv(3,igaus)
        gpgru(1) = 0.0_rp
        gpgru(2) = 0.0_rp
        gpgru(3) = 0.0_rp
        do inode=1,pnode                                         ! grad(u)
           gpgru(1) = gpgru(1)+gpcar(1,inode,igaus)*elunk(inode)
           gpgru(2) = gpgru(2)+gpcar(2,inode,igaus)*elunk(inode)
           gpgru(3) = gpgru(3)+gpcar(3,inode,igaus)*elunk(inode)
        end do
        grnor=sqrt(   gpgru(1)*gpgru(1) &                        ! | grad(u) | 
             &      + gpgru(2)*gpgru(2) &
             &      + gpgru(3)*gpgru(3) ) 

        if( grnor>umbra ) then
           vepar     = gpres/grnor
           CD(igaus) = factt*0.5_rp*shock*chale*vepar-gpdif(igaus)
           if( CD(igaus)>0.0_rp .and. gpve2>1.0d-6 ) then
              cdv2      = CD(igaus)/gpve2
              CF(igaus) = max(0.0_rp,cdv2-gpden(igaus)*gpstp(igaus)*gpden(igaus)) - cdv2
              CF(igaus) = facta*CF(igaus)
           end if
        end if

     end do

  end if

  do igaus=1,pgaus
     CD(igaus)=max(CD(igaus),0.0_rp)
     CF(igaus)=max(CF(igaus),0.0_rp)
  end do

end subroutine shocel
