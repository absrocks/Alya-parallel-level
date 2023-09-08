subroutine elmadr(&
     itask,ielem,pnode,pgaus,ptopo,plapl,pelty,ndofn,nddif,lnods,kfl_shock,&
     kfl_taust,kfl_stabi,kfl_sgsco,kfl_sgsti,kfl_limit,&
     staco,gpdif,gprea,gpden,gpvel,gppro,gpvol,gpgrd,gprhs,&
     gpsha,gpcar,gphes,elunk,elcod,chale, tau, dtinv,&
     shock,bemol,gpdiv,gplap,gpsgs,elmat,elrhs)
  !-----------------------------------------------------------------------
  !****f* mathru/elmadx
  ! NAME
  !    elmadx
  ! DESCRIPTION
  !
  !    Elemental operations for the ADR equations:
  !
  !    ITASK = 1 ............. Assemble matrix and RHS
  !            4 ............. Compute SGS and projection
  !            5 ............. Postprocess projection
  !
  !    The stabilizations are:
  !
  !    KFL_STABI = -2 ....... Galerkin
  !                -1 ....... SU
  !                 0 ....... ASGS
  !                 1 ....... Full OSS
  !                 2 ....... OSS (only convection)
  !                 3 ....... OSS (only convection) with limiter
  !
  ! OUTPUT 
  !    ELMAT ... Matrix
  !    ELRHS ... RHS
  ! USES
  ! USED BY
  !    *_elmope 
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_elmtyp, only       :  ELEXT
  use def_master, only       :  rhsid
  use def_domain, only       :  ndime,ntens,mnode,mgaus,lelch
  use def_elmtyp, only       :  TRI03
  use mod_tauadr, only       :  tauadr
  implicit none

  integer(ip), intent(in)    :: itask,ielem,pnode,pgaus,ptopo,plapl,pelty,ndofn
  integer(ip), intent(in)    :: nddif
  integer(ip), intent(in)    :: lnods(pnode)
  integer(ip), intent(in)    :: kfl_shock,kfl_taust,kfl_stabi
  integer(ip), intent(in)    :: kfl_sgsco,kfl_sgsti,kfl_limit
  real(rp),    intent(in)    :: staco(3)
  real(rp),    intent(in)    :: gpdif(nddif,pgaus)
  real(rp),    intent(inout) :: gprea(pgaus)
  real(rp),    intent(in)    :: gpden(pgaus)
  real(rp),    intent(in)    :: gpvel(ndime,pgaus)
  real(rp),    intent(inout) :: gppro(pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: gpgrd(ndime,pgaus)
  real(rp),    intent(inout) :: gprhs(ndofn,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(in)    :: elunk(pnode,*)
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(in)    :: chale(2)
  real(rp),    intent(inout) :: tau(pgaus)
  real(rp),    intent(in)    :: dtinv
  real(rp),    intent(in)    :: shock
  real(rp),    intent(in)    :: bemol
  real(rp),    intent(in)    :: gpdiv(pgaus)     ! is only used if bemol /= 0    
  real(rp),    intent(out)   :: gplap(pnode,pgaus)
  real(rp),    intent(inout) :: gpsgs(pgaus,*)   ! component 2 is intent in - component 1 is out
  real(rp),    intent(out)   :: elmat(pnode,pnode)
  real(rp),    intent(out)   :: elrhs(ndofn,pnode)

  integer(ip)                :: inode,jnode,igaus,idofn,idime,ireac
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
  real(rp)                   :: c1 ,c2, c3, c4, c5, c6, c7, c8, c9, c10, c11,c12
  real(rp)                   :: dkdx,dkdy,dkdz,rho,div,bemo1,beta,alpha
  real(rp)                   :: s,k,e1,e2,e3,ktau,rtau,t1st2
  real(rp)                   :: axCF,ayCF,azCF,rhodts,c10CD,gpres
  real(rp)                   :: a1x,a1y,a1z,a2x,a2y,a2z,gpnve,x1,g1,g2,g3,g4
  real(rp)                   :: ax,ay,az,b,rb,rb1,rrtau,dummr
  real(rp)                   :: atx,aty,atz,gpre2(mgaus)
  real(rp)                   :: kx,ky,kz,dkxdx,dkydy,dkzdz

  !-------------------------------------------------------------------
  !
  ! If SUPG stabilization method
  !
  !-------------------------------------------------------------------
 
  


  !-------------------------------------------------------------------
  !
  ! Add time derivative to RHS: f = f + rho * u^n / dt
  !
  !-------------------------------------------------------------------

  if( dtinv /= 0.0_rp ) then
     do igaus = 1,pgaus
        c1 = 0.0_rp
        do inode = 1,pnode
           c1 = c1 + elunk(inode,2) * gpsha(inode,igaus)
        end do
        gprhs(1,igaus) = gprhs(1,igaus) + gpden(igaus) * dtinv * c1
     end do
  end if

  !-------------------------------------------------------------------
  !
  ! Redefine advection for ASGS or Full OSS
  !
  !-------------------------------------------------------------------

  if( kfl_stabi ==  0 .or. kfl_stabi == 1 ) then
     do igaus=1,pgaus
        do idime=1,ndime
           gpadv(idime,igaus) = gpvel(idime,igaus) !! DMM - gpgrd(idime,igaus)/gpden(igaus)
        end do
     end do
  else
     do igaus=1,pgaus
        do idime=1,ndime
           gpadv(idime,igaus) = gpvel(idime,igaus)
        end do
     end do
  end if

  !-------------------------------------------------------------------
  !
  ! Compute stabilization
  !
  !-------------------------------------------------------------------

  if( kfl_stabi /= -2 .or. kfl_shock /= 0_ip ) then

     if( nddif == ndime ) then
        do igaus = 1,pgaus
           k = 0.0_rp
           do idime = 1,ndime
              k = k + gpdif(idime,igaus) 
           end do
           k = k / real(ndime,rp)
           call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip) 
           gpnve = gpden(igaus)*gpnve
           call tauadr(&
                kfl_taust,staco,gpnve,k,gprea(igaus),&
                chale(1),chale(2),tau(igaus))          
           if( kfl_sgsti /= 0 .and. tau(igaus) /= 0.0_rp ) &
                tau(igaus) = 1.0_rp / ( gpden(igaus) * dtinv + 1.0_rp / tau(igaus) )
        end do
     else 
        do igaus = 1,pgaus
           call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip) 
           gpnve = gpden(igaus)*gpnve
           call tauadr(&
                kfl_taust,staco,gpnve,gpdif(1,igaus),gprea(igaus),&
                chale(1),chale(2),tau(igaus))
           if( kfl_sgsti /= 0 .and. tau(igaus) /= 0.0_rp ) &
                tau(igaus) = 1.0_rp / ( gpden(igaus) * dtinv + 1.0_rp / tau(igaus) )
        end do
     end if

  else

     do igaus = 1,pgaus
        tau(igaus) = 0.0_rp
     end do

  end if

  !-------------------------------------------------------------------
  !
  ! Laplacian matrix: GPLAP
  !
  !-------------------------------------------------------------------

  if( plapl /= 2 .and. ( kfl_stabi == 0 .or. kfl_stabi == 1 .or. kfl_shock /= 0 ) ) then

     if( plapl == 1 ) then
        if( ndime == 1 ) then
           do igaus=1,pgaus
              do inode=1,pnode
                 gplap(inode,igaus)=&
                      gphes(1,inode,igaus)
              end do
           end do
        else if( ndime == 2 ) then
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

  if( itask /= 5 ) then

     if( kfl_stabi /= 0 ) then

        !-------------------------------------------------------------------
        !
        ! OSS: Limiter
        !
        ! Soto:
        ! || rho*tau1'*a.grad(u) - P || 
        ! / ( || rho*tau1'*a.grad(u) || + || P || )
        !
        !-------------------------------------------------------------------

        if( kfl_limit == -1 .or. kfl_stabi == -1 .or. kfl_stabi == -2 ) then

           do igaus = 1,pgaus
              gppro(igaus) = 0.0_rp
           end do

        else if( kfl_limit > 0 ) then

           do igaus = 1,pgaus
              c1 = 0.0_rp
              c2 = gppro(igaus)
              do inode = 1,pnode
                 c4 = 0.0_rp
                 do idime = 1,ndime
                    c4 = c4 + gpvel(idime,igaus) * gpcar(idime,inode,igaus)
                 end do
                 c1 = c1 + c4 * elunk(inode,1)
              end do
              c1 = gpden(igaus) * tau(igaus) * c1 
              c3 = abs(c1) + abs(c2) 
              if( c3 > 1.0e-6_rp ) then
                 beta  = abs( c1 - c2 ) / c3
              else
                 beta  = 0.0_rp
              end if
              if( kfl_limit == 1 ) then
                 alpha = min(1.0_rp,2.0_rp*(1.0_rp-beta))            ! SOTO
                 !alpha = 1.0_rp-beta
              else if( kfl_limit == 2 ) then
                 alpha = 0.5_rp*(tanh(20.0_rp*(beta-0.8_rp))+1.0_rp) ! Very diffusive
              end if
              gppro(igaus) = alpha * gppro(igaus)
           end do

        end if

        if( kfl_sgsti /= 0 ) then
           !
           ! alpha * Pi + tau1' *rho * u' / dt
           !
           do igaus = 1,pgaus
              gppro(igaus) = gppro(igaus) &
                   + tau(igaus) * gpden(igaus) * dtinv * gpsgs(igaus,2)
           end do
        end if

     end if
  end if

  if( itask == 1 ) then

     !-------------------------------------------------------------------
     !
     ! Shock capturing: CF and CD
     !
     !-------------------------------------------------------------------

     do igaus = 1,pgaus
        CD(igaus) = 0.0_rp
        CF(igaus) = 0.0_rp
     end do
     if( kfl_shock /= 0 ) then
        call elmshc(&
             pnode,pgaus,ptopo,pelty,nddif,kfl_shock,shock,dtinv,gpden,&
             gpvel,gprhs, tau,gpdif,gprea,gpgrd,elunk,gpsha,gpcar,&
             gplap,gphes,chale,CD,CF)
        if( kfl_stabi == -2 ) then
           do igaus = 1,pgaus
              tau(igaus) = 0.0_rp
           end do
        end if
     end if

     !-------------------------------------------------------------------
     !
     ! Initialization: ELMAT and ELRHS
     !
     !-------------------------------------------------------------------

     do jnode = 1,pnode
        do idofn = 1,ndofn
           elrhs(idofn,jnode) = 0.0_rp
        end do
        do inode = 1,pnode
           elmat(inode,jnode) = 0.0_rp
        end do
     end do

     bemo1 = 1.0_rp - bemol

     !-------------------------------------------------------------------
     !
     ! Use Lumped mass to compute reaction
     !
     !-------------------------------------------------------------------

     ireac = 0
     if( ireac == 1 ) then
        do igaus = 1,pgaus
           gpre2(igaus) = gprea(igaus)
           gprea(igaus) = 0.0_rp
        end do
     end if

     !-------------------------------------------------------------------
     !
     ! ASGS + Full OSS: Assembly in 2D for P1 elements and 1 Gauss point
     !
     !-------------------------------------------------------------------

     if( kfl_stabi == 0 .or. kfl_stabi == 1 ) then

        !if( pgaus == 1 .and. pelty == TRI03 ) then
        !
        !   if( nddif == 1 ) then
        !   else
        !   end if

        if( ndime == 1 ) then

           !-------------------------------------------------------------------
           !
           ! ASGS: Assembly in 1D
           !
           !-------------------------------------------------------------------

           do igaus = 1,pgaus

              k          = gpdif(1,igaus)                      ! k
              s          = gprea(igaus)                        ! s
              rho        = gpden(igaus)                        ! rho

              ax         = gpvel(1,igaus)
              a1x        = gpadv(1,igaus)                      ! Residual convection
              a2x        = gpvel(1,igaus)+gpgrd(1,igaus)/rho   ! Adjoint convection
              dkdx       = gpgrd(1,igaus)                      ! grad(k)

              ktau       = tau(igaus)*k                        ! tau*k
              rtau       = tau(igaus)*rho                      ! tau*rho
              t1st2      = 1.0_rp - tau(igaus)*s               ! 1 - tau*s
              rhodts     = rho*dtinv + s                       ! rho/dt + s

              c1         = t1st2*rhodts - bemol*gpdiv(igaus)   ! (1 - s*tau)*(rho/dt + s) - bemol*div(rho*a)
              c2         = ktau*rhodts                         ! k*tau*(rho/dt + s)
              c3         = rtau*rhodts                         ! rho*tau*(rho/dt + s)
              c4         = rho*t1st2                           ! rho*(1 - tau*s)
              c5         = rho*ktau                            ! rho*k*tau
              c6         = rtau*rho                            ! rho*rho*tau
              c7         = k*s*tau(igaus)                      ! k*s*tau
              c8         = -ktau*k                             ! -k^2*tau
              c9         = -c5                                 ! -rho*k*tau
              c10        = k                                   ! k
              c11        = 1.0_rp                              ! 1
              c12        = -bemol                              ! -bemol

              g2         = gprhs(1,igaus)
              g3         = gppro(igaus)
              g4         = g3 + tau(igaus)*g2              

              e1         = a1x*c6 
              axCF       = CF(igaus)*a1x
              c10CD      = c10 + CD(igaus)

              A1(igaus)  = c1
              A2(igaus)  = -bemol*rho*ax + a2x*c3
              A5(igaus)  =  bemo1*rho*ax - a1x*s*rtau
              A8(igaus)  = a2x*e1 + c10CD  + axCF*a1x

              D1(igaus)  = g2-s*g4
              D2(igaus)  = a2x*rho*g4
              D5(igaus)  = k*g4

              B1(igaus)  = c2
              B2(igaus)  = c7
              B3(igaus)  = a1x*c5

              B6(igaus)  = c8
              B7(igaus)  = a2x*c9

           end do

           do igaus=1,pgaus
              do inode=1,pnode
                 elrhs(1,inode)=elrhs(1,inode) + gpvol(igaus)&
                      &            *(   D1(igaus)*gpsha(inode,igaus)    &
                      &                +D2(igaus)*gpcar(1,inode,igaus)  &
                      &                +D5(igaus)*gplap(inode,igaus)    )
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

        else if( ndime == 2 ) then

           !-------------------------------------------------------------------
           !
           ! ASGS: Assembly in 2D
           !
           !-------------------------------------------------------------------

           do igaus=1,pgaus

              k          = gpdif(1,igaus)                        ! k
              s          = gprea(igaus)                        ! s
              rho        = gpden(igaus)                        ! rho

              ax         = gpvel(1,igaus)                      ! Convection
              ay         = gpvel(2,igaus)
              a1x        = gpadv(1,igaus)                      ! Residual convection
              a1y        = gpadv(2,igaus)
              a2x        = gpvel(1,igaus)+gpgrd(1,igaus)/rho   ! Adjoint convection
              a2y        = gpvel(2,igaus)+gpgrd(2,igaus)/rho 
              dkdx       = gpgrd(1,igaus)                      ! grad(k)
              dkdy       = gpgrd(2,igaus)

              ktau       = tau(igaus)*k                        ! tau*k
              rtau       = tau(igaus)*rho                      ! tau*rho
              t1st2      = 1.0_rp - tau(igaus)*s               ! 1 - tau*s
              rhodts     = rho*dtinv + s                       ! rho/dt + s

              c1         = t1st2*rhodts - bemol*gpdiv(igaus)   ! (1 - s*tau)*(rho/dt + s) - bemol*div(rho*a)
              c2         = ktau*rhodts                         ! k*tau*(rho/dt + s)
              c3         = rtau*rhodts                         ! rho*tau*(rho/dt + s)
              c4         = rho*t1st2                           ! rho*(1 - tau*s)
              c5         = rho*ktau                            ! rho*k*tau
              c6         = rtau*rho                            ! rho*rho*tau
              c7         = k*s*tau(igaus)                      ! k*s*tau
              c8         = -ktau*k                             ! -k^2*tau
              c9         = -c5                                 ! -rho*k*tau
              c10        = k                                   ! k
              c11        = 1.0_rp                              ! 1
              c12        = -bemol                              ! -bemol

              g2         = gprhs(1,igaus)
              g3         = gppro(igaus)
              g4         = g3 + tau(igaus)*g2              

              e1         = a1x*c6 
              e2         = a1y*c6 
              axCF       = CF(igaus)*a1x
              ayCF       = CF(igaus)*a1y
              c10CD      = c10 + CD(igaus)

              A1(igaus)  = (c1)
              A2(igaus)  = -bemol*rho*ax + a2x*c3
              A3(igaus)  = -bemol*rho*ay + a2y*c3
              A5(igaus)  =  bemo1*rho*ax - a1x*s*rtau
              A6(igaus)  =  bemo1*rho*ay - a1y*s*rtau              
              A8(igaus)  =  a2x*e1 + c10CD  + axCF*a1x
              A9(igaus)  =  a2x*e2          + axCF*a1y
              A11(igaus) =  a2y*e1          + ayCF*a1x
              A12(igaus) =  a2y*e2 + c10CD  + ayCF*a1y

              B1(igaus)  = c2
              B2(igaus)  = c7
              B3(igaus)  = a1x*c5
              B4(igaus)  = a1y*c5

              B6(igaus)  = c8
              B7(igaus)  = a2x*c9
              B8(igaus)  = a2y*c9

              D1(igaus)  = g2-s*g4
              D2(igaus)  = a2x*rho*g4
              D3(igaus)  = a2y*rho*g4
              D5(igaus)  = k*g4

           end do

           do igaus=1,pgaus
              do inode=1,pnode
                 do idofn=1,ndofn
                    elrhs(idofn,inode)=elrhs(idofn,inode) + gpvol(igaus)&
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
           ! ASGS: Assembly in 3D
           !
           !-------------------------------------------------------------------

           do igaus = 1,pgaus

              k          = gpdif(1,igaus)
              kx         = gpdif(1,igaus)
              ky         = gpdif(1,igaus)
              kz         = gpdif(1,igaus)
              !kx         = gpdif(1,igaus)
              !ky         = gpdif(2,igaus)
              !kz         = gpdif(3,igaus)
              s          = gprea(igaus)
              rho        = gpden(igaus)

              ax         = gpvel(1,igaus)
              ay         = gpvel(2,igaus)
              az         = gpvel(3,igaus)
              a1x        = gpadv(1,igaus)
              a1y        = gpadv(2,igaus)
              a1z        = gpadv(3,igaus)
              a2x        = gpvel(1,igaus)+gpgrd(1,igaus)/rho
              a2y        = gpvel(2,igaus)+gpgrd(2,igaus)/rho
              a2z        = gpvel(3,igaus)+gpgrd(3,igaus)/rho
              dkdx       = gpgrd(1,igaus) 
              dkdy       = gpgrd(2,igaus)
              dkdz       = gpgrd(3,igaus)

              ktau       = tau(igaus)*k                        ! tau*k
              rtau       = tau(igaus)*rho                      ! tau*rho
              t1st2      = 1.0_rp - tau(igaus)*s               ! 1 - tau*s
              rhodts     = rho*dtinv + s                       ! rho/dt + s

              c1         = t1st2*rhodts - bemol*gpdiv(igaus)   ! (1 - s*tau)*(rho/dt + s) - bemol*div(rho*a)
              c2         = ktau*rhodts                         ! k*tau*(rho/dt + s)
              c3         = rtau*rhodts                         ! rho*tau*(rho/dt + s)
              c4         = rho*t1st2                           ! rho*(1 - tau*s)
              c5         = rho*ktau                            ! rho*k*tau
              c6         = rtau*rho                            ! rho*rho*tau
              c7         = k*s*tau(igaus)                      ! k*s*tau
              c8         = -ktau*k                             ! -k^2*tau
              c9         = -c5                                 ! -rho*k*tau
              c10        = k                                   ! k
              c11        = 1.0_rp                              ! 1
              c12        = -bemol                              ! -bemol

              g2         = gprhs(1,igaus)
              g3         = gppro(igaus)
              g4         = g3 + tau(igaus)*g2              

              e1         = a1x*c6 
              e2         = a1y*c6 
              e3         = a1z*c6 
              axCF       = CF(igaus)*a1x
              ayCF       = CF(igaus)*a1y
              azCF       = CF(igaus)*a1z
              c10CD      = CD(igaus)

              A1(igaus)  = (c1)

              A2(igaus)  = -bemol*rho*ax + a2x*c3
              A3(igaus)  = -bemol*rho*ay + a2y*c3
              A4(igaus)  = -bemol*rho*az + a2z*c3

              A5(igaus)  =  bemo1*rho*ax - a1x*s*rtau
              A6(igaus)  =  bemo1*rho*ay - a1y*s*rtau
              A7(igaus)  =  bemo1*rho*az - a1z*s*rtau

              A8(igaus)  = a2x*e1 + c10CD  + axCF*a1x + kx
              A9(igaus)  = a2x*e2          + axCF*a1y
              A10(igaus) = a2x*e3          + axCF*a1z
              A11(igaus) = a2y*e1          + ayCF*a1x
              A12(igaus) = a2y*e2 + c10CD  + ayCF*a1y + ky
              A13(igaus) = a2y*e3          + ayCF*a1z
              A14(igaus) = a2z*e1          + azCF*a1x
              A15(igaus) = a2z*e2          + azCF*a1y
              A16(igaus) = a2z*e3 + c10CD  + azCF*a1z + kz


              B1(igaus)  = c2
              B2(igaus)  = c7
              B3(igaus)  = a1x*c5
              B4(igaus)  = a1y*c5
              B5(igaus)  = a1z*c5

              B6(igaus)  = c8
              B7(igaus)  = a2x*c9
              B8(igaus)  = a2y*c9
              B9(igaus)  = a2z*c9

              D1(igaus)  = g2-s*g4
              D2(igaus)  = a2x*rho*g4
              D3(igaus)  = a2y*rho*g4
              D4(igaus)  = a2z*rho*g4
              D5(igaus)  = k*g4

           end do

           do igaus=1,pgaus
              do inode=1,pnode
                 do idofn=1,ndofn
                    elrhs(idofn,inode)=elrhs(idofn,inode) + gpvol(igaus)   &
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

     else 

        !-------------------------------------------------------------------
        !
        ! OSS: Assembly
        !
        !-------------------------------------------------------------------

        if( ndime == 1 ) then

           do igaus=1,pgaus

              k          = gpdif(1,igaus)                        ! k
              s          = gprea(igaus)                        ! s
              rho        = gpden(igaus)                        ! rho
              div        = gpdiv(igaus)                        ! div(rho*a)
              a1x        = rho * gpvel(1,igaus)  

              g1         = rho * dtinv + s - bemol * div
              g2         = gprhs(1,igaus)
              g3         = gppro(igaus)
              axCF       = CF(igaus) * a1x

              A1(igaus)  =   g1
              A2(igaus)  = - bemol * a1x
              A5(igaus)  =   bemo1 * a1x 
              A8(igaus)  = k + a1x * a1x * tau(igaus) + CD(igaus)  + axCF * a1x
              D1(igaus)  = g2
              D2(igaus)  = a1x * g3

           end do

           do igaus = 1,pgaus
              do inode = 1,pnode
                 do idofn = 1,ndofn
                    elrhs(idofn,inode) = elrhs(idofn,inode) + gpvol(igaus)   &
                         &              *(  D1(igaus) * gpsha(inode,igaus)   &
                         &                + D2(igaus) * gpcar(1,inode,igaus) )

                 end do
              end do
           end do

           do igaus = 1,pgaus
              do inode = 1,pnode
                 do jnode = 1,pnode
                    elmat(inode,jnode) = elmat(inode,jnode) + gpvol(igaus) * (        &
                         & +gpsha(inode,igaus)*(    A1(igaus)*gpsha(jnode,igaus)      &  ! A1  Na Nb
                         &                       +  A5(igaus)*gpcar(1,jnode,igaus)  ) &  ! A5  Na dNb/dx
                         & +gpsha(jnode,igaus)*(    A2(igaus)*gpcar(1,inode,igaus)  ) &  ! A2  Nb dNa/dx
                         & +gpcar(1,inode,igaus)*(  A8(igaus)*gpcar(1,jnode,igaus)  ) )  ! A8  dNa/dx dNb/dx
                 end do
              end do
           end do

        else if( ndime == 2 ) then

           do igaus=1,pgaus

              k          = gpdif(1,igaus)                      ! k
              s          = gprea(igaus)                        ! s
              rho        = gpden(igaus)                        ! rho
              div        = gpdiv(igaus)                        ! div(rho*a)
              ax         = gpvel(1,igaus)  
              ay         = gpvel(2,igaus)
              rrtau      = rho * rho * tau(igaus) + CF(igaus) 
              rb         = rho * bemol
              rb1        = rho * bemo1
              g1         = rho * dtinv + s - bemol * div 
              g2         = gprhs(1,igaus)
              g3         = gppro(igaus) 
              CD(igaus)  = CD(igaus) + k

              A1(igaus)  =   g1
              A2(igaus)  = - rb * ax
              A3(igaus)  = - rb * ay

              A5(igaus)  =   rb1 * ax 
              A6(igaus)  =   rb1 * ay 

              A8(igaus)  = ax * ax * rrtau + CD(igaus) 
              A9(igaus)  = ax * ay * rrtau

              A11(igaus) = ax * ay * rrtau 
              A12(igaus) = ay * ay * rrtau + CD(igaus) 

              D1(igaus)  = g2
              D2(igaus)  = rho * ax * g3
              D3(igaus)  = rho * ay * g3

           end do

           do igaus = 1,pgaus
              do inode = 1,pnode
                 do idofn = 1,ndofn
                    elrhs(idofn,inode) = elrhs(idofn,inode) + gpvol(igaus)   &
                         &              *(  D1(igaus) * gpsha(inode,igaus)   &
                         &                + D2(igaus) * gpcar(1,inode,igaus) &
                         &                + D3(igaus) * gpcar(2,inode,igaus) )

                 end do
              end do
           end do

           do igaus = 1,pgaus
              do inode = 1,pnode
                 do jnode = 1,pnode
                    elmat(inode,jnode) = elmat(inode,jnode) + gpvol(igaus) * (        &
                         & +gpsha(inode,igaus)*(    A1(igaus)*gpsha(jnode,igaus)      &  ! A1  Na Nb
                         &                       +  A5(igaus)*gpcar(1,jnode,igaus)    &  ! A5  Na dNb/dx
                         &                       +  A6(igaus)*gpcar(2,jnode,igaus)  ) &  ! A6  Na dNb/dy
                         & +gpsha(jnode,igaus)*(    A2(igaus)*gpcar(1,inode,igaus)    &  ! A2  Nb dNa/dx
                         &                       +  A3(igaus)*gpcar(2,inode,igaus)  ) &  ! A3  Nb dNa/dy
                         & +gpcar(1,inode,igaus)*(  A8(igaus)*gpcar(1,jnode,igaus)    &  ! A8  dNa/dx dNb/dx
                         &                       +  A9(igaus)*gpcar(2,jnode,igaus)  ) &  ! A9  dNa/dx dNb/dy
                         & +gpcar(2,inode,igaus)*( A11(igaus)*gpcar(1,jnode,igaus)    &  ! A11 dNa/dy dNb/dx
                         &                       + A12(igaus)*gpcar(2,jnode,igaus)  ) )  ! A12 dNa/dy dNb/dy
                 end do
              end do
           end do

        else if( ndime == 3 ) then

           do igaus=1,pgaus

              k          = gpdif(1,igaus)                      ! k
              s          = gprea(igaus)                        ! s
              rho        = gpden(igaus)                        ! rho
              div        = gpdiv(igaus)                        ! div(rho*a)
              ax         = gpvel(1,igaus)  
              ay         = gpvel(2,igaus)
              az         = gpvel(3,igaus)
              rrtau      = rho * rho * tau(igaus) + CF(igaus) 
              rb         = rho * bemol
              rb1        = rho * bemo1
              g1         = rho * dtinv + s - bemol * div 
              g2         = gprhs(1,igaus)
              g3         = gppro(igaus)              
              CD(igaus)  = CD(igaus) + k

              A1(igaus)  =   g1
              A2(igaus)  = - rb * ax
              A3(igaus)  = - rb * ay
              A4(igaus)  = - rb * az

              A5(igaus)  =   rb1 * ax 
              A6(igaus)  =   rb1 * ay 
              A7(igaus)  =   rb1 * az 

              A8(igaus)  = ax * ax * rrtau + CD(igaus) 
              A9(igaus)  = ax * ay * rrtau
              A10(igaus) = ax * az * rrtau

              A11(igaus) = ay * ax * rrtau 
              A12(igaus) = ay * ay * rrtau + CD(igaus) 
              A13(igaus) = ay * az * rrtau 

              A14(igaus) = az * ax * rrtau 
              A15(igaus) = az * ay * rrtau 
              A16(igaus) = az * az * rrtau + CD(igaus) 

              D1(igaus)  = g2
              D2(igaus)  = rho * ax * g3
              D3(igaus)  = rho * ay * g3
              D4(igaus)  = rho * az * g3

           end do

           do igaus = 1,pgaus
              do inode = 1,pnode
                 do idofn = 1,ndofn
                    elrhs(idofn,inode) = elrhs(idofn,inode) + gpvol(igaus)   &
                         &              *(  D1(igaus) * gpsha(inode,igaus)   &
                         &                + D2(igaus) * gpcar(1,inode,igaus) &
                         &                + D3(igaus) * gpcar(2,inode,igaus) &
                         &                + D4(igaus) * gpcar(3,inode,igaus) )

                 end do
              end do
           end do

           do igaus=1,pgaus
              do inode=1,pnode
                 do jnode=1,pnode
                    elmat(inode,jnode)=elmat(inode,jnode)+gpvol(igaus)*(&
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
                         &                       + A16(igaus)*gpcar(3,jnode,igaus)  ) )
                 end do
              end do
           end do
        end if

     end if
     !
     ! Add lumped reaction
     !
     if( ireac == 1 ) then
        call elmlum(pnode,pgaus,pelty,elcod,gpre2,gpcar,elmat)
     end if
     !
     ! Extension elements
     !
     if( lelch(ielem) == ELEXT ) then
        call elmext(&
             4_ip,1_ip,pnode,dummr,dummr,dummr,dummr,elmat,&
             dummr,elrhs,elunk)
     end if
     
  else if ( itask == 4 ) then

     !----------------------------------------------------------------------
     !
     ! Projection: 
     ! 1. - tau1' * R(u) 
     ! 2.   tau1' * rho * [ a.grad(u) ]
     ! 
     !----------------------------------------------------------------------

     if( kfl_stabi >= 1 ) then

        do inode = 1,pnode
           elrhs(1,inode) = 0.0_rp
        end do

        if( kfl_stabi == 1 ) then

           do igaus = 1,pgaus
              gpres  =  -gprhs(1,igaus) 
              c1     =   gpden(igaus) * dtinv + gprea(igaus)
              do inode=1,pnode
                 gpres = gpres + elunk(inode,1)*(&
                      &  + c1 * gpsha(inode,igaus)   &
                      &  - gpdif(1,igaus) * gplap(inode,igaus))
                 do idime = 1,ndime
                    gpres = gpres + elunk(inode,1) * &
                         ( gpden(igaus) * gpvel(idime,igaus) - gpgrd(idime,igaus) ) &
                         * gpcar(idime,inode,igaus)
                 end do
              end do
              gpres = gpvol(igaus) * tau(igaus) * gpres
              do inode = 1,pnode
                 elrhs(1,inode) = elrhs(1,inode) + gpres * gpsha(inode,igaus) 
              end do
           end do

        else 

           do igaus = 1,pgaus
              gpres =  0.0_rp
              do inode = 1,pnode
                 c1 = 0.0_rp
                 do idime = 1,ndime
                    c1 = c1 + gpvel(idime,igaus) * gpcar(idime,inode,igaus)
                 end do
                 gpres = gpres + c1 * elunk(inode,1)
              end do
              gpres = gpvol(igaus) * tau(igaus) * gpden(igaus) * gpres
              do inode = 1,pnode
                 elrhs(1,inode) = elrhs(1,inode) + gpres * gpsha(inode,igaus)
              end do
           end do

        end if
        !
        ! Extension elements
        !
        if( lelch(ielem) == ELEXT ) then
           call elmext(&
                3_ip,1_ip,pnode,dummr,dummr,dummr,dummr,elmat,&
                dummr,elrhs,elunk)
        end if
        !
        ! Assemble
        !
        call assrhs(1_ip,pnode,lnods,elrhs,rhsid)

     end if

     !----------------------------------------------------------------------
     !
     ! Update SGS:
     ! 1. u' = tau1' * [ u'n/dt + R(u) ] + Pi 
     ! 2. u' = tau1' * [ u'n/dt - rho * a.grad(u) ] + Pi 
     !
     ! The temporal term was already added at the beginning
     ! 
     !----------------------------------------------------------------------

     if( kfl_sgsti /= 0 .or. kfl_sgsco /= 0 ) then

        do igaus = 1,pgaus
           gpsgs(igaus,1) = 0.0_rp
        end do

        if( kfl_stabi >= 2 ) then
           do igaus = 1,pgaus
              do inode = 1,pnode
                 do idime = 1,ndime
                    c1 =  gpden(igaus) * gpadv(idime,igaus) 
                    gpsgs(igaus,1) = gpsgs(igaus,1) &
                         - c1 * elunk(inode,1) * gpcar(idime,inode,igaus)
                 end do
              end do
           end do
        else
           do igaus = 1,pgaus
              gpsgs(igaus,1) = gpsgs(igaus,1) + gprhs(1,igaus) 
              c1 = - gpden(igaus) * dtinv - gprea(igaus)
              do inode=1,pnode
                 gpsgs(igaus,1) = gpsgs(igaus,1) + elunk(inode,1)*(&
                      &  + c1 * gpsha(inode,igaus)   &
                      &  + gpdif(1,igaus) * gplap(inode,igaus))
                 do idime = 1,ndime
                    gpsgs(igaus,1) = gpsgs(igaus,1) - elunk(inode,1) * &
                         gpden(igaus) * gpadv(idime,igaus) * gpcar(idime,inode,igaus)
                 end do
              end do
           end do
        end if
        do igaus = 1,pgaus
           gpsgs(igaus,1)  = tau(igaus) * gpsgs(igaus,1) 
        end do

        if( kfl_stabi >= 1 ) then
           !
           ! Add projection (+ time derivative of SGS)
           !
           do igaus = 1,pgaus
              gpsgs(igaus,1)  = gpsgs(igaus,1) + gppro(igaus) 
           end do
        end if

     end if

  else if( itask == 5 .and. kfl_limit /= 0 ) then

     !-------------------------------------------------------------------
     !
     ! Postprocess limiter
     ! alpha = 0 in sharp region
     !       = 1 in smooth region
     !
     !-------------------------------------------------------------------

     do inode = 1,pnode
        elrhs(1,inode) = 0.0_rp
     end do

     do igaus = 1,pgaus
        c1 = 0.0_rp
        c2 = gppro(igaus)
        do inode = 1,pnode
           c4 = 0.0_rp
           do idime = 1,ndime
              c4 = c4 + gpvel(idime,igaus) * gpcar(idime,inode,igaus)
           end do
           c1 = c1 + c4 * elunk(inode,1)
        end do
        c1 = gpden(igaus) * tau(igaus) * c1 
        c3 = abs(c1) + abs(c2) 
        if( c3 > 1.0e-6_rp ) then
           beta = abs( c1 - c2 ) / c3
        else
           beta = 0.0_rp
        end if
        if( kfl_limit == 1 ) then
           alpha = min(1.0_rp,2.0_rp*(1.0_rp-beta))
        else if( kfl_limit == 2 ) then
           alpha = 0.5_rp*(tanh(20.0_rp*(beta-0.8_rp))+1.0_rp)
        else
           alpha = 0.0_rp
        end if
        do inode = 1,pnode
           elrhs(1,inode) =  elrhs(1,inode) + alpha * gpvol(igaus) * gpsha(inode,igaus)
        end do
     end do
     !
     ! Extension elements
     !
     if( lelch(ielem) == ELEXT ) then
        call elmext(&
             3_ip,1_ip,pnode,dummr,dummr,dummr,dummr,elmat,&
             dummr,elrhs,elunk)
     end if
     !
     ! Assemble
     !
     call assrhs(1_ip,pnode,lnods,elrhs,rhsid)

  else if( itask == 6 .and. kfl_sgsti /= 0 ) then

     !-------------------------------------------------------------------
     !
     ! Postprocess SGS
     !
     !-------------------------------------------------------------------

     do inode = 1,pnode
        elrhs(1,inode) = 0.0_rp
     end do

     do igaus = 1,pgaus
        do inode = 1,pnode
           elrhs(1,inode) =  elrhs(1,inode) + gpvol(igaus) * gpsha(inode,igaus) * gpsgs(igaus,1)
        end do
     end do
     !
     ! Extension elements
     !
     if( lelch(ielem) == ELEXT ) then
        call elmext(&
             3_ip,1_ip,pnode,dummr,dummr,dummr,dummr,elmat,&
             dummr,elrhs,elunk)
     end if
     !
     ! Assemble
     !
     call assrhs(1_ip,pnode,lnods,elrhs,rhsid)

  else if( itask == 7 ) then

     !-------------------------------------------------------------------
     !
     ! Postprocess Projection
     !
     !-------------------------------------------------------------------

     do inode = 1,pnode
        elrhs(1,inode) = 0.0_rp
     end do

     do igaus = 1,pgaus
        do inode = 1,pnode
           elrhs(1,inode) =  elrhs(1,inode) + gpvol(igaus) * gpsha(inode,igaus) * gppro(igaus)
        end do
     end do
     !
     ! Extension elements
     !
     if( lelch(ielem) == ELEXT ) then
        call elmext(&
             3_ip,1_ip,pnode,dummr,dummr,dummr,dummr,elmat,&
             dummr,elrhs,elunk)
     end if
     !
     ! Assemble
     !
     call assrhs(1_ip,pnode,lnods,elrhs,rhsid)

  end if

end subroutine elmadr
