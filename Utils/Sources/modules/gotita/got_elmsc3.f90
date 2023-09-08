subroutine got_elmsc3(&
     pnode,pgaus,pevat,ndofn,gpgcd,gpvdr,gpcdr,gpst2,&
     gprhs,resic,gpvol,gpcar,elcdr,dtinv_got,chale,elmat)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmshc
  ! NAME
  !   got_elmshc
  ! DESCRIPTION
  !   This routine computes the contribution from the shock capturing
  !   term to the turbulence equation      
  ! OUTPUT
  !    ELMAT ... Element matrix
  ! USES
  ! USED BY
  !    got_elmope
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_gotita, only       :  kfl_shocc_got,ndofn_got,shock_got,&
       &                        kfl_timec_got,kfl_weigc_got
  implicit none  
  integer(ip), intent(in)    :: pnode,pgaus,pevat,ndofn
  real(rp),    intent(in)    :: gpgcd(ndime,pgaus),gpvdr(ndime,pgaus)
  real(rp),    intent(in)    :: gpcdr(pgaus,2),gpst2(pgaus)
  real(rp),    intent(in)    :: gprhs(pgaus),resic(pnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus),gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: elcdr(pnode),dtinv_got,chale(*)
  real(rp),    intent(inout) :: elmat(pevat,pevat)                    
  real(rp)                   :: scdif,sldif,gpres
  real(rp)                   :: gpve2,grnor,addit,gpv12,fact1
  integer(ip)                :: idime,inode,jnode,kdime,ldime,igaus

  do igaus=1,pgaus
     !
     ! Continuity residual |R(u)|
     ! |R(u)| = |f - alpha/(theta*dt)-div(alpha*u)|
     ! f = alpha^n/(theta*dt)
     !
     gpres=gprhs(igaus)
     do inode=1,pnode
        gpres=gpres-resic(inode,igaus)*elcdr(inode)
     end do
     if(kfl_timec_got==1.and.kfl_weigc_got==0)&
          gpres=gpres-dtinv_got*(gpcdr(igaus,1)-gpcdr(igaus,2))
     gpres=abs(gpres)
     !
     ! Square velocity norm u^2
     !
     gpve2=0.0_rp
     do idime=1,ndime
        gpve2=gpve2+gpvdr(idime,igaus)*gpvdr(idime,igaus)
     end do
     !
     ! Alpha gradient norm |grad(alpha)|
     !
     call vecnor(gpgcd(1,igaus),ndime,grnor,2_ip) 

     if(gpve2>0.0_rp.and.grnor>1.0e-6) then

        scdif = 0.5_rp*shock_got*chale(2)*gpres/grnor          ! scdif = kiso = 1/2*C*h*[|R|/|grad(alpha)|]
        !
        ! For P1=1.5. For Q1=0.7
        !
        if((pnode>3.and.ndime==2).or.(pnode>4.and.ndime==3)) then
           !scdif=1.5_rp*scdif
        else
           !scdif=0.7_rp*scdif        
        end if
        !
        ! Compute diffusion introduced along the streamlines
        !
        if(kfl_shocc_got==1) then                              ! Isotropic SC:
           sldif=scdif                                         ! sldif = kiso
        else                                                   ! Anisotropic SC:
           sldif=max(0.0_rp,scdif-gpst2(igaus)*gpve2)          ! sldif = <kiso-k'> = max(0,kiso-(tau)*a^2)
        endif
        fact1=(sldif-scdif)/gpve2

        do inode=1,pnode
           do jnode=1,pnode 
              addit=0.0_rp
              do kdime=1,ndime
                 addit=addit+scdif*gpcar(kdime,inode,igaus)&   ! kiso*grad(Ni).grad(Nj)
                      *gpcar(kdime,jnode,igaus)                 
                 do ldime=1,ndime                              ! (<kiso-k'>-kiso)/u^2*
                    addit=addit+fact1&                         ! grad(Ni).(u x u).grad(Nj)
                         *gpcar(kdime,inode,igaus)&
                         *gpvdr(kdime,igaus)&
                         *gpvdr(ldime,igaus)&
                         *gpcar(ldime,jnode,igaus)
                 end do
              end do
              elmat(inode,jnode)=elmat(inode,jnode)&
                   +addit*gpvol(igaus)
           end do
        end do
     end if
  end do

end subroutine got_elmsc3
!------------------------------------------------------------------------
! NOTES
!
! Shock capturing for continuity equation:
! 
! d alpha
! ------- + div(alpha*u) = 0
!    dt
!
! R    = Residual of the equation   [kg/(m^3*s)]
! tau  = Stabilization parameter    [s]
! kiso = Isotropic SC diffusion     [m^2/s]
! k'   = Anisotropic SC diffusion   [m^2/s]
! C    = Shock capturing constant   
!
!        1           |R|
! kiso = - C*h  ------------- , k'=tau_2*u^2
!        2      |grad(alpha)|
!
! Isotropic:     kiso*grad(alpha).grad(v)  
!                                             u x u
! Anisotropic:   (<kiso-k'>-kiso)*grad(alpha).-----.grad(v)
!                                              u^2
!***
!------------------------------------------------------------------------
