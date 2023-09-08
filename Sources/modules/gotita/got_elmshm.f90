subroutine got_elmshm(&
     pnode,pgaus,pevat,ndofn,gpgvd,gpvdr,gpcdr,gpst1,&
     gprhs,resim,gpvol,gpcar,elcdr,elvdr,dtinv_got,&
     chale,elmat)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmshc
  ! NAME
  !   got_elmshc
  ! DESCRIPTION
  !   This routine computes the contribution from the shock capturing
  !   term to the turbulence equation. 
  ! OUTPUT
  !    ELMAT ... Element matrix
  ! USES
  ! USED BY
  !    got_elmope
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mgaus,mnode
  use def_gotita, only       :  kfl_shocm_got,ndofn_got,shock_got,&
       &                        kfl_probl_got,kfl_coupl_got,penal_got,&
       &                        kfl_timem_got,kfl_weigm_got
  implicit none  
  integer(ip), intent(in)    :: pnode,pgaus,pevat,ndofn
  real(rp),    intent(in)    :: gpgvd(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpvdr(ndime,pgaus,2),gpcdr(pgaus)
  real(rp),    intent(in)    :: gpst1(pgaus)
  real(rp),    intent(in)    :: gprhs(ndofn_got(3),pgaus)
  real(rp),    intent(in)    :: resim(ndofn_got(3),pnode,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: elvdr(ndime,pnode),elcdr(pnode)
  real(rp),    intent(in)    :: dtinv_got,chale(*)
  real(rp),    intent(inout) :: elmat(pevat,pevat)                    
  real(rp)                   :: scdif,sldif,gpres
  real(rp)                   :: gpve2,grnor,addit,gpv12
  integer(ip)                :: idime,inode,jnode,kdime,ldime,igaus
  integer(ip)                :: idofn,jdofn,mdime

  do igaus=1,pgaus
     !
     ! m-momentum residual |R(u)|
     !
     do mdime=1,ndime
        gpres=gprhs(mdime,igaus)
        do inode=1,pnode
           gpres=gpres&
                -resim(      1,inode,igaus)*elvdr(mdime,inode)&
                -resim(1+mdime,inode,igaus)*elcdr(inode)
        end do
        if(kfl_timem_got==1.and.kfl_weigm_got==0) then
           if(kfl_coupl_got==0) then
              gpres=gpres-dtinv_got*&
                   (gpvdr(mdime,igaus,1)-gpvdr(mdime,igaus,2))
           else
              gpres=gpres-dtinv_got*(gpcdr(igaus)+penal_got)&
                   *(gpvdr(mdime,igaus,1)-gpvdr(mdime,igaus,2))
           end if
        end if
        gpres=abs(gpres)
        !
        ! Square velocity norm u^2 and 
        !
        gpve2=0.0_rp
        do idime=1,ndime
           gpve2=gpve2+gpvdr(idime,igaus,1)*gpvdr(idime,igaus,1)
        end do
        !
        ! Velocity m-component gradient norm |grad(um)|
        !
        call vecnor(gpgvd(1,mdime,igaus),ndime,grnor,2_ip) 

        if(gpve2>0.0_rp.and.grnor>1.0e-6) then

           gpv12 = 1.0_rp/gpve2
           scdif = 0.5_rp*shock_got*chale(2)*gpres/grnor    ! scdif = kiso = 1/2*C*h*[|R|/|grad(um)|]
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
           if(kfl_shocm_got==1) then                       ! Isotropic SC:
              sldif=scdif                                  ! sldif = kiso
           else                                            ! Anisotropic SC:
              if(kfl_coupl_got==1.and.kfl_probl_got==1) then
                 sldif=max(0.0_rp,scdif-gpst1(igaus)&      ! sldif = <kiso-k'> = max(0,kiso-alpha^2*(tau_1)*u^2)
                      *gpve2*(gpcdr(igaus)+penal_got))     !                   = max(0,kiso-alpha*(tau_1')*u^2)
              else
                 sldif=max(0.0_rp,scdif-gpst1(igaus)*gpve2)
              end if
           end if

           do inode=1,pnode
              do jnode=1,pnode 
                 addit=0.0_rp
                 do kdime=1,ndime
                    addit=addit+scdif&                     ! kiso*grad(Ni).grad(Nj)
                         *gpcar(kdime,inode,igaus)&   
                         *gpcar(kdime,jnode,igaus)                 
                    do ldime=1,ndime                       ! (<kiso-k'>-kiso)/u^2*
                       addit=addit+(sldif-scdif)&          ! grad(Ni).(u x u).grad(Nj)
                            *gpcar(kdime,inode,igaus)&
                            *gpvdr(kdime,igaus,1)&
                            *gpvdr(ldime,igaus,1)*gpv12&
                            *gpcar(ldime,jnode,igaus)
                    end do
                 end do
                 idofn=(inode-1)*ndofn+mdime
                 jdofn=(jnode-1)*ndofn+mdime
                 elmat(idofn,jdofn)=elmat(idofn,jdofn)&
                      +addit*gpvol(igaus)
              end do
           end do
        end if
     end do
  end do

end subroutine got_elmshm
!------------------------------------------------------------------------
! NOTES
!
! Shock capturing for the m-momentum equations (m=1,2,3):
! 
! Coupled equations:
!       du
! alpha -- + alpha*(u.grad)u + sig*alpha*(u-ua) = 0
!       dt
!
! For uncoupled equations, take alpha=1.
!
! R    = Residual of the equation   [kg/(m^3*s)]
! tau  = Stabilization parameter    [s]
! kiso = Isotropic SC diffusion     [m^2/s]
! k'   = Anisotropic SC diffusion   [m^2/s]
! C    = Shock capturing constant   
!
!        1         |R|
! kiso = - C*h  ---------- , k'=alpha^2*tau_1*u^2=alpha*tau_1'*u^2
!        2      |grad(um)|
!
! Isotropic:     kiso*grad(um).grad(v)  
!                                         u x u
! Anisotropic:   (<kiso-k'>-kiso)*grad(um).-----.grad(v)
!                                          u^2
!***
!------------------------------------------------------------------------
