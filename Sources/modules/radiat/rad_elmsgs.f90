subroutine rad_elmsgs(&
     pnode,pgaus,gpvel,gpdif,gprea,gpden,gpcod,chale,&
     rtemp,gprhs,gpgrt,eltem,gpsta,gpstt,gpstp,gpsgs,&
     gpres)
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_elmsgs
  ! NAME
  !   rad_elmsgs
  ! DESCRIPTION
  !    Compute the stabilization parameters, the subgrid scale and the
  !    residual
  ! OUTPUT
  !    GPSTP ... Stabilization parameter tau' [T/rho]=[T*L^3/M]
  !    GPSTT ... tau'*tau^{-1}
  !    GPSGS ... T'=tau'*(R(u)+r)
  !              R(u)=f-rho*u/(theta*dt)-rho*a.grad(u)+div[k*grad(u)]-s*u
  !              f=Q+rho*u^n/(theta*dt)
  !              r=u'^n/(theta*dt)
  !    GPRES ... Coarse residual R(u)+r
  ! USES
  ! USED BY
  !    rad_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mgaus,kfl_naxis
  use def_radiat, only     :  kfl_sgsno_rad,kfl_taust_rad,&
       &                      kfl_bubbl_rad,dtinv_rad,kfl_ellen_rad,&
       &                      staco_rad
  use mod_tauadr, only       :  tauadr
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  real(rp),    intent(in)  :: gpvel(ndime,pgaus)
  real(rp),    intent(in)  :: gpdif(pgaus),gprea(pgaus),gpden(pgaus)
  real(rp),    intent(in)  :: gpcod(ndime,pgaus)
  real(rp),    intent(in)  :: eltem(pnode,*)
  real(rp),    intent(in)  :: chale(2)
  real(rp),    intent(in)  :: rtemp(pnode,pgaus),gprhs(pgaus)
  real(rp),    intent(in)  :: gpgrt(ndime,pgaus)
  real(rp),    intent(out) :: gpsta(pgaus),gpstt(pgaus),gpstp(pgaus)
  real(rp),    intent(out) :: gpsgs(pgaus,2),gpres(mgaus)
  integer(ip)              :: igaus,inode,idime
  real(rp)                 :: gpnve,gpmve(3)

  if(kfl_taust_rad==0) then
     !
     ! No stabilization of Bubble
     !
     return
     gpsta=0.0_rp
     gpstp=0.0_rp
     gpstt=1.0_rp

  else if(kfl_bubbl_rad==1) then
     !
     ! Bubble: make TTEMP to be the adjoint
     !
     gpsta= 0.0_rp
     gpstt= 0.0_rp
     gpstp=-1.0_rp

  else
     !
     ! Stabilization parameter GPSTA=tau
     ! 
     do igaus=1,pgaus

        if(kfl_ellen_rad==-1) then
           !
           ! a_eff=rho*Cp*|a.grad(T)|/|grad(T)|
           !
           gpmve(1) = 0.0_rp
           gpnve    = 0.0_rp
           do idime=1,ndime
              gpmve(1) = gpmve(1) + gpvel(idime,igaus)*gpgrt(idime,igaus)
              gpnve    = gpnve    + gpgrt(idime,igaus)*gpgrt(idime,igaus)
           end do
           if(gpnve/=0.0_rp) gpnve=gpden(igaus)*abs(gpmve(1))/sqrt(gpnve)
        else
           if(kfl_naxis==0) then
              call vecnor(gpvel(1,igaus),ndime,gpnve,2_ip) 
              gpnve=gpden(igaus)*gpnve
           else
              do idime=1,ndime
                 gpmve(idime)=gpden(igaus)*gpvel(idime,igaus)
              end do
              gpmve(1)=gpmve(1)-gpdif(igaus)/gpcod(1,igaus)
              call vecnor(gpmve,ndime,gpnve,2_ip) 
           end if           
        end if
        
        call tauadr(&
             kfl_taust_rad,staco_rad,gpnve,gpdif(igaus),gprea(igaus),&
             chale(1),chale(2),gpsta(igaus))
     end do
 
     if(kfl_sgsno_rad==1) then
        !
        ! Residual GPRES=R(u)+r
        !
        do igaus=1,pgaus
           gpres(igaus)=gprhs(igaus)
           do inode=1,pnode
              gpres(igaus)=gpres(igaus)&
                   -rtemp(inode,igaus)*eltem(inode,1)
           end do
        end do

        ! 
        ! Subgrid scale GPSGS=tau*(R(u)+r)
        !
        do igaus=1,pgaus
           gpsgs(igaus,1)=gpsta(igaus)*gpres(igaus)
           gpstp(igaus)=gpsta(igaus)
           gpstt(igaus)=1.0_rp 
        end do

     else

        do igaus=1,pgaus
           gpstp(igaus)=gpsta(igaus)
           gpstt(igaus)=1.0_rp 
        end do
     end if
  end if

end subroutine rad_elmsgs
