subroutine tem_elmpre(&
     ielem,pnode,pgaus,pmate,gpden,gpsph,gpsgv,&
     gpsha,gpcar,gphes,elvel,eltem,elcod,elmsh,&
     gpvel,gptem,gprhs,gpcod,gpgrt,lnods,gpmsh)
  !-----------------------------------------------------------------------
  !****f* Temper/tem_elmpre
  ! NAME
  !   tem_elmpre
  ! DESCRIPTION
  !    Compute some Gauss values
  ! OUTPUT 
  !    GPGRD(NDIME) ... grad(k) coefficient
  !    GPTEM .......... Temperature of previous iterations and time steps
  !    GPVEL .......... Advection a
  !    GPADV .......... Advection term a.grad(Ni)
  !    GPRGD .......... Thermal conductivity gradient grad(k+kt)
  ! USES
  ! USED BY
  !    tem_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  dpthe,cutim,kfl_coupl,kcond,vesgs,velom
  use def_domain, only       :  mnode,mgaus,ndime,ntens,kfl_naxis,&
       &                        xfiel
  use def_temper, only       :  kfl_advec_tem,ADR_tem,&
       &                        kfl_sourc_tem,kfl_sgsve_tem,&
       &                        kfl_exacs_tem,kfl_regim_tem,&
       &                        prtur_tem,&
       &                        kfl_ellen_tem,grtem_tem
  use mod_ADR,    only       :  BDF
  use def_kermod, only       :  gasco 
  use mod_ker_space_time_function
  implicit none 
  integer(ip), intent(in)    :: ielem
  integer(ip), intent(in)    :: pnode,pgaus,pmate,lnods(pnode)
  real(rp),    intent(in)    :: gpsph(pgaus),gpsgv(ndime,pgaus)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)    :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(in)    :: elvel(ndime,pnode)
  real(rp),    intent(in)    :: eltem(pnode,*)
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(in)    :: elmsh(ndime,pnode)
  real(rp),    intent(inout) :: gpden(pgaus)
  real(rp),    intent(out)   :: gpvel(ndime,pgaus),gpcod(ndime,pgaus)
  real(rp),    intent(out)   :: gprhs(pgaus)
  real(rp),    intent(out)   :: gptem(pgaus,*)
  real(rp),    intent(out)   :: gpgrt(ndime,pgaus)
  real(rp),    intent(out)   :: gpmsh(ndime,pgaus)
  integer(ip)                :: idime,inode,igaus,itime,ipoin,mxdim
  real(rp)                   :: dummr, fact0, ugrat
  !
  ! Temperature: GPTEM
  !
  if( ADR_tem % kfl_time_integration /= 0 ) then
     do igaus = 1,pgaus
        gptem(igaus,2) = 0.0_rp
        do inode = 1,pnode
           gptem(igaus,2) = gptem(igaus,2) + gpsha(inode,igaus) * eltem(inode,2)
        end do
     end do
     if( ADR_tem % kfl_time_scheme == BDF ) then
        do itime = 3,ADR_tem % kfl_time_order + 1
           do igaus = 1,pgaus
              gptem(igaus,itime) = 0.0_rp
              do inode = 1,pnode
                 gptem(igaus,itime) = gptem(igaus,itime) &
                      + eltem(inode,itime) * gpsha(inode,igaus)
              end do
           end do
        end do
     end if
  end if
  !
  ! Density: GPDEN=rho*Cp
  !
  if( kfl_regim_tem == 1 ) then
     do igaus = 1,pgaus
        gpden(igaus) = gpden(igaus) * (gpsph(igaus)-gasco)     ! rho*Cv=rho*(Cp-R)
     end do
  else if (kfl_regim_tem /= 4) then
     do igaus = 1,pgaus
        gpden(igaus) = gpden(igaus) * gpsph(igaus)             ! rho*Cp
     end do
  end if
  !
  ! Coordinates
  !
  do igaus = 1,pgaus
     do idime = 1,ndime
        gpcod(idime,igaus) = 0.0_rp
     end do
     do inode = 1,pnode
        do idime = 1,ndime
           gpcod(idime,igaus) = gpcod(idime,igaus)&
                + gpsha(inode,igaus) * elcod(idime,inode)
        end do
     end do
  end do
  !
  ! Velocity GPVEL=a and advection term GPADV=a.grad(Ni)-(k/r,0).grad(Ni)
  !
  do igaus = 1,pgaus
     do idime = 1,ndime
        gpvel(idime,igaus) = 0.0_rp
     end do
  end do

  if( kfl_advec_tem /= 0 ) then

     if( kfl_advec_tem == 1 ) then
        do igaus = 1,pgaus
           do inode = 1,pnode
              do idime = 1,ndime
                 gpvel(idime,igaus) = gpvel(idime,igaus) &
                      + gpsha(inode,igaus) * elvel(idime,inode)
              end do
           end do
        end do
        if( kfl_sgsve_tem == 1 ) then 
           do igaus = 1,pgaus
              do idime = 1,ndime
                 gpvel(idime,igaus) = gpvel(idime,igaus) +vesgs(ielem) % a(idime,igaus,1)
              end do
           end do
        end if

     else if( kfl_advec_tem >= 2 ) then
        call tem_velfun(pgaus,gpcod,gpvel)

     else if( kfl_advec_tem == -1 ) then
        do idime = 1,ndime
           gpgrt(idime,igaus) = 0.0_rp
        end do
        do inode = 1,pnode
           do idime = 1,ndime
              gpgrt(idime,igaus) = gpgrt(idime,igaus)&
                   + gpcar(idime,inode,igaus)*eltem(inode,1)
           end do
        end do
        dummr = 0.0_rp
        do idime = 1,ndime
           dummr = dummr +gpgrt(idime,igaus) * gpgrt(idime,igaus)
        end do
        dummr = sqrt(dummr)
        if( dummr /= 0.0_rp ) then
           do idime = 1,ndime
              gpvel(idime,igaus) = -gpgrt(idime,igaus) / dummr
           end do
        end if
     end if

     if( associated(velom) ) then
        gpmsh = 0.0_rp
        do igaus = 1,pgaus
           do inode = 1,pnode
              do idime = 1,ndime
                 gpmsh(idime,igaus) = gpmsh(idime,igaus)&
                      + gpsha(inode,igaus) * elmsh(idime,inode)
              end do
           end do
           gpvel(1:ndime,igaus) = gpvel(1:ndime,igaus) - gpmsh(1:ndime,igaus)
        end do        
     end if
     
  end if
  !
  ! Source term: GPRHS
  !
  if( kfl_sourc_tem > 0 ) then
     mxdim = min(2_ip,ndime)
     do igaus = 1,pgaus
        call ker_space_time_function(&
             kfl_sourc_tem,gpcod(1,igaus),gpcod(mxdim,igaus),gpcod(ndime,igaus),cutim,gprhs(igaus))
     end do
  else if( kfl_sourc_tem < 0 ) then
     do igaus = 1,pgaus
        gprhs(igaus) = xfiel(-kfl_sourc_tem) % a(1,ielem,1)
     end do
  else
     do igaus = 1,pgaus
        gprhs(igaus) = 0.0_rp     
     end do
  end if
  !
  ! Low-Mach: dp0/dt
  ! prefactor alpha * T ~ 1 
  !
  if( kfl_regim_tem >= 3 ) then
     do igaus = 1,pgaus
        gprhs(igaus) = gprhs(igaus) + dpthe 
     end do

  end if
  !
  ! GPGRT: grad(T)
  !
  if( kfl_ellen_tem == -1 ) then
     do igaus = 1,pgaus
        do idime = 1,ndime
           gpgrt(idime,igaus) = 0.0_rp
        end do
        do inode = 1,pnode
           ipoin = lnods(inode)
           do idime = 1,ndime
              gpgrt(idime,igaus) = gpgrt(idime,igaus) &
                   + gpcar(idime,inode,igaus) * eltem(inode,1)
           end do
        end do
     end do
  end if

end subroutine tem_elmpre
