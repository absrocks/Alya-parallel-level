subroutine nsi_elmpre(&
     pmate,pnode,pgaus,plapl,gpsgt,gpsha,gpcar,gphes,&
     elcod,elpre,elvel,eltem,elkin,elden,elfle,gppre,&
     gpcod,gpvel,gptem,gpgve,gpgpr,gpgrk,gpgde,gpadv,&
     gplap,gpgv2,gpgrt,gpdiv,gpden,gpfle,gpsgs)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmpre
  ! NAME 
  !    nsi_elmpre
  ! DESCRIPTION
  !    Compute some variables at the Gauss points
  ! OUTPUT
  !    GPPRE ... Pressure
  !    GPCOD ... Coordinates
  !    GPVEL ... Velocity
  !    GPTEM ... Temperature
  !    GPVEL ... Velocity gradients
  !    GPGPR ... Pressure gradients
  !    GPGPT ... Temperature gradients
  !    GPGRK ... Turbulence kinetic energy gradients
  !    GPGV2 ... Velocity second order derivative
  ! USES
  ! USED BY
  !    nsi_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mgaus,mnode,ntens
  use def_nastin, only     :  kfl_penal_nsi,kfl_timei_nsi,kfl_linea_nsi,&
       &                      kfl_cotem_nsi,kfl_advec_nsi,kfl_tiacc_nsi,&
       &                      kfl_tisch_nsi,kfl_predi_nsi,kfl_sgste_nsi,&
       &                      kfl_grtur_nsi,kfl_cotur_nsi,kfl_regim_nsi,&
       &                      kfl_visco_nsi,kfl_sgsti_nsi,&
       &                      kfl_surte_nsi,kfl_sgsco_nsi,&
       &                      nbdfp_nsi
  implicit none
  integer(ip), intent(in)  :: pmate,pnode,pgaus,plapl
  real(rp),    intent(in)  :: gpsgt(pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(in)  :: elcod(ndime,pnode),elpre(pnode,*)
  real(rp),    intent(in)  :: elvel(ndime,mnode,*),eltem(pnode,*)
  real(rp),    intent(in)  :: elkin(pnode),elden(pnode,*)
  real(rp),    intent(in)  :: elfle(pnode)
  real(rp),    intent(out) :: gpcod(ndime,pgaus),gppre(pgaus,2)
  real(rp),    intent(out) :: gpvel(ndime,pgaus,*),gptem(pgaus,*)
  real(rp),    intent(out) :: gpgve(ndime,ndime,pgaus) 
  real(rp),    intent(out) :: gpgv2(ntens,ndime,pgaus)
  real(rp),    intent(out) :: gpgrt(ndime,pgaus),gpdiv(pgaus)
  real(rp),    intent(out) :: gpden(pgaus,*),gpfle(pgaus)
  real(rp),    intent(in)  :: gpsgs(ndime,pgaus)
  real(rp),    intent(out) :: gpgpr(ndime,pgaus,2),gpgrk(ndime,pgaus)
  real(rp),    intent(out) :: gpgde(ndime,pgaus),gpadv(ndime,pgaus)
  real(rp),    intent(out) :: gplap(pnode,pgaus)
  integer(ip)              :: idime,jdime,inode,igaus,itime
  real(rp)                 :: fact1,fact2,fact3,gpve2
  !
  ! GPCOD: Coordinates
  !
  do igaus=1,pgaus
     do idime=1,ndime
        gpcod(idime,igaus)=0.0_rp
     end do
     do inode=1,pnode
        do idime=1,ndime
           gpcod(idime,igaus)=gpcod(idime,igaus)&
                +elcod(idime,inode)*gpsha(inode,igaus)
        end do
     end do
  end do
  !
  ! GPPRE: Pressure
  !
  if( kfl_regim_nsi == 1 .or. kfl_penal_nsi >= 1 ) then
     do igaus=1,pgaus
        gppre(igaus,1)=0.0_rp 
        if(kfl_penal_nsi==3) then                         ! Previous pressure : p(n,i-1) 
           gppre(igaus,2)=0.0_rp 
           do inode=1,pnode                  
              gppre(igaus,2)=gppre(igaus,2)&
                   +gpsha(inode,igaus)*elpre(inode,2)
           end do
        end if                                            ! Previous pressure : p(n-1) 
        do inode=1,pnode                                  ! Artifitial compressibility
           gppre(igaus,1)=gppre(igaus,1)&
                +gpsha(inode,igaus)*elpre(inode,1)
        end do
     end do
  end if
  !
  ! GPTEM: Temperature
  !
  !if(  kfl_regim_nsi==1.or.kfl_regim_nsi==2.or.kfl_regim_nsi==3.or.kfl_cotem_nsi==1.or.&
  !     lawvi_nsi(pmate)==1.or.lawde_nsi(pmate)==2) then
  if(  kfl_regim_nsi==1.or.kfl_regim_nsi==2.or.kfl_regim_nsi==3.or.kfl_cotem_nsi==1) then
     do igaus=1,pgaus
        gptem(igaus,1)=0.0_rp 
        do inode=1,pnode                  
           gptem(igaus,1)=gptem(igaus,1)&
                +gpsha(inode,igaus)*eltem(inode,1)
        end do
     end do
     !
     ! GPTEM+GPSGT: Temperature + temperature subgrid scale
     !
     if(kfl_sgste_nsi==1) then
        do igaus=1,pgaus
           gptem(igaus,1)=gptem(igaus,1)+gpsgt(igaus)
        end do
     end if
  end if
  if((kfl_regim_nsi==1.or.kfl_regim_nsi==3).and.kfl_timei_nsi==1) then
     do igaus=1,pgaus
        gptem(igaus,2)=0.0_rp 
        do inode=1,pnode                  
           gptem(igaus,2)=gptem(igaus,2)&
                +gpsha(inode,igaus)*eltem(inode,2)
        end do
     end do     
  end if
  !
  ! GPGRT: Temperature gradient
  !
  if(kfl_regim_nsi==1.or.kfl_regim_nsi==2.or.kfl_regim_nsi==3) then
     do igaus=1,pgaus
        do idime=1,ndime
           gpgrt(idime,igaus)=0.0_rp
        end do
        do inode=1,pnode
           do idime=1,ndime
              gpgrt(idime,igaus)=gpgrt(idime,igaus)&
                   +eltem(inode,1)*gpcar(idime,inode,igaus)
           end do
        end do
     end do     
  end if
  !
  ! GPVEL: Velocity
  !
  if(kfl_timei_nsi==1) then
     do itime=2,nbdfp_nsi
        do igaus=1,pgaus
           do idime=1,ndime
              gpvel(idime,igaus,itime)=0.0_rp
           end do
           do inode=1,pnode
              do idime=1,ndime
                 gpvel(idime,igaus,itime)=gpvel(idime,igaus,itime)&
                      +elvel(idime,inode,itime)*gpsha(inode,igaus)
              end do
           end do
        end do
     end do
  end if
  if(kfl_advec_nsi/=0) then
     do igaus=1,pgaus
        do idime=1,ndime
           gpvel(idime,igaus,1)=0.0_rp
        end do
        do inode=1,pnode
           do idime=1,ndime
              gpvel(idime,igaus,1)=gpvel(idime,igaus,1)&
                   +elvel(idime,inode,1)*gpsha(inode,igaus)
           end do
        end do
     end do
  end if
  !
  ! GPADV: Advection velocity
  !
  if(kfl_advec_nsi/=0) then
     do igaus=1,pgaus
        do idime=1,ndime
           gpadv(idime,igaus)=gpvel(idime,igaus,1)
        end do
     end do
  end if
  !
  ! Convection tracking: uc = u + u'
  !
  if( kfl_sgsco_nsi >= 1 ) then
     do igaus = 1,pgaus
        do idime = 1,ndime
           gpadv(idime,igaus) = gpadv(idime,igaus) + gpsgs(idime,igaus)
        end do
     end do
  end if
  !
  ! GPDIV: div(u)
  !
  if(kfl_regim_nsi==1.or.kfl_regim_nsi==2) then
     do igaus=1,pgaus
        gpdiv(igaus)=0.0_rp
        do inode=1,pnode
           do idime=1,ndime
              gpdiv(igaus)=gpdiv(igaus)&
                   +elvel(idime,inode,1)*gpcar(idime,inode,igaus)
           end do
        end do
     end do
  end if
  !
  ! GPPRE: Pressure for compressible regime and transient
  !
  if(kfl_timei_nsi==1.and.kfl_regim_nsi==1) then
     do igaus=1,pgaus
        gppre(igaus,2)=0.0_rp
        do inode=1,pnode
           gppre(igaus,2)=gppre(igaus,2)&
                +elpre(inode,2)*gpsha(inode,igaus)
        end do
     end do
     if(kfl_tisch_nsi==2) then
        do itime=3,1+kfl_tiacc_nsi
           do igaus=1,pgaus
              gppre(igaus,itime)=0.0_rp
              do inode=1,pnode
                 gppre(igaus,itime)=gppre(igaus,itime)&
                      +elpre(inode,itime)*gpsha(inode,igaus)
              end do
           end do
        end do
     end if
  end if
  !
  ! GPDEN: Pressure for compressible regime and transient
  !
  if(kfl_timei_nsi==1.and.kfl_regim_nsi==2) then 
     do igaus=1,pgaus
        gpden(igaus,2)=0.0_rp
        do inode=1,pnode
           gpden(igaus,2)=gpden(igaus,2)&
                +elden(inode,2)*gpsha(inode,igaus)
        end do
     end do
     if(kfl_tisch_nsi==2) then
        do itime=3,1+kfl_tiacc_nsi
           do igaus=1,pgaus
              gpden(igaus,itime)=0.0_rp
              do inode=1,pnode
                 gpden(igaus,itime)=gpden(igaus,itime)&
                      +elpre(inode,itime)*gpsha(inode,igaus)
              end do
           end do
        end do
     end if
  end if
  !
  ! GPGVE: Velocity gradients
  !
  if( (kfl_advec_nsi/=0.and.kfl_linea_nsi==2) .or. (kfl_cotur_nsi==-1.or.kfl_cotur_nsi==-2.or.kfl_cotur_nsi==-3.or.kfl_cotur_nsi==-10)  ) then
     do igaus=1,pgaus
        do idime=1,ndime
           do jdime=1,ndime
              gpgve(jdime,idime,igaus)=0.0_rp
           end do
        end do
        do inode=1,pnode
           do idime=1,ndime
              do jdime=1,ndime
                 gpgve(jdime,idime,igaus)=gpgve(jdime,idime,igaus)&
                      +elvel(idime,inode,1)*gpcar(jdime,inode,igaus)
              end do
           end do
        end do
     end do
  end if
  !
  ! GPGPR: Pressure gradients
  !
  if( kfl_regim_nsi==1 .or. kfl_sgsco_nsi /= 0 .or. kfl_sgsti_nsi /= 0 ) then
     do igaus=1,pgaus
        do idime=1,ndime
           gpgpr(idime,igaus,1)=0.0_rp  
           gpgpr(idime,igaus,2)=0.0_rp  
        end do
     end do
     do inode=1,pnode
        do igaus=1,pgaus
           do idime=1,ndime
              gpgpr(idime,igaus,1)=gpgpr(idime,igaus,1)&
                   +gpcar(idime,inode,igaus)*elpre(inode,1)
           end do
        end do
     end do
  end if
  !
  ! GPGRK: Turbulence kinetic energy gradients
  !
  if(kfl_grtur_nsi/=0) then
     do igaus=1,pgaus
        do idime=1,ndime
           gpgrk(idime,igaus)=0.0_rp  
        end do
     end do
     do inode=1,pnode
        do igaus=1,pgaus
           do idime=1,ndime
              gpgrk(idime,igaus)=gpgrk(idime,igaus)&
                   +gpcar(idime,inode,igaus)*elkin(inode)
           end do
        end do
     end do
  end if
  !
  ! GPGDE: Density gradients in compressible regime
  !
!!$  if( kfl_regim_nsi == 222 ) then
!!$
!!$     if(lawde_nsi(pmate)==1) then
!!$        !
!!$        ! Constant enthalpy:
!!$        !              gam      1                    u.grad(u)
!!$        ! grad(rho) = ----- --------- ( grad(p) + p* --------- )
!!$        !             gam-1 H-0.5*u^2                H-0.5*u^2
!!$        !
!!$        fact1=densi_nsi(1,pmate)/(densi_nsi(1,pmate)-1.0_rp)
!!$        do igaus=1,pgaus
!!$           gpve2=0.0_rp
!!$           do idime=1,ndime
!!$              gpve2=gpve2+gpvel(idime,igaus,1)*gpvel(idime,igaus,1)
!!$           end do
!!$           gpve2=1.0_rp/(densi_nsi(2,pmate)-0.5_rp*gpve2)
!!$           fact2=fact1*gpve2
!!$           do idime=1,ndime
!!$              fact3=0.0_rp
!!$              do jdime=1,ndime
!!$                 fact3=fact3+gpgve(jdime,idime,igaus)*gpvel(jdime,igaus,1)
!!$              end do
!!$              gpgde(idime,igaus)=fact2*(gpgpr(idime,igaus,1)&
!!$                   +gppre(igaus,1)*fact3*gpve2)
!!$           end do
!!$        end do
!!$
!!$     else if(lawde_nsi(pmate)==2) then
!!$        !
!!$        ! Compute via the solution of energy equation
!!$        !
!!$        do igaus=1,pgaus
!!$           do idime=1,ndime
!!$              gpgde(idime,igaus)=0.0_rp
!!$           end do
!!$           do inode=1,pnode
!!$              do idime=1,ndime
!!$                 gpgde(idime,igaus)=gpgde(idime,igaus)&
!!$                      +gpcar(idime,inode,igaus)*elden(inode,1)
!!$              end do
!!$           end do
!!$        end do
!!$
!!$     end if
!!$  end if
  !
  ! GPLAP: Laplacian
  !
  if(kfl_visco_nsi==1.and.plapl==1) then
     do igaus=1,pgaus
        do inode=1,pnode
           gplap(inode,igaus)=0.0_rp
           do idime=1,ndime
              gplap(inode,igaus)=gplap(inode,igaus)+gphes(idime,inode,igaus)
           end do
        end do
     end do
  end if
  !
  ! GPFLE: Level set function 
  !
  if( kfl_surte_nsi /= 0 ) then
     do igaus=1,pgaus
        gpfle(igaus)=0.0_rp
        do inode=1,pnode
           gpfle(igaus)=gpfle(igaus)+gpsha(inode,igaus)*elfle(inode)
        end do
     end do
  end if
  
end subroutine nsi_elmpre
