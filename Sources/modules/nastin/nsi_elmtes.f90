subroutine nsi_elmtes(&
     itask,pnode,pgaus,plapl,gpsha,gpcar,gphes,gpden,&
     gppor,gpvis,gpgvi,gpsp1,gpsp3,gptt1,gptt2,gpadv,&
     gplap,gpdiv,p1vec,p2vec,p2sca) 
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmtes
  ! NAME 
  !    nsi_elmtes
  ! DESCRIPTION
  !    Compute the test function at the Gauss point
  ! USES
  ! USED BY
  !    nsi_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  mnode,ndime,ntens
  use def_nastin, only     :  kfl_advec_nsi,kfl_visco_nsi,kfl_grvis_nsi,&
       &                      corio_nsi,fvela_nsi,fvins_nsi,kfl_regim_nsi
  implicit none
  integer(ip), intent(in)  :: itask,pnode,pgaus,plapl
  real(rp),    intent(in)  :: gpden(pgaus),gpvis(pgaus)
  real(rp),    intent(in)  :: gppor(pgaus),gpgvi(ndime,pgaus)
  real(rp),    intent(in)  :: gpadv(ndime,pgaus),gplap(pnode,pgaus)
  real(rp),    intent(in)  :: gpdiv(pgaus)
  real(rp),    intent(in)  :: gpsp1(pgaus),gpsp3(pgaus)
  real(rp),    intent(in)  :: gptt1(pgaus),gptt2(pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(out) :: p1vec(ndime,ndime,pnode,pgaus)
  real(rp),    intent(out) :: p2vec(ndime,pnode,pgaus),p2sca(pnode,pgaus)
  integer(ip)              :: idime,inode,jdime,igaus
  real(rp)                 :: fact1,fact2,factx,facty,factz

  if(itask==1.or.itask==2) then
     !
     ! P1VEC: normal assembly needed for momentum equation 
     !
     do igaus=1,pgaus                
        do inode=1,pnode
           do idime=1,ndime
              do jdime=1,ndime
                 p1vec(jdime,idime,inode,igaus)=0.0_rp
              end do
           end do
        end do
     end do
     !
     ! p1vec= (tau1'/tau1)*v - tau1'*( -rho*(uc.grad)v - 2*rho*(w x v) + sig*v )
     !      = v*[ (tau1'/tau1) - tau1'*sig ] + tau1'*( rho*(uc.grad)v + 2*rho*(w x v) )
     !
     do igaus=1,pgaus
        do inode=1,pnode
           fact1=gpsha(inode,igaus)&                      ! v*[(tau1'/tau1)-tau1'*sig]
                *(gptt1(igaus)-gpsp1(igaus)*gppor(igaus)) 
           do idime=1,ndime           
              p1vec(idime,idime,inode,igaus)=fact1
           end do
        end do
     end do
     if(kfl_advec_nsi==1) then ! OJO
        !
        ! tau1'*rho*(uc.grad)v
        !
        do igaus=1,pgaus
           fact1=gpsp1(igaus)*gpden(igaus)
           do inode=1,pnode
              fact2=0.0_rp
              do idime=1,ndime
                 fact2=fact2+gpadv(idime,igaus)*gpcar(idime,inode,igaus)
              end do
              fact2=fact1*fact2
              do idime=1,ndime
                 p1vec(idime,idime,inode,igaus)=p1vec(idime,idime,inode,igaus)+fact2
              end do
           end do
        end do
     end if
     if(corio_nsi>1.0e-10_rp) then
        !
        ! tau1'*2*rho*(w x v)
        ! x-equation: v=(v,0,0) => w x v = (    0 , wz*v , -wy*v)     
        ! y-equation: v=(0,v,0) => w x v = (-wz*v ,    0 ,  wx*v)     
        ! z-equation: v=(0,0,v) => w x v = ( wy*v ,-wx*v ,     0)     
        !
        if(ndime==2) then
           do igaus=1,pgaus           
              factz=gpsp1(igaus)*2.0_rp*gpden(igaus)*fvela_nsi(3)
              do inode=1,pnode 
                 p1vec(1,2,inode,igaus)=p1vec(1,2,inode,igaus)+factz*gpsha(inode,igaus)
                 p1vec(2,1,inode,igaus)=p1vec(2,1,inode,igaus)-factz*gpsha(inode,igaus)
              end do
           end do
        else
           do igaus=1,pgaus           
              factx=gpsp1(igaus)*2.0_rp*gpden(igaus)*fvela_nsi(1)
              facty=gpsp1(igaus)*2.0_rp*gpden(igaus)*fvela_nsi(2)
              factz=gpsp1(igaus)*2.0_rp*gpden(igaus)*fvela_nsi(3)
              do inode=1,pnode
                 p1vec(1,2,inode,igaus)=p1vec(1,2,inode,igaus)+factz*gpsha(inode,igaus)
                 p1vec(1,3,inode,igaus)=p1vec(1,3,inode,igaus)-facty*gpsha(inode,igaus)
                 p1vec(2,1,inode,igaus)=p1vec(2,1,inode,igaus)-factz*gpsha(inode,igaus)
                 p1vec(2,3,inode,igaus)=p1vec(2,3,inode,igaus)+factx*gpsha(inode,igaus)
                 p1vec(3,1,inode,igaus)=p1vec(3,1,inode,igaus)+facty*gpsha(inode,igaus)
                 p1vec(3,2,inode,igaus)=p1vec(3,2,inode,igaus)-factx*gpsha(inode,igaus)  
              end do
           end do

        end if
     end if
     !
     ! Viscous term: 2*tau1'*mu*div[eps(v)] 
     !
     if(kfl_visco_nsi==1.and.plapl==1) then
        do igaus=1,pgaus
           fact1=gpsp1(igaus)*gpvis(igaus)
           do inode=1,pnode
              !
              ! tau1'*mu*( d^2v/dxk^2 )
              !
              fact2=gplap(inode,igaus)*fact1
              do idime=1,ndime
                 p1vec(idime,idime,inode,igaus)=p1vec(idime,idime,inode,igaus)+fact2
              end do
              !
              ! tau1'*mu*( d^2v/dxi^2 ) not computed
              !
              if(fvins_nsi>10.9_rp) then
                 do idime=1,ndime
                    p1vec(idime,idime,inode,igaus)=p1vec(idime,idime,inode,igaus)&
                         +fact1*gphes(idime,inode,igaus)
                 end do
              end if
           end do
        end do
     end if
     !
     ! Viscous term: tau1'*grad(mu).eps(u): tau'*(grad(mu).grad(v)+dmu/dxj*dvj/dxi)
     !
     if(kfl_visco_nsi==1.and.kfl_grvis_nsi==1) then
        if(ndime==2) then
           do igaus=1,pgaus
              do inode=1,pnode
                 fact1=0.0_rp
                 do jdime=1,ndime
                    fact1=fact1+gpgvi(jdime,igaus)*gpcar(jdime,inode,igaus)
                 end do
                 fact1=fact1*gpsp1(igaus)
                 p1vec(1,1,inode,igaus)=p1vec(1,1,inode,igaus)&
                      +fact1+gpsp1(igaus)*gpgvi(1,igaus)*gpcar(1,inode,igaus) 
                 p1vec(2,2,inode,igaus)=p1vec(2,2,inode,igaus)&
                      +fact1+gpsp1(igaus)*gpgvi(2,igaus)*gpcar(2,inode,igaus) 
              end do
           end do
        else
           do igaus=1,pgaus
              do inode=1,pnode
                 fact1=0.0_rp
                 do jdime=1,ndime
                    fact1=fact1+gpgvi(jdime,igaus)*gpcar(jdime,inode,igaus)
                 end do
                 fact1=fact1*gpsp1(igaus)
                 p1vec(1,1,inode,igaus)=p1vec(1,1,inode,igaus)&
                      +fact1+gpsp1(igaus)*gpgvi(1,igaus)*gpcar(1,inode,igaus) 
                 p1vec(2,2,inode,igaus)=p1vec(2,2,inode,igaus)&
                      +fact1+gpsp1(igaus)*gpgvi(2,igaus)*gpcar(2,inode,igaus) 
                 p1vec(3,3,inode,igaus)=p1vec(3,3,inode,igaus)&
                      +fact1+gpsp1(igaus)*gpgvi(3,igaus)*gpcar(3,inode,igaus) 
              end do
           end do
        end if
     end if
  end if

  if(itask==1.or.itask==3) then
     !
     ! P2VEC, P2SCA: needed for continuity equation 
     !
     !
     ! P2VEC: tau1'*rho*grad(q)
     !
     do igaus=1,pgaus
        !fact1=gpsp1(igaus)*gpden(igaus)
        fact1=gpsp1(igaus) ! OJO CAMBIAR
        do inode=1,pnode
           do idime=1,ndime
              p2vec(idime,inode,igaus)=fact1*gpcar(idime,inode,igaus)
           end do
        end do
     end do
     !
     ! P2SCA: (tau2^{-1}*tau2')*q
     !
     do igaus=1,pgaus
        do inode=1,pnode
           p2sca(inode,igaus)=gptt2(igaus)*gpsha(inode,igaus)
        end do
     end do
     if(kfl_regim_nsi==1.or.kfl_regim_nsi==2) then
        !
        ! P2SCA: P2SCA + tau3*(u.grad(q) - div(u)*q)
        !
        do igaus=1,pgaus
           do inode=1,pnode
              fact2=0.0_rp
              do idime=1,ndime
                 fact2=fact2+gpadv(idime,igaus)*gpcar(idime,inode,igaus)
              end do
              p2sca(inode,igaus)=p2sca(inode,igaus)&
                   +gpsp3(igaus)*(fact2-gpdiv(igaus)*gpsha(inode,igaus))
           end do
        end do
     end if
  end if

end subroutine nsi_elmtes
