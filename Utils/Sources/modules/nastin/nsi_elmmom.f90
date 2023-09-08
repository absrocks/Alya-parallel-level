subroutine nsi_elmmom(&
     pnode,pevat,pgaus,plapl,ndofn,elpre,gpvis,gpgvi,&
     gptem,gpgrt,gprhs,gpsp2,gpden,rmomu,rcont,gpsha,gpcar,&
     gphes,gplap,p1vec,wmatr,wrhsi,wmatc)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_elmmom
  ! NAME 
  !    nsi_elmmom
  ! DESCRIPTION
  !    Compute the mementum equation at Gauss point
  !
  !    Three forms of the viscous term are considered.
  !    1. Laplacian form:  mu*grad(u):grad(v)
  !    2. Divergence form: 2*mu*eps(u):grad(v)=2*mu*eps(u):eps(v)
  !    3. Complete form:   2*mu*eps'(u):grad(v)
  !    eps(u) =1/2[grad(u)+grad(u)^t]
  !    eps'(u)=1/2[grad(u)+grad(u)^t]-1/3*div(u)I
  !    
  !    <-------- Laplacian
  !    <---------------------------- Divergence
  !    <-------------------------------------------------------- Complete
  !    mu*(dui/dxj*dvi/dxj)+mu*(dui/dxj*dvi/dxi)-2/3*mu*(dvi/dxi*dui/dxi)
  !
  ! USES
  ! USED BY
  !    nsi_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode,ntens
  use def_nastin, only     :  kfl_visco_nsi,kfl_grvis_nsi,kfl_savco_nsi,&
       &                      fvins_nsi,kfl_intpr_nsi,&
       &                      kfl_regim_nsi,grnor_nsi,gravi_nsi
  implicit none
  integer(ip), intent(in)  :: pnode,pevat,pgaus,plapl,ndofn
  real(rp),    intent(in)  :: elpre(pnode,2),gprhs(ndime+1,pgaus)
  real(rp),    intent(in)  :: gpvis(pgaus)
  real(rp),    intent(in)  :: gpgvi(ndime,pgaus),gptem(pgaus)
  real(rp),    intent(in)  :: gpgrt(ndime,pgaus),gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),    intent(in)  :: gphes(ntens,mnode,pgaus)
  real(rp),    intent(in)  :: gplap(pnode,pgaus)
  real(rp),    intent(in)  :: p1vec(ndime,ndime,pnode,pgaus)
  real(rp),    intent(in)  :: rmomu(ndime,ndime,pnode,pgaus)
  real(rp),    intent(in)  :: rcont(ndime+1,pnode,pgaus)
  real(rp),    intent(in)  :: gpsp2(pgaus),gpden(pgaus)
  real(rp),    intent(out) :: wmatr(pevat,pevat,pgaus)
  real(rp),    intent(out) :: wmatc(pevat,pevat,pgaus)
  real(rp),    intent(out) :: wrhsi(pevat,pgaus)
  integer(ip)              :: idofv,jdofv,jdofp,igaus,ievat,jevat
  integer(ip)              :: idime,inode,jnode,jdime,kdime
  real(rp)                 :: fact1,fact2

  if(kfl_savco_nsi/=2) then
     !
     ! Viscous term: ( mu*grad(u) , grad(v) ) or ( 2*mu*eps(u) , eps(v) )
     !
     if(kfl_visco_nsi==1) then
        !
        ! ( mu*dui/dxk , dv/dxk ) 
        !
        if(ndime==2) then
           do igaus=1,pgaus
              do inode=1,pnode
                 idofv=(inode-1)*ndofn
                 do idime=1,ndime
                    idofv=idofv+1
                    do jnode=inode+1,pnode
                       jdofv=(jnode-1)*ndofn+idime
                       fact1=0.0_rp
                       fact1=fact1+gpvis(igaus)*gpcar(1,inode,igaus)&
                            *gpcar(1,jnode,igaus)
                       fact1=fact1+gpvis(igaus)*gpcar(2,inode,igaus)&
                            *gpcar(2,jnode,igaus)
                       wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)+fact1
                       wmatr(jdofv,idofv,igaus)=wmatr(jdofv,idofv,igaus)+fact1
                    end do
                    wmatr(idofv,idofv,igaus)=wmatr(idofv,idofv,igaus)&
                         +gpvis(igaus)*gpcar(1,inode,igaus)&
                         *gpcar(1,inode,igaus)  
                    wmatr(idofv,idofv,igaus)=wmatr(idofv,idofv,igaus)&
                         +gpvis(igaus)*gpcar(2,inode,igaus)&
                         *gpcar(2,inode,igaus)  
                 end do
              end do
           end do
        else
           do igaus=1,pgaus
              do inode=1,pnode
                 idofv=(inode-1)*ndofn
                 do idime=1,ndime
                    idofv=idofv+1
                    do jnode=inode+1,pnode
                       jdofv=(jnode-1)*ndofn+idime
                       fact1=0.0_rp
                       fact1=fact1+gpvis(igaus)*gpcar(1,inode,igaus)&
                            *gpcar(1,jnode,igaus)
                       fact1=fact1+gpvis(igaus)*gpcar(2,inode,igaus)&
                            *gpcar(2,jnode,igaus)
                       fact1=fact1+gpvis(igaus)*gpcar(3,inode,igaus)&
                            *gpcar(3,jnode,igaus)
                       wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)+fact1
                       wmatr(jdofv,idofv,igaus)=wmatr(jdofv,idofv,igaus)+fact1
                    end do
                    wmatr(idofv,idofv,igaus)=wmatr(idofv,idofv,igaus)&
                         +gpvis(igaus)*gpcar(1,inode,igaus)&
                         *gpcar(1,inode,igaus)  
                    wmatr(idofv,idofv,igaus)=wmatr(idofv,idofv,igaus)&
                         +gpvis(igaus)*gpcar(2,inode,igaus)&
                         *gpcar(2,inode,igaus)  
                    wmatr(idofv,idofv,igaus)=wmatr(idofv,idofv,igaus)&
                         +gpvis(igaus)*gpcar(3,inode,igaus)&
                         *gpcar(3,inode,igaus)  
                 end do
              end do
           end do
        end if
     end if
     if(kfl_visco_nsi==1.and.fvins_nsi>0.9_rp) then
        !
        ! ( mu*duj/dxi , dv/dxj ) (only div form)
        !
        if(ndime==2) then
           do igaus=1,pgaus
              do inode=1,pnode
                 idofv=(inode-1)*ndofn
                 do idime=1,ndime
                    idofv=idofv+1
                    do jnode=1,pnode
                       jdofv=(jnode-1)*ndofn
                       fact1=gpvis(igaus)*gpcar(idime,jnode,igaus)     
                       jdofv=jdofv+1
                       wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                            +fact1*gpcar(1,inode,igaus)
                       jdofv=jdofv+1
                       wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                            +fact1*gpcar(2,inode,igaus)
                    end do
                    if(fvins_nsi==2.0_rp) then
                       fact1=-2.0_rp/3.0_rp*gpvis(igaus)*gpcar(idime,inode,igaus)
                       do jnode=1,pnode
                          jdofv=(jnode-1)*ndofn   
                          jdofv=jdofv+1
                          wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                               +fact1*gpcar(1,jnode,igaus)
                          jdofv=jdofv+1
                          wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                               +fact1*gpcar(2,jnode,igaus)
                       end do
                    end if
                 end do
              end do
           end do
        else
           do igaus=1,pgaus
              do inode=1,pnode
                 idofv=(inode-1)*ndofn
                 do idime=1,ndime
                    idofv=idofv+1
                    do jnode=1,pnode
                       jdofv=(jnode-1)*ndofn
                       fact1=gpvis(igaus)*gpcar(idime,jnode,igaus)     
                       jdofv=jdofv+1
                       wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                            +fact1*gpcar(1,inode,igaus)
                       jdofv=jdofv+1
                       wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                            +fact1*gpcar(2,inode,igaus)
                       jdofv=jdofv+1
                       wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                            +fact1*gpcar(3,inode,igaus)
                    end do
                    if(fvins_nsi==2.0_rp) then
                       fact1=-2.0_rp/3.0_rp*gpvis(igaus)*gpcar(idime,inode,igaus)
                       do jnode=1,pnode
                          jdofv=(jnode-1)*ndofn   
                          jdofv=jdofv+1
                          wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                               +fact1*gpcar(1,jnode,igaus)
                          jdofv=jdofv+1
                          wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                               +fact1*gpcar(2,jnode,igaus)
                          jdofv=jdofv+1
                          wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                               +fact1*gpcar(3,jnode,igaus)
                       end do
                    end if
                 end do
              end do
           end do
        end if
     end if
     !
     ! Pressure: -( p , div(v) )
     !  
     if(kfl_intpr_nsi==1) then
        !
        ! Contribution in LHS
        !
        if(ndime==2) then
           do igaus=1,pgaus
              fact1=1.0_rp
              do jnode=1,pnode
                 jdofp=jnode*ndofn
                 fact2=fact1*gpsha(jnode,igaus)
                 do inode=1,pnode
                    idofv=(inode-1)*ndofn
                    idofv=idofv+1
                    wmatr(idofv,jdofp,igaus)=wmatr(idofv,jdofp,igaus)&
                         -fact2*gpcar(1,inode,igaus)
                    idofv=idofv+1
                    wmatr(idofv,jdofp,igaus)=wmatr(idofv,jdofp,igaus)&
                         -fact2*gpcar(2,inode,igaus)
                 end do
              end do
           end do
        else
           do igaus=1,pgaus
              fact1=1.0_rp
              do jnode=1,pnode
                 jdofp=jnode*ndofn
                 fact2=fact1*gpsha(jnode,igaus)
                 do inode=1,pnode
                    idofv=(inode-1)*ndofn
                    idofv=idofv+1
                    wmatr(idofv,jdofp,igaus)=wmatr(idofv,jdofp,igaus)&
                         -fact2*gpcar(1,inode,igaus)
                    idofv=idofv+1
                    wmatr(idofv,jdofp,igaus)=wmatr(idofv,jdofp,igaus)&
                         -fact2*gpcar(2,inode,igaus)
                    idofv=idofv+1
                    wmatr(idofv,jdofp,igaus)=wmatr(idofv,jdofp,igaus)&
                         -fact2*gpcar(3,inode,igaus)
                 end do
              end do
           end do
        end if
     end if
     !
     ! Viscous term: ( -2*mu*div[eps(u)] , -v )
     !
     if(kfl_visco_nsi==1.and.plapl==1) then
        do igaus=1,pgaus
           do inode=1,pnode 
              idofv=(inode-1)*ndofn
              fact1=gpsha(inode,igaus)*gpvis(igaus)
              do idime=1,ndime
                 idofv=idofv+1
                 do jnode=1,pnode
                    !
                    ! ( mu*d^2ui/dxk^2, vj )
                    !
                    jdofv=(jnode-1)*ndofn+idime
                    wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                         +gplap(jnode,igaus)*fact1
                    !
                    ! ( mu*duk/(dxk dxi), vi ) (only div form)
                    !  It is not computed if we assume that duk/dxk=0, i.e. div u=0
                    !
                    if(fvins_nsi>0.9_rp) then
                       jdofv=(jnode-1)*ndofn
                       if(ndime==2) then
                          if(idime==1) then
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(1,jnode,igaus)
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(3,jnode,igaus)
                          else if(idime==2) then
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(3,jnode,igaus)
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(2,jnode,igaus) 
                          end if
                       else
                          if(idime==1) then
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(1,jnode,igaus)
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(4,jnode,igaus)
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(5,jnode,igaus)
                          else if(idime==2) then
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(4,jnode,igaus)
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(2,jnode,igaus) 
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(6,jnode,igaus) 
                          else if(idime==3) then
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(5,jnode,igaus)
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(6,jnode,igaus) 
                             jdofv=jdofv+1 
                             wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                                  +fact1*gphes(3,jnode,igaus) 
                          end if
                       end if
                    end if
                 end do
              end do
           end do
        end do
     end if
     !
     ! Viscous term: ( -2*grad(mu).eps(u) , -v )
     !
     if(kfl_visco_nsi==1.and.kfl_grvis_nsi==1) then
        do igaus=1,pgaus
           do inode=1,pnode
              idofv=(inode-1)*ndofn
              do idime=1,ndime
                 idofv=idofv+1
                 do jnode=1,pnode
                    jdofv=(jnode-1)*ndofn+idime
                    do kdime=1,ndime                                        ! dmu/dxk*dui/dxk*v
                       wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                            +gpgvi(kdime,igaus)*gpcar(kdime,jnode,igaus)&
                            *gpsha(inode,igaus)
                    end do
                    jdofv=(jnode-1)*ndofn
                    do jdime=1,ndime                                        ! dmu/dxk*duk/dxi*v
                       jdofv=jdofv+1
                       wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                            +gpgvi(jdime,igaus)*gpcar(idime,jnode,igaus)&
                            *gpsha(inode,igaus)
                    end do
                 end do
              end do
           end do
        end do
     end if
  end if
  if(kfl_savco_nsi==1) then
     do igaus=1,pgaus
        do ievat=1,pevat
           do jevat=1,pevat
              wmatc(jevat,ievat,igaus)=wmatr(jevat,ievat,igaus)
              wmatr(jevat,ievat,igaus)=0.0_rp
           end do
        end do
     end do
  end if
  !
  ! Velocity: div(u)*tau2'*div(v)
  !
  if(ndime==2) then
     do igaus=1,pgaus
        do inode=1,pnode
           idofv=(inode-1)*ndofn
           do idime=1,ndime
              idofv=idofv+1
              fact1=gpsp2(igaus)*gpcar(idime,inode,igaus)
              do jnode=1,pnode              
                 jdofv=(jnode-1)*ndofn   
                 jdofv=jdofv+1
                 fact2=fact1*rcont(1,jnode,igaus)
                 wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)+fact2
                 jdofv=jdofv+1
                 fact2=fact1*rcont(2,jnode,igaus)
                 wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)+fact2
              end do
           end do
        end do
     end do
  else
     do igaus=1,pgaus
        do inode=1,pnode
           idofv=(inode-1)*ndofn
           do idime=1,ndime
              idofv=idofv+1
              fact1=gpsp2(igaus)*gpcar(idime,inode,igaus)
              do jnode=1,pnode              
                 jdofv=(jnode-1)*ndofn   
                 jdofv=jdofv+1
                 fact2=fact1*rcont(1,jnode,igaus)
                 wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)+fact2
                 jdofv=jdofv+1
                 fact2=fact1*rcont(2,jnode,igaus)
                 wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)+fact2
                 jdofv=jdofv+1
                 fact2=fact1*rcont(3,jnode,igaus)
                 wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)+fact2
              end do
           end do
        end do
     end do
  end if
  !
  ! Velocity: rho*(rho/T*DT/Dt-rho/p*dp/dt)*tau2'*div(v)
  !
  if(kfl_regim_nsi>=1) then
     do igaus=1,pgaus
        do inode=1,pnode
           idofv=(inode-1)*ndofn
           do idime=1,ndime
              idofv=idofv+1
              wrhsi(idofv,igaus)=wrhsi(idofv,igaus)&
                   +gpsp2(igaus)*gpcar(idime,inode,igaus)&
                   *gprhs(ndime+1,igaus)
           end do
        end do
     end do
  end if
  !
  ! Velocity: Other terms ( rmomu(u) , p1vec(v) )
  !
  if(ndime==2) then
     do igaus=1,pgaus
        do inode=1,pnode
           idofv=(inode-1)*ndofn
           do idime=1,ndime
              idofv=idofv+1
              do jnode=1,pnode
                 jdofv=(jnode-1)*ndofn
                 do jdime=1,ndime
                    jdofv=jdofv+1
                    wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                         +p1vec(idime,    1,inode,igaus)&
                         *rmomu(    1,jdime,jnode,igaus)
                    wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                         +p1vec(idime,    2,inode,igaus)&
                         *rmomu(    2,jdime,jnode,igaus)
                 end do
              end do
           end do
        end do
     end do
  else
     do igaus=1,pgaus
        do inode=1,pnode
           idofv=(inode-1)*ndofn
           do idime=1,ndime
              idofv=idofv+1
              do jnode=1,pnode
                 jdofv=(jnode-1)*ndofn
                 do jdime=1,ndime
                    jdofv=jdofv+1
                    wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                         +p1vec(idime,    1,inode,igaus)&
                         *rmomu(    1,jdime,jnode,igaus)
                    wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                         +p1vec(idime,    2,inode,igaus)&
                         *rmomu(    2,jdime,jnode,igaus)
                    wmatr(idofv,jdofv,igaus)=wmatr(idofv,jdofv,igaus)&
                         +p1vec(idime,    3,inode,igaus)&
                         *rmomu(    3,jdime,jnode,igaus)
                 end do
              end do
           end do
        end do
     end do
  end if
  !
  ! Pressure: ( grad(p) , p1vec(v)-v )
  !
  if(kfl_intpr_nsi==1) then
     fact1=1.0_rp
  else
     fact1=0.0_rp
  end if

  do igaus=1,pgaus
     do inode=1,pnode
        idofv=(inode-1)*ndofn
        do idime=1,ndime
           idofv=idofv+1
           do jnode=1,pnode
              jdofp=jnode*ndofn
              wmatr(idofv,jdofp,igaus)=wmatr(idofv,jdofp,igaus)&
                   -fact1*gpsha(inode,igaus)*gpcar(idime,jnode,igaus)
              do jdime=1,ndime
                 wmatr(idofv,jdofp,igaus)=wmatr(idofv,jdofp,igaus)&
                      +p1vec(idime,jdime,inode,igaus)*gpcar(jdime,jnode,igaus)
              end do
           end do
        end do
     end do
  end do
  !
  ! External force: (f,v)
  !
  do igaus=1,pgaus
     do inode=1,pnode
        idofv=(inode-1)*ndofn
        do idime=1,ndime
           idofv=idofv+1
           do kdime=1,ndime
              wrhsi(idofv,igaus)=wrhsi(idofv,igaus)&
                   +p1vec(idime,kdime,inode,igaus)*gprhs(kdime,igaus)
           end do
        end do
     end do
  end do

end subroutine nsi_elmmom
