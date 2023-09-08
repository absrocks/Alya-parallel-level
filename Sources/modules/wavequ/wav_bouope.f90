subroutine wav_bouope(itask)
  !-----------------------------------------------------------------------
  !****f* Wavequ/wav_bouope
  ! NAME 
  !    wav_bouope
  ! DESCRIPTION
  !    Compute the contribution of the boundary condition
  !
  !    g0*Mu^{n+1} + (A+g1*M)u^{n-1} + g2*Mu^{n-1} = f^n
  !
  !    - Open boundary condition: flux c^2*grad(u).n is assembled.
  !      ------------------------
  !                     +-
  !       LHS <= LHS -  |  c^2*v*grad(u^n).n ds
  !                    -+ S
  !
  !    - Absorbing boundary condition: grad(u).n + alpha*du/dt = 0
  !      -----------------------------
  !
  !      In the explicit case the time derivative du/dt en n is_ 
  !  
  !      Backward: du/dt=(u^{n}   - u^{n-1}}/dt
  !      Centered: du/dt=(u^{n+1} - u^{n-1}}/(2*dt)
  !      Forward:  du/dt=(u^{n+1} - u^{n}}/dt
  !
  !      - Option 1: du/dt is computed as du/dt=(u^{n}-u^{n-1}}/dt
  !        ---------      
  !                      +-
  !        LHS <= LHS +  |   alpha*c^2/dt*Ni*Nj ds
  !                    -+ S
  !                      +-
  !        RHS <= RHS +  |   alpha*c^2/dt*Ni*u^{n-1} ds
  !                     -+ S
  !
  !      - Option 2: du/dt is computed as du/dt=(u^{n+1}-u^{n-1}}/(2*dt)
  !        ---------
  !
  !        The following integrals are computed using a closed rule,
  !        so that the mass matrix is modified.
  !                      +-
  !        LHS <= LHS +  |   alpha*c^2/(2*dt)*Ni*Nj ds
  !                    -+ S
  !                      +-
  !        RHS <= RHS +  |   alpha*c^2/(2*dt)*Ni*u^{n-1} ds
  !                     -+ S
  !
  ! USES
  ! USED BY
  !    wav_solite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_wavequ
  implicit none
  integer(ip), intent(in) :: itask
  real(rp)                :: elmat(mnode,mnode)
  real(rp)                :: elcod(ndime,mnode)
  real(rp)                :: bowav(mnodb)
  real(rp)                :: bocod(ndime,mnodb)
  real(rp)                :: gbbas(ndime,ndime)
  real(rp)                :: gbcar(ndime,mnode)
  real(rp)                :: gpcar(ndime,mnode,mgaus)
  integer(ip)             :: ielem,inode,ipoin,jnode,idime,knode
  integer(ip)             :: igaus,igaub,iboun,inodb,pblty,jnodb
  integer(ip)             :: pnodb,pmate,pnode,pelty,pgaus,pgaub,itime
  real(rp)                :: gbden,gbkap,gpden,gpkap,gbdet,gbsur,gpwav
  real(rp)                :: fact0,fact1,fact2,gpdet,xjaci(9),xjacm(9)

  if(kfl_onbou_wav==1) then

     select case(itask)

     case(1)
        !
        ! Left-hand side
        !
        do iboun=1,nboun

           if(  kfl_fixbo_wav(iboun)==2.or.&
                kfl_fixbo_wav(iboun)==3.or.&
                kfl_fixbo_wav(iboun)==4.or.&
                kfl_fixbo_wav(iboun)==5) then
              !
              ! Element characteristics
              !
              pblty=ltypb(iboun)
              pnodb=nnode(pblty)
              pgaub=ngaus(pblty)
              ielem=lelbo(iboun)
              pelty=ltype(ielem)
              pnode=nnode(pelty)
              pmate=lmate_wav(ielem)
              pgaus=ngaus(pelty)
              !
              ! Initialization
              !
              do inode=1,pnode
                 do jnode=1,pnode
                    elmat(inode,jnode)=0.0_rp
                 end do
              end do
              !
              ! Gather operations BOCOD and ELCOD
              !
              do inodb=1,pnodb
                 ipoin=lnodb(inodb,iboun)
                 do idime=1,ndime
                    bocod(idime,inodb)=coord(idime,ipoin)
                 end do
              end do
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 do idime=1,ndime
                    elcod(idime,inode)=coord(idime,ipoin)
                 end do
              end do

              if(kfl_timet_wav==2) then
                 !
                 ! Implicit
                 !
                 !
                 ! Cartesian derivatives GPCAR
                 !
                 do igaus=1,pgaus
                    call elmder(&
                         pnode,ndime,elmar(pelty)%deriv(1,1,igaus),& 
                         elcod,gpcar(1,1,igaus),gpdet,xjacm,xjaci)
                 end do
                 ! 
                 ! Properties GPDEN and GPKAP
                 !
                 gpden=densi_wav(1,pmate)
                 gpkap=kappa_wav(1,pmate)
                 !
                 ! Loop over Gauss points
                 !
                 gauss_points: do igaub=1,pgaub
                    !
                    ! Jacobian
                    !
                    call bouder(&
                         pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&
                         bocod,gbbas,gbdet)
                    gbsur=elmar(pblty)%weigp(igaub)*gbdet 
                    !
                    ! Properties at IGAUB
                    !
                    gbden=0.0_rp
                    gbkap=0.0_rp
                    do igaus=1,pgaus
                       do inodb=1,pnodb                  
                          fact0=elmar(pelty)%shaga(igaus,lboel(inodb,iboun))&
                               *elmar(pblty)%shape(inodb,igaub)
                          gbden=gbden+gpden*fact0
                          gbkap=gbkap+gpkap*fact0
                       end do
                    end do

                    if(kfl_fixbo_wav(iboun)==2) then
                       !
                       ! Open boundary condition
                       !
                       !
                       ! Ensure normal GBBAS(:,NDIME) is exterior
                       !
                       call chenor(pnode,gbbas,bocod,elcod)    
                       !
                       ! Cartesian derivatives GBCAR at IGAUB 
                       !
                       do inode=1,pnode                                       
                          do idime=1,ndime                                    
                             gbcar(idime,inode)=0.0_rp
                             do inodb=1,pnodb
                                knode=lboel(inodb,iboun)
                                do igaus=1,pgaus 
                                   gbcar(idime,inode)=gbcar(idime,inode)&
                                        +elmar(pelty)%shaga(igaus,knode)&
                                        *elmar(pblty)%shape(inodb,igaub)&
                                        *gpcar(idime,inode,igaus)
                                end do
                             end do
                          end do
                       end do
                       !
                       ! Element matrix
                       !
                       fact0=gbsur*gbkap/gbden
                       do inodb=1,pnodb
                          fact1=fact0*elmar(pblty)%shape(inodb,igaub)
                          inode=lboel(inodb,iboun)
                          do jnode=1,pnode
                             fact2=0.0_rp
                             do idime=1,ndime
                                fact2=fact2+gbcar(idime,jnode)*gbbas(idime,ndime)
                             end do
                             elmat(inode,jnode)=elmat(inode,jnode)-fact1*fact2
                          end do
                       end do

                    else if(kfl_fixbo_wav(iboun)==3.or.kfl_fixbo_wav(iboun)==4) then
                       !
                       ! Absorbing boundary condition 
                       !
                       fact0=gbsur*gbkap/gbden*dtinv_wav
                       if(kfl_fixbo_wav(iboun)==4) fact0=0.5_rp*fact0
                       if(kfl_tisch_wav==3) fact0=fact0*nebet_wav
                       do inodb=1,pnodb
                          inode=lboel(inodb,iboun)
                          fact1=fact0*elmar(pblty)%shape(inodb,igaub)
                          do jnodb=1,pnodb
                             jnode=lboel(jnodb,iboun)
                             elmat(inode,jnode)=elmat(inode,jnode)&
                                  +fact1*elmar(pblty)%shape(jnodb,igaub)
                          end do
                       end do

                    end if

                 end do gauss_points
                 !
                 ! Boundary condition
                 !
                 if(kfl_onnod_wav==1.and.kfl_timet_wav==2) then
                    do inode=1,pnode
                       ipoin=lnods(inode,ielem)
                       if(kfl_fixno_wav(ipoin)==1) then
                          do jnode=1,pnode
                             elmat(inode,jnode)=0.0_rp
                          end do
                       end if
                    end do
                 end if
                 !
                 ! Assemble matrix
                 !
                 call assmat(&
                      solve(1)%ndofn,pnode,pnode,npoin,solve(1)%kfl_algso,&
                      ielem,lnods(1,ielem),elmat,amatr)

              else
                 !
                 ! Explicit 
                 !
                 if(kfl_fixbo_wav(iboun)==4.or.kfl_fixbo_wav(iboun)==5) then
                    !
                    ! Absorbing boundary condition
                    !
                    ! 
                    ! Properties GPDEN and GPKAP
                    !
                    gpden=densi_wav(1,pmate)
                    gpkap=kappa_wav(1,pmate)
                    !
                    ! Loop over Gauss points
                    !
                    do igaub=1,pnodb
                       !
                       ! Jacobian
                       !
                       call bouder(&
                            pnodb,ndime,ndimb,elmar(pblty)%deric(1,1,igaub),&
                            bocod,gbbas,gbdet)
                       gbsur=elmar(pblty)%weigc(igaub)*gbdet 
                       !
                       ! Properties at IGAUB
                       !
                       gbden=0.0_rp
                       gbkap=0.0_rp
                       do igaus=1,pgaus            
                          fact0=elmar(pelty)%shaga(igaus,lboel(igaub,iboun))
                          gbden=gbden+gpden*fact0
                          gbkap=gbkap+gpkap*fact0
                       end do
                       fact0=gbsur*gbkap/gbden*dtinv_wav
                       if(kfl_fixbo_wav(iboun)==4) fact0=0.5_rp*fact0
                       ipoin=lnodb(igaub,iboun)
                       vmass_wav(ipoin)=vmass_wav(ipoin)+fact0

                    end do

                 end if

              end if

           end if

        end do

     case(2)
        !
        ! Right-hand side
        !
        do iboun=1,nboun

           if(  kfl_fixbo_wav(iboun)==3.or.&
                kfl_fixbo_wav(iboun)==4.or.&
                kfl_fixbo_wav(iboun)==5) then

              pblty=ltypb(iboun)
              pnodb=nnode(pblty)
              pgaub=ngaus(pblty)
              ielem=lelbo(iboun)
              pelty=ltype(ielem)
              pnode=nnode(pelty)
              pmate=lmate_wav(ielem)
              pgaus=ngaus(pelty)
              !
              ! Gather BOCOD
              !
              do inodb=1,pnodb
                 ipoin=lnodb(inodb,iboun)
                 do idime=1,ndime
                    bocod(idime,inodb)=coord(idime,ipoin)
                 end do
              end do              
              !
              ! Gather operations BOWAV
              !
              if(kfl_fixbo_wav(iboun)==3) then
                 !
                 ! Backward: BOWAV=u^{n}-u^{n-1}
                 !
                 do inodb=1,pnodb
                    ipoin=lnodb(inodb,iboun)
                    bowav(inodb)=wavam(ipoin,3)-wavam(ipoin,4)
                 end do

              else if(kfl_fixbo_wav(iboun)==4) then
                 !
                 ! Centered: BOWAV=-u^{n-1}/2.0
                 !
                 do inodb=1,pnodb
                    ipoin=lnodb(inodb,iboun)
                    bowav(inodb)=-0.5_rp*wavam(ipoin,4)
                 end do

              else if(kfl_fixbo_wav(iboun)==5) then
                 !
                 ! Forward:  BOWAV=-u^{n}
                 !
                 do inodb=1,pnodb
                    ipoin=lnodb(inodb,iboun)
                    bowav(inodb)=-wavam(ipoin,3)   
                 end do
              end if
              ! 
              ! Properties GPDEN and GPKAP
              !
              gpden=densi_wav(1,pmate)
              gpkap=kappa_wav(1,pmate)

              if(kfl_timet_wav==2) then
                 !
                 ! Implicit absorbing (1st and 2nd order time derivative)
                 !
                 do igaub=1,pgaub
                    !
                    ! Jacobian GBDET
                    !
                    call bouder(&
                         pnodb,ndime,ndimb,elmar(pblty)%deriv(1,1,igaub),&
                         bocod,gbbas,gbdet)
                    gbsur=elmar(pblty)%weigp(igaub)*gbdet 
                    !
                    ! Properties GBDEN and GBKAP at IGAUB
                    !
                    gbden=0.0_rp
                    gbkap=0.0_rp
                    do igaus=1,pgaus
                       do inodb=1,pnodb                  
                          fact0=elmar(pelty)%shaga(igaus,lboel(inodb,iboun))&
                               *elmar(pblty)%shape(inodb,igaub)
                          gbden=gbden+gpden*fact0
                          gbkap=gbkap+gpkap*fact0
                       end do
                    end do
                    !
                    ! Absorbing boundary condition:
                    !               +-
                    ! RHS <= RHS +  |   c^2*alpha/dt*Ni*u^{n-1} ds
                    !              -+ S
                    !
                    gpwav=0.0_rp
                    do inodb=1,pnodb
                       gpwav=gpwav+elmar(pblty)%shape(inodb,igaub)*bowav(inodb)
                    end do
                    fact0=-gbsur*gbkap/gbden*dtinv_wav*gpwav
                    if(kfl_fixbo_wav(iboun)==4) fact0=0.5_rp*fact0
                    if(kfl_tisch_wav==3) fact0=fact0*nebet_wav
                    do inodb=1,pnodb
                       ipoin=lnodb(inodb,iboun)
                       rhsid(ipoin)=rhsid(ipoin)&
                            +fact0*elmar(pblty)%shape(inodb,igaub)
                    end do

                 end do

              else
                 !
                 ! Explicit absorbing (1st and 2nd order time derivative)
                 !
                 do igaub=1,pnodb ! Use a closed rule: Gauss points are nodes
                    !
                    ! Jacobian
                    !
                    call bouder(&
                         pnodb,ndime,ndimb,elmar(pblty)%deric(1,1,igaub),&
                         bocod,gbbas,gbdet)
                    gbsur=elmar(pblty)%weigc(igaub)*gbdet 
                    !
                    ! Properties GBDEN and GBKAP at IGAUB
                    !
                    gbden=0.0_rp
                    gbkap=0.0_rp
                    do igaus=1,pgaus                 
                       fact0=elmar(pelty)%shaga(igaus,lboel(igaub,iboun))
                       gbden=gbden+gpden*fact0
                       gbkap=gbkap+gpkap*fact0
                    end do
                    !
                    ! Absorbing boundary condition: 
                    !               +-
                    ! RHS <= RHS +  |   c^2*alpha/(dt)*Ni*u^{n} ds
                    !              -+ S
                    !               +-
                    ! RHS <= RHS +  |   c^2*alpha/(2*dt)*Ni*u^{n-1} ds
                    !              -+ S
                    !
                    fact0=-gbsur*gbkap/gbden*dtinv_wav*bowav(igaub)
                    ipoin=lnodb(igaub,iboun)
                    rhsid(ipoin)=rhsid(ipoin)+fact0

                 end do

              end if

           end if
           
        end do

     end select

  end if

end subroutine wav_bouope
