subroutine ibm_forces()
  !-----------------------------------------------------------------------
  !****f* ibm_forces/ibm_forces
  ! NAME
  !    ibm_forces
  ! DESCRIPTION
  !    This routines solves the Euler and Newton equations for rigid bodies
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_kermod
  use def_domain
  use def_immbou
  implicit none
  integer(ip)       :: iimbo,idime
  real(rp)          :: grafo,buofo,accel(3)
  real(rp), pointer :: F(:,:),Fp(:),Fv(:)     ! Linear motion
  real(rp), pointer :: T(:,:),Tp(:),Tv(:)     ! Angular motion

  integer(ip)       :: inodb,pnodb,igaub,pgaub,pblty,iboun,ipoin
  real(rp)          :: x(3),F1(3),T1(3),gbpos(3),gbsur,massa
  real(rp)          :: bocod(ndime,mnodb),eucta,surfa,dummr
  real(rp)          :: baloc(ndime,ndime),bouno(3),xx(3)
  real(rp)          :: auxii(3),auxim(3)

  !----------------------------------------------------------------------
  !
  ! Initialize
  ! 
  !----------------------------------------------------------------------

  do iimbo = 1,nimbo
     F  => imbou(iimbo) % force
     T  => imbou(iimbo) % torqu
     do idime = 1,3
        F(idime,1) = 0.0_rp
        T(idime,1) = 0.0_rp
     end do
  end do

  if( kfl_extfo_ibm /= 0 ) then
     if( INOTMASTER ) then
        do iimbo = 1,nimbo
           xx = 0.0_rp
           surfa = 0.0_rp
           boundaries: do iboun=1,imbou(iimbo)%nboib
              !
              ! Pointers and element properties
              !
              F     => imbou(iimbo)%force
              T     => imbou(iimbo)%torqu
              pblty =  imbou(iimbo)%ltyib(iboun)
              pnodb =  nnode(pblty)
              pgaub =  ngaib(pblty)
              !
              ! BOCOD: Gather operations
              !
              do inodb = 1,pnodb 
                 ipoin = imbou(iimbo)%lnoib(inodb,iboun)
                 do idime = 1,ndime
                    bocod(idime,inodb) = imbou(iimbo)%cooib(idime,ipoin)
                 end do
              end do
              !
              ! BOUNO: Exterior normal
              !
              call exteib(iimbo,pnodb,iboun,bouno,0_ip)
              do igaub = 1,pgaub
                 !
                 ! Jacobian
                 !
                 call bouder(&
                      pnodb,ndime,ndimb,elmar(pblty)%derib(1,1,igaub),& 
                      bocod,baloc,eucta) 
                 gbsur = elmar(pblty)%weiib(igaub)*eucta 
                 !
                 ! Gauss point coordinates
                 !
                 do idime = 1,ndime
                    gbpos(idime) = 0.0_rp
                    do inodb = 1,pnodb
                       gbpos(idime) = gbpos(idime) &
                            + elmar(pblty)%shaib(inodb,igaub) * bocod(idime,inodb)
                    end do
                 end do
                 do idime = 1,ndime
                    xx(idime) = xx(idime) + gbpos(idime) * gbsur
                 end do
                 surfa = surfa + gbsur
                 !accel(1) = 0.0_rp
                 !accel(2) = 0.0_rp
                 !accel(3) = 0.0_rp
                 !call extefo(& 
                 !     4_ip,gbpos,&
                 !     imbou(iimbo) % densi,spher_ibm,denme,visme,cutim,accel)
                 !do idime = 1,3
                 !   F1(idime) = accel(1) * bouno(idime)
                 !end do
                 !
                 ! T1: Compute torque= F1 x X
                 !
                 !do idime = 1,ndime
                 !   x(idime) = gbpos(idime) - imbou(iimbo)%posil(idime,1) 
                 !end do
                 !call vecasi(3_ip,x,T1)                     ! T1 = X-Xg
                 !call vecpro(T1,F1,T1,3_ip)                 ! T1 = (X-Xg) x F1
                 !
                 ! F, T: Actualize force and torque 
                 !
                 !do idime = 1,3
                 !   F(idime) = F(idime) + F1(idime) * gbsur
                 !   T(idime) = T(idime) + T1(idime) * gbsur
                 !end do
              end do
           end do boundaries
           do idime = 1,ndime
              xx(idime) = xx(idime) / surfa
           end do
           accel(1) = 0.0_rp
           accel(2) = 0.0_rp
           accel(3) = 0.0_rp
           call extefo(& 
                kfl_extfo_ibm,xx,&
                imbou(iimbo) % densi,spher_ibm,denme,visme,cutim,accel)         
           do idime = 1,3
              F(idime,1) = F(idime,1) + imbou(iimbo) % densi * imbou(iimbo) % volum * accel(idime)
           end do
        end do
     end if
  end if


!!$  if( kfl_extfo_ibm /= 0 ) then               
!!$     accel(1:3) = 0.0_rp
!!$     F1(1:3)    = 0.0_rp
!!$     x(1:3)     = 0.0_rp
!!$
!!$     do iimbo = 1,nimbo
!!$        if(      imbou(iimbo) % kfl_typeb == IBM_EMBEDDED .or. &
!!$             & ( imbou(iimbo) % kfl_typeb >  0 .and. INOTMASTER ) ) then
!!$
!!$           surfa =  0.0_rp
!!$           massa =  imbou(iimbo) % densi * imbou(iimbo) % volum 
!!$           F     => imbou(iimbo) % force
!!$           T     => imbou(iimbo) % torqu
!!$
!!$           boundaries: do iboun = 1,imbou(iimbo) % nboib
!!$              !
!!$              ! Element properties
!!$              !
!!$              pblty =  imbou(iimbo) % ltyib(iboun)
!!$              pnodb =  nnode(pblty)
!!$              pgaub =  ngaib(pblty)
!!$              !
!!$              ! BOCOD: Gather operations
!!$              !
!!$              do inodb = 1,pnodb 
!!$                 ipoin = imbou(iimbo) % lnoib(inodb,iboun)
!!$                 do idime = 1,ndime
!!$                    bocod(idime,inodb) = imbou(iimbo) % cooib(idime,ipoin)
!!$                 end do
!!$              end do
!!$              !
!!$              ! Loop of Gauss points
!!$              !
!!$              do igaub = 1,pgaub
!!$                 !
!!$                 ! Jacobian
!!$                 !
!!$                 call bouder(&
!!$                      pnodb,ndime,ndimb,elmar(pblty) % derib(1,1,igaub),& 
!!$                      bocod,baloc,eucta)
!!$                 call vecuni(ndime,baloc(1:ndime,ndime),dummr)
!!$                 
!!$                 gbsur = elmar(pblty) % weiib(igaub) * eucta 
!!$                 !
!!$                 ! Gauss point coordinates
!!$                 !
!!$                 do idime = 1,ndime
!!$                    gbpos(idime) = 0.0_rp
!!$                    do inodb = 1,pnodb
!!$                       gbpos(idime) = gbpos(idime) &
!!$                            + elmar(pblty) % shaib(inodb,igaub) * bocod(idime,inodb)
!!$                    end do
!!$                 end do
!!$                 call extefo(& 
!!$                      kfl_extfo_ibm,gbpos,&
!!$                      imbou(iimbo) % densi,spher_ibm,denme,visme,cutim,accel)
!!$                 do idime = 1,ndime
!!$                    F1(idime) =  accel(idime) * baloc(idime,ndime)
!!$                 end do
!!$                 !
!!$                 ! T1: Compute torque= F1 x X
!!$                 !
!!$                 do idime = 1,ndime
!!$                    x(idime) = gbpos(idime) - imbou(iimbo) % posil(idime,1) 
!!$                 end do
!!$                 call vecasi(3_ip,x,T1)                     ! T1 = X-Xg
!!$                 call vecpro(x,F1,T1,3_ip)                  ! T1 = (X-Xg) x F1  (pressure) 
!!$                 !
!!$                 ! F, T: Actualize force and torque 
!!$                 !
!!$                 do idime = 1,3
!!$                    F(idime,1) = F(idime,1) + F1(idime) * gbsur
!!$                    T(idime,1) = T(idime,1) + T1(idime) * gbsur
!!$                 end do
!!$
!!$              end do
!!$           end do boundaries
!!$        end if
!!$     end do
!!$     !
!!$     ! Sum forces between slaves ony for body fitted-like particles
!!$     ! 
!!$     call ibm_defini(8_ip)
!!$  end if

  !----------------------------------------------------------------------
  !
  ! Add force and torque due to coupling with other modules
  ! 
  !----------------------------------------------------------------------
  if( coupling('IMMBOU','NASTIN') >= 1 ) then

     do iimbo = 1,nimbo
        F   => imbou(iimbo) % force
        Fp  => imbou(iimbo) % pforce
        Fv  => imbou(iimbo) % vforce
        T   => imbou(iimbo) % torqu
        Tp  => imbou(iimbo) % ptorqu
        Tv  => imbou(iimbo) % vtorqu
        do idime = 1,3
           F(idime,1) = F(idime,1) + Fp(idime) + Fv(idime) 
           T(idime,1) = T(idime,1) + Tp(idime) + Tv(idime) 
        end do
     end do

  else if( kfl_cofor == 1 .or. kfl_cofor == 2 .or. kfl_cofor == 3 ) then

     do iimbo = 1,nimbo
        F   => imbou(iimbo) % force
        Fp  => imbou(iimbo) % pforce
        Fv  => imbou(iimbo) % vforce
        T   => imbou(iimbo) % torqu
        Tp  => imbou(iimbo) % ptorqu
        Tv  => imbou(iimbo) % vtorqu
        if( kfl_cofor == 1 ) then
           do idime = 1,3
              F(idime,1) = F(idime,1) + Fp(idime) + Fv(idime) 
              T(idime,1) = T(idime,1) + Tp(idime) + Tv(idime) 
           end do
        else if( kfl_cofor == 2 ) then
           do idime = 1,3
              F(idime,1) = F(idime,1) + Fv(idime) 
              T(idime,1) = T(idime,1) + Tv(idime) 
           end do
        else if( kfl_cofor == 3 ) then
           do idime = 1,3
              F(idime,1) = F(idime,1) + Fp(idime)
              T(idime,1) = T(idime,1) + Tp(idime)
           end do
        end if
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! Drag force
  ! 
  !----------------------------------------------------------------------

  if( kfl_drafo_ibm /= 0 ) then

     call ibm_forsph()
     do iimbo = 1,nimbo
        F   => imbou(iimbo) % force
        Fv  => imbou(iimbo) % vforce
        do idime = 1,ndime
           F(idime,1) = F(idime,1) + Fv(idime) 
        end do
     end do

  end if

  !----------------------------------------------------------------------
  !
  ! Gravity and Buoyancy
  ! 
  !----------------------------------------------------------------------

  grafo = real( kfl_grafo_ibm , rp )  ! Gravity  force = 1.0
  buofo = real( kfl_buofo_ibm , rp )  ! Buoyancy force = 1.0

  do iimbo = 1,nimbo
     F  => imbou(iimbo) % force
     do idime = 1,ndime
        F(idime,1) = F(idime,1) + grnor * gravi(idime) &
             &     * imbou(iimbo) % volum &
             &     * ( grafo * imbou(iimbo) % densi - buofo * denme )
     end do
  end do

  !do idime = 1,ndiem
  !F(idime) = Fp(idime) - grnor * gravi(idime) * imbou(iimbo) % volum * densi_nsi(1,1)

  !  print *, "grnor: ",grnor 
  !  print *, "gravi: ",gravi(:) 
  !  print *, "volun: ",imbou(1) % volum 
  !  print *, "grafo: ",grafo 
  !  print *, "densi: ",imbou(1) % densi
  !  print *, "buofo: ",buofo 
  !  print *, "denme: ",denme
  !  pause
  !  print*,F(2);stop
  !iimbo=1
  !print*,' '
  !print*,'b=',imbou(1) % force(2)

  !----------------------------------------------------------------------
  !
  ! Catamaran force
  ! 
  !----------------------------------------------------------------------

  if( kfl_catfo_ibm /= 0 ) then ! when I move it to ale put it more general and erase this
     
     F   => imbou(1) % force
     T   => imbou(1) % torqu
     F1(1) = 775.0_rp   ! without - because we have x and y axis oposite to CHE
     F1(2) = 1600.0_rp  ! without - because we have x and y axis oposite to CHE         
     F1(3) = 200.0_rp
     x(1) =  0.0_rp   ! x -X_cog
     x(2) =  0.0_rp
     x(3) =  7.0_rp

     call vecpro(x,F1,T1,3_ip)                  ! T1 = (X) x F1  (pressure) 

     do idime = 1,3
        F(idime,1) = F(idime,1) + F1(idime)
        T(idime,1) = T(idime,1) + T1(idime)
     end do
 
  end if

end subroutine ibm_forces
