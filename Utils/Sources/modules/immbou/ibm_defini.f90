subroutine ibm_defini(itask)
  !-----------------------------------------------------------------------
  !****f* master/ibm_defini
  ! NAME
  !    ibm_defini
  ! DESCRIPTION
  !    This routines creates the IB structure
  ! USED BY
  !    Turnon
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_immbou
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iimbo,ii,jj,inoib,ipoib,ipoin,iboib,idime,jdime
  integer(ip)             :: ielem,inode,kk
  real(rp),    pointer    :: F(:,:), Fv(:), Fp(:)
  real(rp),    pointer    :: T(:,:), Tv(:), Tp(:)

  if( itask == 0 ) then

     !-------------------------------------------------------------------
     !
     ! Initialization
     !
     !-------------------------------------------------------------------

     call ibm_memall(2_ip)
     do iimbo = 1,nimbo         

        imbou(iimbo) % kfl_insid  =  0
        imbou(iimbo) % kfl_typeb  =  0
        imbou(iimbo) % kfl_coupl  =  0
        imbou(iimbo) % kfl_model  =  0

        imbou(iimbo) % bobox      =  0.0_rp
        imbou(iimbo) % massa      = -1.0_rp
        imbou(iimbo) % densi      = -1.0_rp 
        imbou(iimbo) % volum      =  0.0_rp
        imbou(iimbo) % momin      = -1.0_rp
        imbou(iimbo) % posgr      = -1.0e12_rp

        imbou(iimbo) % npoib      =  0
        imbou(iimbo) % nboib      =  0
        imbou(iimbo) % npoin      =  0
        imbou(iimbo) % nelem      =  0

        imbou(iimbo) % posil      =  -1.0e12_rp
        imbou(iimbo) % velol      =  0.0_rp
        imbou(iimbo) % accel      =  0.0_rp
        imbou(iimbo) % force      =  0.0_rp
        imbou(iimbo) % vforce     =  0.0_rp
        imbou(iimbo) % pforce     =  0.0_rp

        imbou(iimbo) % torqu      =  0.0_rp
        imbou(iimbo) % vtorqu     =  0.0_rp
        imbou(iimbo) % ptorqu     =  0.0_rp

        imbou(iimbo) % posia      =  0.0_rp
        imbou(iimbo) % veloa      =  0.0_rp
        imbou(iimbo) % accea      =  0.0_rp
        imbou(iimbo) % rotac      =  0.0_rp
        imbou(iimbo) % quate(1,1) =  1.0_rp
        imbou(iimbo) % quate(1,2) =  1.0_rp

        do idime = 1,3
           imbou(iimbo) % quate(idime+1,1) = 0.0_rp
           imbou(iimbo) % quate(idime+1,2) = 0.0_rp                      
           do jdime = 1,3
              imbou(iimbo) % rotac(idime,jdime) = 0.0_rp
           end do
        end do

        do idime = 1,3
           imbou(iimbo) % rotac(idime,idime) = 1.0_rp      
        end do

     end do
     ncoll_ibm =  0

  else if( itask == 1 ) then

     !-------------------------------------------------------------------
     !
     ! Used for Parall service
     !
     !-------------------------------------------------------------------
     !
     ! Read in ibm_reaphy
     !
     !
     ! Read in ibm_reaphy
     !
     call iexcha( kfl_timei_ibm )
     call iexcha( kfl_rotib_ibm )
     call iexcha( kfl_linib_ibm )
     do ii = 1,3
        call iexcha( nstro_ibm(ii) )
        call iexcha( nstli_ibm(ii) )
     end do
     call iexcha( kfl_mvext_ibm )
     call iexcha( kfl_ralei_ibm )
     call iexcha( nstra_ibm )
     call iexcha( nenra_ibm )
     call iexcha( kfl_staib_ibm )
     call iexcha( kfl_colli_ibm )
     call iexcha( kfl_coibm )
     call iexcha( kfl_cofor )
     call iexcha( kfl_grafo_ibm )             ! Gravity force
     call iexcha( kfl_catfo_ibm )             ! Catamaran force
     call iexcha( kfl_buofo_ibm )             ! Buoyancy force
     call iexcha( kfl_drafo_ibm )             ! Drag force
     call iexcha( kfl_extfo_ibm )             ! External force 

     call rexcha( spher_ibm )                 ! Sphericity 
     call rexcha( staib_ibm )

     call rexcha( ralei_ibm )
     !
     ! Read in ibm_reanut
     !
     call iexcha( kfl_inter_ibm )
     call iexcha( kfl_nforc_ibm )
     call iexcha( kfl_advec )

     call rexcha( safet_ibm )
     call rexcha( beta_ibm  )
     call rexcha( gamma_ibm )
     !
     ! Read in ibm_reageo
     !
     do iimbo = 1,nimbo
        call iexcha( imbou(iimbo) % npoib )
        call iexcha( imbou(iimbo) % nboib )
        call iexcha( imbou(iimbo) % npoin )
        call iexcha( imbou(iimbo) % nelem )
        call iexcha( imbou(iimbo) % kfl_insid )
        call iexcha( imbou(iimbo) % kfl_typeb )
        call iexcha( imbou(iimbo) % kfl_coupl )
        do ii = 1,3
           do jj = 1,2
              call rexcha( imbou(iimbo) % bobox(ii,jj) )
           end do
        end do
        call rexcha( imbou(iimbo) % massa )
        call rexcha( imbou(iimbo) % densi )
        call rexcha( imbou(iimbo) % volum )
        do ii = 1,6
           call rexcha( imbou(iimbo) % momin(ii) )
        end do
        do ii = 1,3
           call rexcha( imbou(iimbo) % posgr(ii) ) 
        end do
        !
        ! Linear motion
        !
        do ii = 1,3
           do jj = 1,2
              call rexcha( imbou(iimbo) % posil(ii,jj) )
           end do
        end do
        do ii = 1,3
           do jj = 1,2
              call rexcha( imbou(iimbo) % velol(ii,jj) )
           end do
        end do
        do ii = 1,3
           do jj = 1,3
              call rexcha( imbou(iimbo) % accel(ii,jj) )
           end do
        end do
        do ii = 1,3
           do jj = 1,2
              call rexcha( imbou(iimbo) % force(ii,jj) )
           end do
        end do
        do ii = 1,3
           call rexcha( imbou(iimbo) % vforce(ii) )
        end do
        do ii = 1,3
           call rexcha( imbou(iimbo) % pforce(ii) )
        end do
        !
        ! Angular motion
        !
        do ii = 1,3
           do jj = 1,2
              call rexcha( imbou(iimbo) % posia(ii,jj) )
           end do
        end do
        do ii = 1,3
           do jj = 1,2
              call rexcha( imbou(iimbo) % veloa(ii,jj) )
           end do
        end do
        do ii = 1,3
           do jj = 1,3
              call rexcha( imbou(iimbo) % accea(ii,jj) )
           end do
        end do
        do jj = 1,3
           do ii = 1,3
              call rexcha( imbou(iimbo) % rotac(ii,jj) )
           end do
        end do
        do ii = 1,2
           do jj = 1,4
              call rexcha( imbou(iimbo) % quate(jj,ii) )
           end do
        end do
        do ii = 1,3
           do jj = 1,2
              call rexcha( imbou(iimbo) % torqu(ii,jj) )
           end do
        end do
        do ii = 1,3
           call rexcha( imbou(iimbo) % vtorqu(ii) )
        end do
        do ii = 1,3
           call rexcha( imbou(iimbo) % ptorqu(ii) )
        end do
     end do

  else if( itask == 5 ) then

     !-------------------------------------------------------------------
     !
     ! Point to coordinates and connectivity
     !
     !-------------------------------------------------------------------

     !do iimbo = 1,nimbo
     !   imbou(iimbo) % cooib     => cooib
     !   imbou(iimbo) % lnoib     => lnoib
     !end do

  else if( itask == 6 ) then

     !-------------------------------------------------------------------
     !
     ! Reduce sum on particles pressure+viscous forces and torques
     !
     !-------------------------------------------------------------------

     if( IPARALL ) then
        call memgen(zero,nimbo*12,zero)
        ii = 0
        do iimbo = 1,nimbo
           Fv        => imbou(iimbo)%vforce
           Fp        => imbou(iimbo)%pforce
           Tv        => imbou(iimbo)%vtorqu
           Tp        => imbou(iimbo)%ptorqu
           ii        =  ii + 1
           gesca(ii) =  Fv(1)
           ii        =  ii + 1
           gesca(ii) =  Fv(2)
           ii        =  ii + 1
           gesca(ii) =  Fv(3)
           ii        =  ii + 1
           gesca(ii) =  Fp(1)
           ii        =  ii + 1
           gesca(ii) =  Fp(2)
           ii        =  ii + 1
           gesca(ii) =  Fp(3)
           ii        =  ii + 1
           gesca(ii) =  Tv(1)
           ii        =  ii + 1
           gesca(ii) =  Tv(2)
           ii        =  ii + 1
           gesca(ii) =  Tv(3)
           ii        =  ii + 1
           gesca(ii) =  Tp(1)
           ii        =  ii + 1
           gesca(ii) =  Tp(2)
           ii        =  ii + 1
           gesca(ii) =  Tp(3)
        end do

        nparr =  nimbo*12
        parre => gesca
        call Parall(9_ip)

        ii = 0
        do iimbo = 1,nimbo
           Fv    => imbou(iimbo)%vforce
           Fp    => imbou(iimbo)%pforce
           Tv    => imbou(iimbo)%vtorqu
           Tp    => imbou(iimbo)%ptorqu
           ii    =  ii + 1
           Fv(1) =  gesca(ii)
           ii    =  ii + 1
           Fv(2) =  gesca(ii) 
           ii    =  ii + 1
           Fv(3) =  gesca(ii) 
           ii    =  ii + 1
           Fp(1) =  gesca(ii)
           ii    =  ii + 1
           Fp(2) =  gesca(ii) 
           ii    =  ii + 1
           Fp(3) =  gesca(ii) 
           ii    =  ii + 1
           Tv(1) =  gesca(ii) 
           ii    =  ii + 1
           Tv(2) =  gesca(ii) 
           ii    =  ii + 1
           Tv(3) =  gesca(ii) 
           ii    =  ii + 1
           Tp(1) =  gesca(ii) 
           ii    =  ii + 1
           Tp(2) =  gesca(ii) 
           ii    =  ii + 1
           Tp(3) =  gesca(ii) 
        end do
        call memgen(two,nimbo*12,zero)
     end if

  else if( itask == 7 ) then

     !-------------------------------------------------------------------
     !
     ! Used for Parall service
     !
     !-------------------------------------------------------------------

     do iimbo = 1,nimbo
        !
        ! Boundary mesh
        !
        do ipoib = 1,imbou(iimbo) % npoib
           do idime = 1,ndime
              call rexcha( imbou(iimbo)%cooin(idime,ipoib) ) ! COOIN
           end do
        end do
        do ipoib = 1,imbou(iimbo) % npoib
           do idime = 1,ndime
              call rexcha( imbou(iimbo)%cooib(idime,ipoib) ) ! COOIB
           end do
        end do
        do ipoib = 1,imbou(iimbo) % npoib
           do idime = 1,ndime
              call rexcha( imbou(iimbo)%cooi2(idime,ipoib) ) ! COOI2
           end do
        end do
        do iboib = 1,imbou(iimbo) % nboib
           do inoib = 1,mnoib
              call iexcha( imbou(iimbo)%lnoib(inoib,iboib) ) ! LNOIB
           end do
        end do
        do iboib = 1,imbou(iimbo) % nboib
           call iexcha( imbou(iimbo)%ltyib(iboib) )          ! LTYIB
        end do
        !
        ! Volume mesh
        !
        do ipoin = 1,imbou(iimbo) % npoin
           do idime = 1,ndime
              call rexcha( imbou(iimbo)%coord(idime,ipoin) ) ! COORD
           end do
        end do
        do ielem = 1,imbou(iimbo) % nelem
           do inode = 1,mnodi
              call iexcha( imbou(iimbo)%lnods(inode,ielem) ) ! LNODS
           end do
        end do
        do ielem = 1,imbou(iimbo) % nelem
           call iexcha( imbou(iimbo)%ltype(ielem) )          ! LTYPE
        end do
        

     end do

  else if( itask == 8 ) then

     !-------------------------------------------------------------------
     !
     ! Reduce sum on particles pressure+viscous forces and torques
     !
     !-------------------------------------------------------------------

     if( IPARALL ) then
        call memgen(zero,nimbo*6,zero)
        ii = 0
        kk = 0
        do iimbo = 1,nimbo
           if( imbou(iimbo) % kfl_typeb > 0 ) then
              kk        =  kk + 1
              F         => imbou(iimbo) % force
              T         => imbou(iimbo) % torqu
              ii        =  ii + 1
              gesca(ii) =  F(1,1)
              ii        =  ii + 1
              gesca(ii) =  F(2,1)
              ii        =  ii + 1
              gesca(ii) =  F(3,1)
              ii        =  ii + 1
              gesca(ii) =  T(1,1)
              ii        =  ii + 1
              gesca(ii) =  T(2,1)
              ii        =  ii + 1
              gesca(ii) =  T(3,1)
           end if
        end do
        
        if( kk > 0 ) then
           nparr =  nimbo*6
           parre => gesca
           call Parall(9_ip)
           ii = 0
           do iimbo = 1,nimbo
              if( imbou(iimbo) % kfl_typeb > 0 ) then
                 F      => imbou(iimbo) % force
                 T      => imbou(iimbo) % torqu
                 ii     =  ii + 1
                 F(1,1) =  gesca(ii)
                 ii     =  ii + 1
                 F(2,1) =  gesca(ii) 
                 ii     =  ii + 1
                 F(3,1) =  gesca(ii) 
                 ii     =  ii + 1
                 T(1,1) =  gesca(ii) 
                 ii     =  ii + 1
                 T(2,1) =  gesca(ii) 
                 ii     =  ii + 1
                 T(3,1) =  gesca(ii)
              end if
           end do
        end if

        call memgen(two,nimbo*6,zero)
     end if

  else if( itask == 9 ) then

     !-------------------------------------------------------------------
     !
     ! Reduce sum on particles viscous forces
     !
     !-------------------------------------------------------------------

     if( IPARALL ) then
        call memgen(zero,nimbo*3,zero)
        ii = 0
        do iimbo = 1,nimbo
           Fv        => imbou(iimbo)%vforce
           ii        =  ii + 1
           gesca(ii) =  Fv(1)
           ii        =  ii + 1
           gesca(ii) =  Fv(2)
           ii        =  ii + 1
           gesca(ii) =  Fv(3)
        end do

        nparr =  nimbo*3
        parre => gesca
        call Parall(9_ip)

        ii = 0
        do iimbo = 1,nimbo
           Fv    => imbou(iimbo)%vforce
           ii    =  ii + 1
           Fv(1) =  gesca(ii)
           ii    =  ii + 1
           Fv(2) =  gesca(ii) 
           ii    =  ii + 1
           Fv(3) =  gesca(ii) 

        end do
        call memgen(two,nimbo*3,zero)
     end if

  end if

end subroutine ibm_defini
