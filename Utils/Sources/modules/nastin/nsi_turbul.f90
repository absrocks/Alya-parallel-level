subroutine nsi_turbul(&
     itask,jtask,pnode,pgaus,igaui,igauf,kfl_cotur,&
     gpsha,gpcar,elvel,gpden,gpvis,gpmut, &
     gpgvi,grvis,gpgve,ielem,kfl_kemod)
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_turbul
  ! NAME 
  !    nsi_turbul
  ! DESCRIPTION
  !    Compute viscosity and its gradient due to turbulence
  !    GBMUT ....... mut
  !    JTASK = 1 ... Compute gradient GRVIS = grad(mut)
  !            0 ... Do not compute gradient GRVIS = grad(mut)
  !    ITASK = 1 ... GPVIS <= GPVIS + GPMUT
  !            0 ... Do not change GPVIS
  !    GPGVI = Viscosity gradient grad(mu+mut)
  ! USES
  !
  ! USED BY
  !    nsi_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp,nvart  
  use def_domain, only       :  mnode,ndime, lnods, walld
  use def_master, only       :  zeror, untur, postp
  use def_nastin, only       :  turbu_nsi,zensi,turmu_nsi,kfl_grvir_nsi
  use def_kermod, only       :  cmu_st, turmu_ker, kfl_logva
  use def_master, only       :  TUR_K_EPS_STD
  use mod_ker_regularization, only : regul_k, regul_e, kfl_regularization
  implicit none 
  integer(ip),  intent(in)   :: itask,jtask,pnode,pgaus,igaui,igauf
  integer(ip),  intent(in)   :: kfl_cotur, kfl_kemod, ielem
  real(rp),     intent(in)   :: gpsha(pnode,pgaus)
  real(rp),     intent(in)   :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)   :: elvel(ndime,pnode)
  real(rp),     intent(in)   :: gpden(pgaus)
  real(rp),     intent(out)  :: gpvis(pgaus)
  real(rp),     intent(inout):: gpmut(pgaus)
  real(rp),     intent(inout):: gpgvi(ndime,pgaus)
  real(rp),     intent(out)  :: grvis(ndime,pgaus)
  real(rp),     intent(out)  :: gpgve(ndime,ndime,pgaus)
  integer(ip)                :: igaus,inode,idime,jdime,kdime, ipoin
  real(rp)                   :: uaste,eltur(2, mnode),seci4
  real(rp)                   :: A0, As, W, phi, cmu, divve, simgr(3,3), kinen 
  real(rp)                   :: epsil, gpwal, omega, gpvor, a1, F2, Cr, f0, f02, Wsqr6
  real(rp)                   :: regue, reguk, sigmr
  if( kfl_cotur /= 0 ) then 

     if (turmu_ker % kfl_exist==0) then   ! Temporal, it should go to kernel
        !
        ! RANS models open integration
        !
        if( kfl_cotur == 1 ) then
           ! 
           ! K Epsilon model (open intergration) 
           !
           if (TUR_K_EPS_STD) then   ! open intergration points for K Epsilon model 
              do igaus = igaui,igauf          
                 kinen = 0.0_rp
                 epsil = 0.0_rp
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    kinen = kinen + gpsha(inode, igaus)*untur(1,ipoin,1)
                    epsil = epsil + gpsha(inode, igaus)*untur(2,ipoin,1)             
                 end do
                 gpmut(igaus) = max (cmu_st*gpden(igaus)*kinen*kinen/epsil, gpvis(igaus))
                 if (kfl_regularization) then
                    reguk = regul_k(kinen)
                    regue = regul_e(epsil)
                    gpmut(igaus)  = cmu_st*gpden(igaus)*reguk*reguk/regue
                 end if
                 if (kfl_logva==1)  gpmut(igaus)  = cmu_st*gpden(igaus)*exp(2.0_rp*kinen-epsil)
              end do
           end if


           !
           ! k eps REALIZABLE or kefp  model
           !        
           if ((kfl_kemod==2.or.kfl_kemod==3).and.itask==1) then
              A0 = (1.0_rp - 3.0_rp*sqrt(0.5_rp*cmu_st))/cmu_st
             
              do igaus = igaui,igauf
                 kinen = 0.0_rp
                 epsil = 0.0_rp
                 do inode = 1,pnode
                    ipoin = lnods(inode,ielem)
                    kinen = kinen + gpsha(inode, igaus)*untur(1,ipoin,1)
                    epsil = epsil + gpsha(inode, igaus)*untur(2,ipoin,1)             
                 end do
                 do idime = 1,ndime
                    do jdime = 1,ndime
                       gpgve(jdime,idime,igaus) = 0.0_rp
                    end do
                 end do
                 do inode = 1,pnode
                    do idime = 1,ndime
                       do jdime = 1,ndime
                          gpgve(jdime,idime,igaus) = gpgve(jdime,idime,igaus) &
                               + gpcar(jdime,inode,igaus) * elvel(idime,inode)
                       end do
                    end do
                 end do
                 divve =0.0_rp
                 do idime = 1,ndime
                    divve = divve + gpgve(idime, idime, igaus)
                    do jdime = 1,ndime
                       simgr(idime,jdime) = 0.5_rp*(gpgve(idime, jdime,igaus) &
                            +  gpgve(jdime, idime,igaus))
                    end do
                 end do
                 seci4 = 0.0_rp
                 W=0.0_rp
                 uaste =0.0_rp
                 do idime = 1,ndime
                    simgr(idime, idime ) = simgr(idime, idime) - divve/3.0_rp
                 end do
                 do idime = 1,ndime
                    do jdime = 1,ndime
                       seci4 = seci4 + simgr(idime, jdime) &        ! S_ij : S_ij
                            *simgr(idime, jdime)
                       uaste = uaste +gpgve(idime,jdime,igaus) &    ! D_ij : D_ij
                            *(gpgve(idime,jdime,igaus))   

                       do kdime =1, ndime
                          W = W +  simgr(idime,jdime)* &
                               simgr(jdime,kdime)* &
                               simgr(kdime,idime)
                       end do
                    end do
                 end do
                 uaste = sqrt(uaste -divve*divve /3.0_rp)
                 if (seci4.gt.zensi) W = W/sqrt(seci4*seci4*seci4)

                 Wsqr6 = sqrt(6.0_rp)*W     

                 Wsqr6 = min(Wsqr6,  1.0_rp)
                 Wsqr6 = max(Wsqr6, -1.0_rp)

                 phi = 1.0_rp/3.0_rp*acos(Wsqr6)
                 As= sqrt(6.0_rp)*cos(phi)

                 ! calculates turb viscosity
                 if (kfl_kemod==2) then ! REALIZABLE MODEL
                    cmu =1.0_rp/(A0+As*kinen/epsil*uaste)
                 else if (kfl_kemod==3) then ! KEFP
                    Cr  = 4.5_rp
                    f0  = Cr/(Cr-1.0_rp)
                    f02 = 4.0_rp*f0*(f0-1.0_rp)
                    if (kfl_regularization) then
                       reguk  = regul_k(kinen)
                       regue  = regul_e(epsil)
                    else
                       reguk  = kinen
                       regue  = epsil
                    end if
                    sigmr = cmu_st *(uaste*reguk/regue)*(uaste*reguk/regue)
                    cmu   = cmu_st *2.0_rp*f0 /(1.0_rp+ sqrt(1.0_rp + f02*sigmr))
                 end if

                 gpmut(igaus)= max(gpden(igaus)*cmu*kinen*kinen/epsil, gpvis(igaus))

                 if (kfl_regularization) &
                      gpmut(igaus)  = cmu*gpden(igaus)*reguk*reguk/regue
                
!!$              ! vis gradient set to zero
!!$              do idime = 1,ndime
!!$                 grvis(idime,igaus) = 0.0_rp
!!$              end do
              end do
           end if
           !
           ! Gradient of turbulent viscosity for RANS models
           !
           if (kfl_grvir_nsi == 0) then
              grvis = 0.0_rp
           else
              call runend('NSI_turbul: Gradient of turbulent viscosity in RANS not coded in nsi_turbul') 
           endif

        end if ! End RANS models

     end if ! not kerpro
     !
     ! Compute mu_t = rho * nu_t for LES
     !
     if( kfl_cotur < 0 ) gpmut(igaui:igauf) = gpden(igaui:igauf) * gpmut(igaui:igauf)

     !
     ! Compute effective viscosity & gradient of viscosity: 
     !    mu_eff        = mu + mu_t 
     !    grad (mu_eff) = grad(mu) + grad(mu_t)
     !
     if (itask == 1) then
        do igaus = igaui,igauf
           gpvis(igaus) = gpvis(igaus) + gpmut(igaus)
           do idime = 1,ndime
              gpgvi(idime,igaus) = gpgvi(idime,igaus) + grvis(idime,igaus)
           end do
        end do
     end if

     !----------------------------------------------------------------------
     !
     ! Needed for postprocessing turbulent viscosity
     !
     !----------------------------------------------------------------------
     if ( jtask == 1 ) then
        if( postp(1) % npp_stepi(74) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,74)) > zensi &
             .or. postp(1) % npp_stepi(75) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,75)) > zensi &
             .or. postp(1) % npp_stepi(76) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,76)) > zensi &
             .or. postp(1) % npp_stepi(77) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,77)) > zensi &
             .or. postp(1) % npp_stepi(78) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,78)) > zensi &
             .or. postp(1) % npp_stepi(92) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,92)) > zensi &
             .or. postp(1) % npp_stepi(93) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,93)) > zensi &
             .or. postp(1) % npp_stepi(94) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,94)) > zensi &
             .or. postp(1) % npp_stepi(95) /= 0 .or. maxval(postp(1) % pos_times(1:nvart,95)) > zensi ) then
           do igaus = igaui,igauf
              turmu_nsi(ielem)%a(1,igaus,1) = gpmut(igaus)
           end do
        end if
     end if

  end if

end subroutine nsi_turbul
