subroutine tur_clippi()
  !------------------------------------------------------------------------
  !****f* Turbul/tur_clippi
  ! NAME 
  !    tur_clippi
  ! DESCRIPTION
  !    This routine performs the unknown cliiping just after calling 
  !    the solver
  ! USED BY
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_turbul !,  only     :  TUR_K_EPS_STD, TUR_FAMILY_K_OMEGA
  use def_kermod , only     :  kfl_logva, cmu_st
  use mod_ker_proper, only  :  ker_proper
  use mod_ker_regularization, only : kfl_regularization
  implicit none
  integer(ip) :: ipoin,itotn,dummi
  real(rp)    :: zezer(4),Cmu,L,nu,k,eps,phi,f,rho(1),mu(1),turbm(1), dummr
  real(rp)    :: keyam, epsam

  if (kfl_logva==1.or.kfl_regularization) then
     return
     if(kfl_clipp_tur==2) then ! upper limiter whe log variables 
        do ipoin=1,npoin
           unkno(ipoin)=min(13.0_rp ,unkno(ipoin)) 
           unkno(ipoin)=max(-19.0_rp ,unkno(ipoin))
        end do
     end if
     return
  end if
  if(kfl_clipp_tur==-2) then
     !
     ! No clipping
     !
     continue

  else if(kfl_clipp_tur==-1) then

     if( INOTMASTER ) then
        if(kfl_algor_tur==1) then
           do ipoin=1,npoin
              if(unkno(ipoin)==0.0_rp) unkno(ipoin)=zetur          
           end do
        else
           do ipoin=1,npoin
              do iunkn_tur=1,nturb_tur
                 itotn=(ipoin-1)*nturb_tur+iunkn_tur
                 if(kfl_fixno_tur(1,ipoin,iunkn_tur)<=0) &
                      unkno(itotn)=max(zezer(iunkn_tur),unkno(itotn))                 
              end do
           end do
        end if
     end if

  else if(kfl_clipp_tur<=2) then
     !
     ! Minimum prescribed
     !
     if(kfl_clipp_tur==0) then
        zezer = zetur            
        if (kfl_custo.ge.2) then ! ABL flows, ambient values           
            keyam = 1.0e-4_rp    ! tke ambient value
            epsam = 7.208e-8_rp  ! diss ambient value
            zezer(1) = keyam      
            zezer(2) = epsam
!            else if(TUR_FAMILY_K_OMEGA) then ! w = eps/(beta*tke)
!               zezer(2) = epsam/(cmu_st*keyam) ! BAD should be actual tke
!            end if
        end if
     else if(kfl_clipp_tur==1) then
        zezer=0.0_rp
     else if(kfl_clipp_tur==2) then
        !
        ! To be used when using wall-functions
        !
        if( ISLAVE ) call runend('TUR_CLIPPI: VODOM IS NEEDED IN PARALLEL')
        call ker_proper('DENSI','IPOIN',ipoin,dummi,rho(1:))
        call ker_proper('VISCO','IPOIN',ipoin,dummi,mu(1:))
        nu  = mu(1)/rho(1)
        Cmu = param_tur(6)
        L   = vodom**(1.0_rp/real(ndime,rp))
        k   = nu**2.0_rp* 1296.0_rp/L**2.0_rp*sqrt(Cmu) ! k
        eps = nu**3.0_rp*46656.0_rp/L**4.0_rp*Cmu       ! eps
        phi = zetur                                     ! phi
        f   = zetur                                     ! f
        zezer(1) = k
        zezer(2) = eps
        zezer(3) = phi
        zezer(4) = f
        if(TUR_FAMILY_K_OMEGA) zezer(2)=Cmu*eps/k
     else if(kfl_clipp_tur==3) then
        zezer=1.0e-15_rp
     end if
     if( INOTMASTER ) then
        if(kfl_algor_tur==1) then
           if (TUR_FAMILY_K_OMEGA.and.iunkn_tur==2) then
              ! limit diss= cmu*k*w > eps_atm
              do ipoin=1,npoin 
                 if ((kfl_fixno_tur(1,ipoin,iunkn_tur)<=0))&  !.or.(kfl_custo.ge.2).and.walld(ipoin).gt.10.0) & !.and.lpoty(ipoin).le.0 )) & 
                      unkno(ipoin)=max(zezer(iunkn_tur)/(cmu_st*untur(1,ipoin,1)),unkno(ipoin)) 
              end do
           else 
              do ipoin=1,npoin
                 if ((kfl_fixno_tur(1,ipoin,iunkn_tur)<=0))&  !.or.(kfl_custo.ge.2).and.walld(ipoin).gt.10.0) & !.and.lpoty(ipoin).le.0 )) & 
                      unkno(ipoin)=max(zezer(iunkn_tur),unkno(ipoin)) 
              end do
           end if
           if (kfl_custo.ge.2.and.iunkn_tur==2) then   !ABL flow, imposes lower limit to epsilon to obtain a maximum mixing length
              if (inv_l_max.gt.0.0001_rp) then ! with coriolis
                 dummr = inv_l_max
              else
                 dummr = 1.0_rp/5000.0_rp ! limiting max mixing length to 5000
              end if
              if (TUR_K_EPS_STD) then 
                 do ipoin =1, npoin   !limits maximum mixing length stronly (lower limit for epsilon)           
                    if ((kfl_fixno_tur(1,ipoin,iunkn_tur)<=0)) & ! .and.walld(ipoin).gt.10.0) &
                         unkno(ipoin) = max(unkno(ipoin), (untur(1,ipoin,1)*untur(1,ipoin,1)*param_tur(6))**0.75_rp*dummr)
                 end do
!!$            else if (TUR_FAMILY_K_OMEGA) then  ! strong limit of maximum mixing length
!!$                 do ipoin =1, npoin              
!!$                    if ((kfl_fixno_tur(1,ipoin,iunkn_tur)<=0)) & ! .and.walld(ipoin).gt.10.0) &
!!$                         unkno(ipoin) = max(unkno(ipoin), sqrt(untur(1,ipoin,1))*(cmu_st**(-0.25_rp))*dummr)
!!$                 end do
              end if
           end if
        else
           do ipoin=1,npoin
              do iunkn_tur=1,nturb_tur
                 itotn=(ipoin-1)*nturb_tur+iunkn_tur
                 if(kfl_fixno_tur(1,ipoin,iunkn_tur)<=0) &
                      unkno(itotn)=max(zezer(iunkn_tur),unkno(itotn))                 
              end do
           end do
        end if
     end if

  else if(kfl_clipp_tur==4) then
     !
     ! Absolute value
     !
     if( INOTMASTER ) then
        if(kfl_algor_tur==1) then
           do ipoin=1,npoin
              if(kfl_fixno_tur(1,ipoin,iunkn_tur)<=0) &
                   unkno(ipoin)=max(zetur,abs(unkno(ipoin)))
           end do
        else
           do ipoin=1,npoin
              itotn=(ipoin-1)*nturb_tur
              do iunkn_tur=1,nturb_tur
                 itotn=itotn+1
                 if(kfl_fixno_tur(1,ipoin,iunkn_tur)<=0) &
                      unkno(itotn)=max(zetur,abs(unkno(itotn)))
              end do
           end do
        end if
     end if

  else if(kfl_clipp_tur==5) then
     !
     ! Last value
     !
     if( INOTMASTER ) then
        if(kfl_algor_tur==1) then
           do ipoin=1,npoin
              if(unkno(ipoin)<=0.0_rp) unkno(ipoin)=untur(iunkn_tur,ipoin,1)
           end do
        else
           do ipoin=1,npoin
              itotn=(ipoin-1)*nturb_tur
              do iunkn_tur=1,nturb_tur
                 itotn=itotn+1
                 if(unkno(itotn)<=0.0_rp) unkno(itotn)=untur(iunkn_tur,ipoin,1)              
              end do
           end do
        end if
     end if

  else if(kfl_clipp_tur==6) then
     !
     ! Maximum value prescribed on inflow
     !
     if(kfl_algor_tur==1) then
        turbm = -1.0_rp
        if( INOTMASTER ) then 
           do ipoin=1,npoin
              if( kfl_fixno_tur(1,ipoin,iunkn_tur) == 6 .or. kfl_fixno_tur(1,ipoin,iunkn_tur) == 8 ) then
                 if( unkno(ipoin) > turbm(1) ) turbm(1) = unkno(ipoin)
              end if
           end do
        end if
        !
        ! Parall: Look for maximum over all subdomains (turbm)
        !
        call pararr('MAX',0_ip,1_ip,turbm(1:))
        if( INOTMASTER ) then 
           do ipoin=1,npoin
!              unkno(ipoin) = max( unkno(ipoin) , 1.0e-6_rp*turbm(1) )
              unkno(ipoin) = max( unkno(ipoin) , clipfac_tur*turbm(1) )
           end do
        end if

     else

!        call runend('TUR_CLIPPI: CLIPPING NOT CODED')

        turbm = -1.0_rp
        if( INOTMASTER ) then 
           do ipoin=1,npoin
              iunkn_tur=1                                  ! kinetic energy
              itotn=(ipoin-1)*nturb_tur+iunkn_tur
              if( kfl_fixno_tur(1,ipoin,iunkn_tur) == 6 .or. kfl_fixno_tur(1,ipoin,iunkn_tur) == 8 ) then
                 if( unkno(itotn) > turbm(1) ) turbm(1) = unkno(itotn)
              end if
           end do
        end if
        !
        ! Parall: Look for maximum over all subdomains (turbm)
        !
        call pararr('MAX',0_ip,1_ip,turbm(1:))
        if( INOTMASTER ) then 
           do ipoin=1,npoin
              iunkn_tur=1                                  ! kinetic energy
              itotn=(ipoin-1)*nturb_tur+ iunkn_tur
              unkno(itotn) = max( unkno(itotn) , clipfac_tur*turbm(1) )
              iunkn_tur=2                                  ! the other one
              itotn=(ipoin-1)*nturb_tur+ iunkn_tur
              unkno(itotn) = max( unkno(itotn) , clipfac_tur*clipfac_tur )
           end do
        end if

     end if

  end if

end subroutine tur_clippi
