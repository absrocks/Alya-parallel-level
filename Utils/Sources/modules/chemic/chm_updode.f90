subroutine chm_updode()
  !------------------------------------------------------------------------
  !****f* partis/chm_updode
  ! NAME 
  !    chm_updode
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_chemic
  implicit none
  integer(ip) :: ipoin,iodes,maxit,iiter,ispec,itave,istat
  real(rp)    :: eps,toler,numer,denom,nume1,deno1
  real(rp)    :: xlhs,xrhs,Cnew
  real(rp)    :: gprea(nodes_chm),gprhs(nodes_chm)
  integer(ip) :: irea1,irea2,ipro1,jodes
  integer(ip) :: iarea,ireac,iskyl,iode1,iode2,iode3
  real(rp)    :: C1,C2,C3,k1,k2,kp,km,ovkT,k1p,k2p,k1m,k2m,deter,T
  real(rp)    :: a(2,2),b(2),c(2,2),xcoef,time1,time2,time3,time4,time5

  if( nodes_chm == 0 ) return

  time1 = 0.0_rp
  time2 = 0.0_rp
  time3 = 0.0_rp
  time4 = 0.0_rp
  time5 = 0.0_rp
  call cputim(time1)

  if( INOTMASTER ) then

     if(1==1) then
        !
        ! Skyline format 
        !
        call chm_usrtem(T)
        ovkT = 1.0_rp/(boltz_chm*T) 

        do ipoin=1,npoin
           !
           ! Initialization
           !
           do iskyl = 1,nskyl_chm
              amatr_chm(iskyl) = 0.0_rp
           end do
           if( kfl_timei_chm == 0 ) then
              do iodes = 1,nodes_chm
                 rhsid_chm(iodes) = 0.0_rp
              end do
           else
              do iodes = 1,nodes_chm
                 iskyl            = idiag_chm(iodes)
                 amatr_chm(iskyl) = dtinv_chm
                 rhsid_chm(iodes) = dtinv_chm * conce(ipoin,iodes+nclas_chm,3)
              end do
           end if
           !
           ! AMATR_CHM, RHSID_CHM: Fill in matrix and RHS
           !
           do iodes = 1,nodes_chm

              ispec = iodes + nclas_chm

              do iarea =  iarea_chm(ispec),iarea_chm(ispec+1)-1
                 ireac =  jarea_chm(iarea)

                 irea1 =  lreac_chm(ireac)%l(1)                                        ! Reactant 1
                 irea2 =  lreac_chm(ireac)%l(2)                                        ! Reactant 2
                 ipro1 = -lreac_chm(ireac)%l(3)                                        ! Product 1

                 iode1 =  irea1-nclas_chm
                 iode2 =  irea2-nclas_chm
                 iode3 =  ipro1-nclas_chm

                 C1    =  conce(ipoin,irea1,1)                                         ! Value Reactant 1
                 C2    =  conce(ipoin,irea2,1)                                         ! Value Reactant 2
                 if( ipro1 /= 0 ) C3 = conce(ipoin,ipro1,1)                            ! Value Product
                 k1    =  0.0_rp
                 k2    =  0.0_rp
                 if( irea1 <= nclas_chm ) &
                      k1 = diffu_chm(1,irea1) * exp(-diffu_chm(2,irea1)*ovkT)                ! k1 = cst1*exp(-E1/kT)
                 if( irea2 <= nclas_chm ) &
                      k2 = diffu_chm(1,irea2) * exp(-diffu_chm(2,irea2)*ovkT)                ! k2 = cst2*exp(-E2/kT)
                 kp    = 4.0_rp * pi * radiu_chm(ireac) * (k1+k2)                            ! k+ = 4 * pi * (rA+rB)*(kA+kB)
                 km    = denma_chm * kp * react_chm(1,ireac) * exp(-react_chm(2,ireac)*ovkT) ! k- = rho * k+ * exp(-E/kT)

                 if( irea1 == ispec ) then                                  ! I am reactant 1
                    jodes = iode1
                    xcoef = kp*C2
                    call chm_updsky(iodes,jodes,xcoef)
                    if( iode3 <= 0 ) then                                   ! Product is EDP
                       rhsid_chm(iodes) = rhsid_chm(iodes) + km*C3
                    else                                                    ! Product is ODE
                       jodes = iode3
                       xcoef = -km
                       call chm_updsky(iodes,jodes,xcoef)
                    end if
                 end if

                 if( irea2 == ispec ) then                                  ! I am reactant 2
                    jodes = iode2
                    xcoef = kp*C1
                    call chm_updsky(iodes,jodes,xcoef)
                    if( iode3 <= 0 ) then                                   ! Product is EDP
                       rhsid_chm(iodes) = rhsid_chm(iodes) + km*C3
                    else                                                    ! Product is ODE
                       jodes = iode3
                       xcoef = -km
                       call chm_updsky(iodes,jodes,xcoef)
                    end if
                 end if

                 if( ipro1 == ispec ) then                                   ! I am product
                    jodes = iode3
                    xcoef = km
                    call chm_updsky(iodes,jodes,xcoef)
                    if( iode1 <= 0 .and. iode2 <= 0 ) then                   ! Product is EDP
                       rhsid_chm(iodes) = rhsid_chm(iodes) + kp*C1*C2
                    else if( iode1 > 0 ) then                                ! Product is ODE
                       jodes = iode1
                       xcoef = -kp*C2
                       call chm_updsky(iodes,jodes,xcoef)
                    else
                       jodes = iode2
                       xcoef = -kp*C1
                       call chm_updsky(iodes,jodes,xcoef)
                    end if

                 else if( ipro1 == 0 .and. iode1 > 0 .and. iode2 > 0 ) then  ! No specy generated
                    call runend('CHM_UPDODE: NOT CODED')

                 end if

              end do
           end do
           !
           ! CONCE(IPOIN,:,1): Solve system and update concentration
           !
           call cputim(time3)
           call lufact(&
                nodes_chm, nskyl_chm, iskyl_chm, amatr_chm, idiag_chm, istat )
           if(istat/=0) call runend('COULD NOT FACTORIZE ODE MATRIX')
           call lusolv(&
                nodes_chm, nskyl_chm, iskyl_chm, 1_ip , amatr_chm, rhsid_chm,&
                nodes_chm, idiag_chm, istat )
           call cputim(time4)
           time5 = time5 + ( time4 - time3 )
           if(istat/=0) call runend('COULD NOT SOLVE ODE SYSTEM')
           ispec = nclas_chm
           do iodes = 1,nodes_chm
              ispec = ispec + 1
              conce(ipoin,ispec,1) = max(0.0_rp,rhsid_chm(iodes))
           end do

        end do

     else if(1==0) then
        !
        ! Direct (only for 2 reactions)
        !
        call chm_usrtem(T)
        ovkT  =  1.0_rp/(boltz_chm*T) 
        ireac =  1
        irea1 =  lreac_chm(ireac)%l(1)                                           ! Reactant 1
        irea2 =  lreac_chm(ireac)%l(2)                                           ! Reactant 2
        ipro1 = -lreac_chm(ireac)%l(3)                                           ! Product 1
        k1    =  0.0_rp
        k2    =  0.0_rp
        if(irea1<=nclas_chm)&
             k1 = diffu_chm(1,irea1)*exp(-diffu_chm(2,irea1)*ovkT)               ! k1 = cst1*exp(-E1/kT)
        if(irea2<=nclas_chm)&
             k2 = diffu_chm(1,irea2)*exp(-diffu_chm(2,irea2)*ovkT)               ! k2 = cst2*exp(-E2/kT)
        k1p   = 4.0_rp*pi*radiu_chm(ireac)*(k1+k2)                               ! k+ = 4 * pi * (rA+rB)*(kA+kB)
        k1m   = denma_chm*k1p*react_chm(1,ireac)*exp(-react_chm(2,ireac)*ovkT)   ! k- = rho * k+ * exp(-E/kT)

        ireac =  2
        irea1 =  lreac_chm(ireac)%l(1)                                           ! Reactant 1
        irea2 =  lreac_chm(ireac)%l(2)                                           ! Reactant 2
        if( ipro1 /= 0 ) &
             ipro1 = -lreac_chm(ireac)%l(3)                                      ! Product 1
        k1    =  0.0_rp
        k2    =  0.0_rp
        if(irea1<=nclas_chm)&
             k1 = diffu_chm(1,irea1)*exp(-diffu_chm(2,irea1)*ovkT)               ! k1 = cst1*exp(-E1/kT)
        if(irea2<=nclas_chm)&
             k2 = diffu_chm(1,irea2)*exp(-diffu_chm(2,irea2)*ovkT)               ! k2 = cst2*exp(-E2/kT)
        k2p   = 4.0_rp*pi*radiu_chm(ireac)*(k1+k2)                               ! k+ = 4 * pi * (rA+rB)*(kA+kB)
        k2m   = denma_chm*k2p*react_chm(1,ireac)*exp(-react_chm(2,ireac)*ovkT)   ! k- = rho * k+ * exp(-E/kT)

        do ipoin=1,npoin
           C1     =  conce(ipoin,1,1)
           a(1,1) =  dtinv_chm + k1m + k2p*C1
           a(1,2) = -k2m 
           b(1)   =  k1p*C1*C1 + dtinv_chm*conce(ipoin,2,3)
           a(2,1) = -k2p*C1
           a(2,2) =  dtinv_chm + k2m
           b(2)   =  0.0_rp + dtinv_chm*conce(ipoin,3,3)
           call invmtx(a,c,deter,2_ip)
           conce(ipoin,2,1)=c(1,1)*b(1)+c(1,2)*b(2)
           conce(ipoin,3,1)=c(2,1)*b(1)+c(2,2)*b(2)
        end do

     end if

  end if

  call cputim(time2)
  cputi_chm(5) = cputi_chm(5) + time5
  cputi_chm(4) = cputi_chm(4) + (time2-time1)-time5 

end subroutine chm_updode

subroutine chm_updsky(iodes,jodes,xcoef)

  use def_kintyp, only    :  ip,rp
  use def_chemic, only    :  amatr_chm,iskyl_chm,idiag_chm
  integer(ip), intent(in) :: iodes,jodes
  real(rp),    intent(in) :: xcoef
  integer(ip)             :: kskyl

  if(iodes<jodes) then
     kskyl=iskyl_chm(jodes+1)-(jodes-iodes)
     amatr_chm(kskyl)=amatr_chm(kskyl)+ xcoef
  else     
     kskyl=idiag_chm(iodes)-(iodes-jodes)
     amatr_chm(kskyl)=amatr_chm(kskyl)+ xcoef
  end if
  
end subroutine chm_updsky
