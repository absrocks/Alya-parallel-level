subroutine chm_odepro(ipoin,iode1,nodes,gprea,gprhs)
  !------------------------------------------------------------------------
  !****f* partis/chm_odepro
  ! NAME 
  !    chm_odepro
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_parame, only     :  pi
  use def_domain, only     :  npoin
  use def_master, only     :  conce
  use def_chemic, only     :  nodes_chm,iarea_chm,jarea_chm,&
       &                      diffu_chm,radiu_chm,denma_chm,&
       &                      nclas_chm,lreac_chm,react_chm,&
       &                      boltz_chm
  implicit none
  integer(ip), intent(in)  :: ipoin,iode1,nodes
  real(rp),    intent(out) :: gprea(nodes)
  real(rp),    intent(out) :: gprhs(nodes)
  integer(ip)              :: iodes,irea1,irea2,ipro1,kodes
  integer(ip)              :: iarea,ireac,ii,ispec
  real(rp)                 :: C1,C2,C3,k1,k2,kp,km,ovkT,T

  call chm_usrtem(T)
  ovkT = 1.0_rp/(boltz_chm*T) 

  do ii=1,nodes
     gprea(ii) = 0.0_rp
     gprhs(ii) = 0.0_rp
  end do

  do ii=1,nodes

     iodes = ii+iode1                                                 ! Number of the ODE
     ispec = iodes+nclas_chm                                          ! Number of specy

     do iarea =  iarea_chm(ispec),iarea_chm(ispec+1)-1
        ireac =  jarea_chm(iarea)                                     ! Reaction ireac 
        irea1 =  lreac_chm(ireac)%l(1)                                ! Reactant 1
        irea2 =  lreac_chm(ireac)%l(2)                                ! Reactant 2
        ipro1 = -lreac_chm(ireac)%l(3)                                ! Product 1
        C1    =  conce(ipoin,irea1,1)                                 ! Value Reactant 1
        C2    =  conce(ipoin,irea2,1)                                 ! Value Reactant 2
        C3    =  conce(ipoin,ipro1,1)                                 ! Value Reactant 2
        k1    =  0.0_rp
        k2    =  0.0_rp
        if(irea1<=nclas_chm)&
             k1 = diffu_chm(1,irea1)*exp(-diffu_chm(2,irea1)*ovkT)    ! k1 = cst1*exp(-E1/kT)
        if(irea2<=nclas_chm)&
             k2 = diffu_chm(1,irea2)*exp(-diffu_chm(2,irea2)*ovkT)    ! k2 = cst2*exp(-E2/kT)
        kp    = 4.0_rp*pi*radiu_chm(ireac)*(k1+k2)                    ! k+ = 4 * pi * (rA+rB)*(kA+kB)
        km    = denma_chm*kp*react_chm(1,ireac)*exp(-react_chm(2,ireac)*ovkT) ! k- = rho * k+ * exp(-E/kT) 
        !
        ! 3 cases
        ! -------
        !
        ! I1 + I2 <=> I3
        ! dC1/dt + (k+ * C2) * C1 = k- * C3
        !
        ! I1 + I1 <=> I2
        ! dC1/dt + (2 * k+ * C1) * C1 = 2 * k- * C2
        !
        ! I2 + I3 <=> I1
        ! dC1/dt - (k-) * C1 = - k+ * C2 * C3 
        !
        if(irea1==ispec) then                                         ! I am reactant 1
           gprea(ii) = gprea(ii) + kp*C2
           gprhs(ii) = gprhs(ii) + km*C3
        end if
        if(irea2==ispec) then                                         ! I am reactant 2
           gprea(ii) = gprea(ii) + kp*C1
           gprhs(ii) = gprhs(ii) + km*C3
        end if
        if(ipro1==ispec) then                                         ! I am product
           gprea(ii) = gprea(ii) + km
           gprhs(ii) = gprhs(ii) + kp*C1*C2
        end if

     end do

  end do

end subroutine chm_odepro
