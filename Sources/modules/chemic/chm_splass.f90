subroutine chm_splass()
  !-----------------------------------------------------------------------
  !****f* partis/chm_splass
  ! NAME 
  !    chm_splass
  ! DESCRIPTION
  !    This routine assembles matrix and RHS
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_chemic
  use def_domain
  use def_solver
  implicit none
  integer(ip), save :: ipass=0
  integer(ip)       :: ipoin,izdom,iclas,iarea,ireac
  integer(ip)       :: irea1,irea2,ipro1,idime,kzdom
  integer(ip)       :: izmat,kzmat,izrhs,kzrhs,nsmat,nsrhs
  real(rp)          :: k,kp,km,k1,k2,ovkT,C1,C2,C3,xfact,adiag
  real(rp)          :: C0eq1,Eform1,Ceq1,C0eq2,Eform2,Ceq2,T

  if( kfl_assem_chm >= 2 ) then
     !----------------------------------------------------------------------
     !
     ! Save constant matrix
     !
     !----------------------------------------------------------------------

     if( kfl_assem_chm >= 3 ) then

        nsmat = solve_sol(1) % nzmat 
        nsrhs = solve_sol(1) % nzrhs
        kzmat = ( iclas_chm - 1 ) * nsmat 
        kzrhs = ( iclas_chm - 1 ) * nsrhs 

        if( ipass == 0 ) then           
           do izmat = 1,nsmat
              kzmat = kzmat + 1
              smatr_chm(kzmat) = amatr(izmat) 
           end do
           do izrhs = 1,nsrhs
              kzrhs = kzrhs + 1
              shsid_chm(kzrhs) = rhsid(izrhs) 
           end do

        else 
           do izmat = 1,nsmat
              kzmat = kzmat + 1
              amatr(izmat) = smatr_chm(kzmat)
           end do
           do izrhs = 1,nsrhs
              kzrhs = kzrhs + 1
              rhsid(izrhs) = shsid_chm(kzrhs)
           end do
           kfl_assem_chm = 4

        end if
        ipass = ipass + 1

     end if

     !----------------------------------------------------------------------
     !
     ! Constants
     !
     !----------------------------------------------------------------------

     call chm_usrtem(T)
     ovkT  = 1.0_rp/(boltz_chm*T) 
     iclas = iclas_chm

     !----------------------------------------------------------------------
     !
     ! Reaction
     !
     !----------------------------------------------------------------------

     do ipoin = 1,npoin

        izdom = idima_chm(ipoin) 
        do iarea =  iarea_chm(iclas),iarea_chm(iclas+1)-1

           ireac =  jarea_chm(iarea)          ! Reaction ireac 
           irea1 =  lreac_chm(ireac) %l (1)   ! Reactant 1
           irea2 =  lreac_chm(ireac) %l (2)   ! Reactant 2
           ipro1 = -lreac_chm(ireac) %l (3)   ! Product 1

           C1 =  conce(ipoin,irea1,1) 
           C2 =  conce(ipoin,irea2,1) 
           if( ipro1 /= 0 ) C3 = conce(ipoin,ipro1,1) 

           k1 =  0.0_rp
           k2 =  0.0_rp

           if( irea1 <= nclas_chm ) &                                                  ! k1 = k1'*exp(-E1/kT)
                k1 = diffu_chm(1,irea1) * exp(-diffu_chm(2,irea1)*ovkT)
           if( irea2 <= nclas_chm ) &                                                  ! k2 = k2'*exp(-E2/kT)
                k2 = diffu_chm(1,irea2) * exp(-diffu_chm(2,irea2)*ovkT)

           kp = 4.0_rp * pi * radiu_chm(ireac) * (k1+k2)                               ! k+ = 4 * pi * (rA+rB)*(kA+kB)
           km = denma_chm * kp * react_chm(1,ireac) * exp(-react_chm(2,ireac) * ovkT)  ! k- = rho * k+ * exp(-E/kT)
           kp = kp * vmass_chm(ipoin)
           km = km * vmass_chm(ipoin)

           if( ipro1 /= 0 ) then
              !
              ! Reaction has a product: R1 + R2 <=> P1
              !
              if( irea1 == iclas ) then                                  ! I am reactant 1
                 amatr(izdom) = amatr(izdom) + kp * C2 
                 rhsid(ipoin) = rhsid(ipoin) + km * C3
              end if

              if( irea2 == iclas ) then                                  ! I am reactant 2
                 amatr(izdom) = amatr(izdom) + kp * C1
                 rhsid(ipoin) = rhsid(ipoin) + km * C3
              end if

              if( ipro1 == iclas ) then                                  ! I am product
                 amatr(izdom) = amatr(izdom) + km
                 rhsid(ipoin) = rhsid(ipoin) + kp * C1 * C2
              end if

           else if( ipro1 == 0 ) then
              !
              ! Reaction does not have a product: R1 + R2 <=> 0
              !
              C0eq1  = equil_chm(1,irea1)                                ! C0eq
              Eform1 = equil_chm(2,irea1)                                ! E_form: formation energy
              Ceq1   = C0eq1 * exp(-Eform1*ovkT)                         ! Ceq = C0eq * exp(-E_form/kT)

              C0eq2  = equil_chm(1,irea2)                                ! C0eq
              Eform2 = equil_chm(2,irea2)                                ! E_form: formation energy
              Ceq2   = C0eq2 * exp(-Eform2*ovkT)                         ! Ceq = C0eq * exp(-E_form/kT)

              C3     = Ceq1 * Ceq2

              if( irea1 == iclas ) then                                  ! I am reactant 1
                 amatr(izdom) = amatr(izdom) + kp * C2
                 rhsid(ipoin) = rhsid(ipoin) - kp * C3
              end if

              if( irea2 == iclas ) then                                  ! I am reactant 2
                 amatr(izdom) = amatr(izdom) + kp * C1
                 rhsid(ipoin) = rhsid(ipoin) - kp * C3
              end if

           end if

        end do

     end do

     !----------------------------------------------------------------------
     !
     ! Matrix-based time step
     !
     !----------------------------------------------------------------------

     dtmat_chm = 0.0_rp
     do ipoin = 1,npoin
        izdom = idima_chm(ipoin) 
        dtmat_chm = max(dtmat_chm,amatr(izdom)/vmass_chm(ipoin))
     end do
     dtmat_chm = 1.0_rp/dtmat_chm

     !----------------------------------------------------------------------
     !
     ! Time
     !
     !----------------------------------------------------------------------

     if( kfl_timei_chm > 0 ) then
        do ipoin = 1,npoin
           izdom = idima_chm(ipoin) 
           xfact = vmass_chm(ipoin) * dtinv_chm
           amatr(izdom) = amatr(izdom) + xfact
           rhsid(ipoin) = rhsid(ipoin) + xfact * conce(ipoin,iclas,3)
        end do
     end if

     !----------------------------------------------------------------------
     !
     ! Source
     !
     !----------------------------------------------------------------------

     if( kfl_sourc_chm == 1 ) then
        !
        ! Constant source
        !
        do ipoin = 1,npoin
           rhsid(ipoin) = rhsid(ipoin) + vmass_chm(ipoin) * sourc_chm
        end do

     else if( kfl_sourc_chm == 2 ) then
        !
        ! Sphere source
        !
        do ipoin = 1,npoin
           k = 0.0_rp
           do idime = 1,ndime
              k = k + ( coord(idime,ipoin) - coord(idime,ipoin) ) ** 2
           end do
           if( sqrt(k) <= sorad_chm ) then
              rhsid(ipoin) = rhsid(ipoin) + vmass_chm(ipoin) * sourc_chm
           end if
        end do

     else if( kfl_sourc_chm < 0 ) then
        !
        ! Meteo only
        !
        do ipoin = 1,npoin
           rhsid(ipoin) = rhsid(ipoin) + vmass_chm(ipoin) * tmrat_chm(iclas_chm,ipoin)
        end do

     end if

     !----------------------------------------------------------------------
     !
     ! Boundary conditions
     !
     !----------------------------------------------------------------------

     do ipoin = 1,npoin
        if(  kfl_fixno_chm(iclas,ipoin) >= 1 ) then
           kzdom = idima_chm(ipoin) 
           adiag = amatr(kzdom)
           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
              amatr(izdom) = 0.0_rp
           end do
           amatr(kzdom) = adiag
           rhsid(ipoin) = adiag * bvess_chm(iclas,ipoin)
        end if
     end do

  end if

end subroutine chm_splass

