   subroutine solite(iunkn)
     !-----------------------------------------------------------------------
     !****f* LowMac/lom_solite
     ! NAME 
     !    lom_solite
     ! DESCRIPTION
     !    This routine solves an iteration for the iunkn equation
     !
     !  iunkn = 1 : x-Velocity
     !  iunkn = 2 : y-Velocity
     !  iunkn = 3 : kinetic turb
     !  iunkn = 4 : epsilon turb
     !  iunkn = 5 : temper,  only for transient problem (initial solution must be given for all fields)
     ! USED BY
     !    Doiter
     !***
     !-----------------------------------------------------------------------
     
     use def_master
     implicit none

     integer(4), intent(in)::  iunkn
     real(8)               ::  wmatr(nnode, nnode), wrhsi(nnode), elmat(nnode, nnode),elrhs(nnode)
     real(8)               ::  chale, detjm, dvolu, elcod(nnode), elvel(nnode,2,2), elkey(nnode, 3)
     real(8)               ::  eleps(nnode, 3), cartd(nnode), gpkey(3), gpeps(3), gpvel(2,2), gptem(2), eltem(nnode, 2)
     real(8)               ::  gpcod, wvalu, tvalu, bvalu(2), shave(nnode), grvel(2), adiag, hessi(nnode), hevel(2)
     real(8)               ::  gpmut, produ, diffu, react, force, tract, vnorm, tnume,tdeno, prodm, prodt
     real(8)               ::  error, c1p, lm, epsil_aux(npoin), key_aux(npoin), pert, tau, resid
     integer               ::  itera,  maxin(5), inico, finco, index, iz, jz, kz
     real(8)               ::  B, ui, dtrac, accum, grmut, eta, DC, W, uaste, A0, As, grkey, greps
     real(8)               ::  grtem, linke, linre, rfact , C3, zbyLm, F, qwall, sqkey, Cg, alpha, richg, difkw, modve
     real(8)               ::  percn, dracn, LADEN, n, zm, lam, sp, sd, betap, betad, C4, C5, &
                               isone,  ladta(11), codta(11), facto,cmuch,  umtar, gpdhf(nnode), S
     logical               :: flctr
     real(8)               :: phi_m, phi_h, phi_e, logfu, hconv, htflx_can  !, logft
     real(8), save         :: prgra,  umean(2), divhf(4000)
     real(8)               :: gpvis = 1.5e-5
     !NEW variables for GABLS3 case type (implemented by JBR)
     real(8)               ::  elu_adv(nnode), elv_adv(nnode),eltempe_adv(nnode),elu_geo(nnode), elv_geo(nnode)
     real(8)               ::  elu_meso(nnode), elv_meso(nnode), elth_meso(nnode)
     real(8)               ::  gptem_adv,gpu_adv,gpv_adv,gpu_geo,gpv_geo
     real(8)               ::  gpu_meso,gpv_meso, gpth_meso
     real(8)               ::  eta_can = 0.6 ! extintion coefficient
     real(8)               ::  perio, a_damp, leng_damp
     perio = 24.0*3600.0
     htflx_can =200.0/rhocp
   !  htflx_can = 0.2
!     htflx_can = 200.0/rhocp*sin(2.0*pi/perio*ctime)
   !  call interp_radiat(htflx_can) ;  htflx_can = - htflx_can /rhocp
    
   !  call interp_ugeos(ugeos)
     ! print *,  ctime, htflx_can
     !
     !  ***intialization of variables
     !
     flctr = .false. ! flow control for canopy ( imposing pressure gradient) 
     linre = 2.0d0  ! linearization, 2 : N-R, 1 : implicit, 0: explicit
     linke = 2.0d0  ! linearization, 2 : N-R, 1 : implicit, 0: explicit
     error = 100.0d0
     pert  = 0.0d0
     DC    = 0.0d0 ! 
     Cg    = 0.0d0 ! thermal coefficient
     itera = 1  ! inner iteration counter
     maxin(1) = 1
     maxin(2) = 1 ! Linear equations (unless frontline)
     maxin(3) = 5  !maxit(3)
     maxin(4) = 5  !maxit(4) ! Linear equations
     maxin(5) = 1 ! temp equation, linear equation.
     hessi = 0.0d0
     prodt = 0.0d0
     difkw = 0.0d0
     betap = 0.0d0
     betad = 0.0d0
     C4    = 0.0d0
     C5    = 0.0d0
     isone = 1.0d0
     !Thermal stability Monin-Obukhov correction factors   
     phi_e   = 1.0d0
     phi_m   = 1.0d0
     phi_h=mo_prand
     !NEW variables for GABLS3 case type (implemented by JBR)
     eltempe_adv =0.0d0
     elu_adv = 0.0d0
     elv_adv = 0.0d0
     elu_geo = 0.0d0
     elv_geo = 0.0d0
     elu_meso = 0.0d0
     elv_meso = 0.0d0


     if (istep .eq.  1 ) divhf= 0.0d0
     
     if ((.not.flctr).or.(.not.kfl_canop)) &
          prgra =0.0d0 ! driving pressure gradient (Canopy model) 


     umtar = 2.0d0  ! umean target ( valid for control flow for canopy )
     if  (kfl_canop) then 
     !
     ! LOADS CANOPY MODEL
     !    
        if ((istep.eq.1.or..not.flctr).and.iunkn==1) &
        prgra = 0.00d0
        
        if (kfl_canmo == 0) then ! default (Sogachev's)
           betap = 0.0d0          
           betad = 12.0d0*sqrt(cmu)
           C4    = 0.0d0
           C5    = C1-C2          
           isone = 0.0d0 !not in k equation
        else if (kfl_canmo == 1) then ! Svenson
           betap = 1.0d0
           betad = 0.0d0
           C4    = 1.95d0
           C5    = 0.0d0
        else if (kfl_canmo == 2) then ! Sanz
           betap = 1.0d0
           c4    = sigka*(2.0d0/sigep-((c2-c1)*sqrt(cmu)/6.0d0) *(2.0d0/0.05d0)**(2.0d0/3.0d0))
           C5    = C4           
           betad = sqrt(cmu)*((2.0/0.05d0)**(2.0d0/3.0d0))*betap + 3.0d0/sigka
!           C4    = 0.9d0
!           C5    = 0.9d0
!           betad=  5.1d0
!           print *, 'betad, c4', betad, c4
!           stop
        else if (kfl_canmo == 3) then ! Green
           betap = 1.0d0
           betad = 4.0d0
           C4    = 1.5d0
           C5    = 1.5d0
        else if (kfl_canmo == 4) then !Lopes da costa
           betap = 0.17d0
           betad = 3.37d0
           C4    = 0.9d0
           C5    = 0.9d0
        else
           call runend('ERROR:NO CANOPY MODEL')
        end if           

     end if
     !
     if (iunkn.ge.0) then
        do ipoin=1, npoin
           key_aux(ipoin)=keyva(ipoin,1)
           epsil_aux(ipoin)=epsil(ipoin,1)
        end do
     end if
     !
     if (kfl_thmod.eq.1) then
        !
        !Monin-Obukhov length evaluation
        gpmut =  cmu*keyva(1,1)*keyva(1,1)/epsil(1,1)        
        qwall =  rhocp*gpmut*(tempe(2,1)- tempe(1,1))/(coord(2)-coord(1))/sigte ! qwall positive in stable atm
        if (kfl_bouco_vel.eq.2) then ! wall law
           if (.not.kfl_abl2) ustar2=ustar
           qwall = rhocp*ustar2*kar/(mo_prand*log(1.0+dwall/rough))* (tempe(1,1)-tewal)
        end if
        lmoni =  rhocp*teref*ustar*ustar*ustar/(kar*gravi*qwall)
        !
        !Monin-Obukhov length clipping (implemented by JBR)
        !if (lmoni.lt.-1.0d-6.and.lmoni.gt.-1.0d0) lmoni = -1.0d0
        !if (lmoni.gt.1.0d-6.and.lmoni.lt.1.0d0)   lmoni = 1.0d0   
        !
        !Mixmimum mixing length and Ambient dissiaption rate 
        l_max= lenmy
        !epsam =  ((cmu*keyam*keyam)**(0.75))/lenmy !(implemented by JBR) 
     end if
     
     !  open(51, file='Wind55.plot')
     ! Realizable model 
     if (kfl_model==2) then
        A0=(1.0d0-3.0d0*sqrt(0.5d0*cmu0))/cmu0
        cmuch = sqrt(cmu/0.09d0) ! change factor for cmu
     end if
     
     !
     !Iterations for solving the model 
     !
     do while (itera.le.maxin(iunkn).and.error.gt.toler(iunkn))
        !
        ! Construct the system matrix and right-hand-side.
        !
        amatr=0.0d0  ! matrix initialization
        rhsid=0.0d0  ! RHS initialization
        
        ! Loop over elements
        elements: do ielem=1,nelem
           
           ! Initializations
           elmat=0.0d0 !elemental matrix
           elrhs=0.0d0 !elemental rhs      
           gpdhf=0.0d0 !heat flux divergence

           ! Gather operations     
           do inode =1, nnode
              ipoin = lnods(inode, ielem)
              elcod(inode)=coord(ipoin)
              do itime =1, 2  
                 elvel(inode,1,itime) = veloc(ipoin,1,itime)
                 elvel(inode,2,itime) = veloc(ipoin,2,itime)
                 elkey(inode,itime)=keyva(ipoin, itime)
                 eleps(inode,itime)=epsil(ipoin, itime)
                 eltem(inode,itime)=tempe(ipoin, itime)
              end do
              elkey(inode,3) = key_aux(ipoin)
              eleps(inode,3) = epsil_aux(ipoin)
           end do
           ! Mesoscale advections/tendencies (implemented by JBR)  
           if (kfl_momadv) then               
              do inode =1, nnode
                 ipoin = lnods(inode, ielem)
                 elu_adv(inode) = u_adv(ipoin)
                 elv_adv(inode) = v_adv(ipoin)
              end do
           end if
           if (kfl_thadv) then
              do inode =1, nnode
                 ipoin = lnods(inode, ielem)
                 eltempe_adv(inode) = tempe_adv(ipoin)
              end do
           end if
           if (kfl_pressgr) then           
              do inode =1, nnode
                 ipoin = lnods(inode, ielem)
                 elu_geo(inode) = u_geo(ipoin)
                 elv_geo(inode) = v_geo(ipoin)
              end do
           end if
           !Nudging to mesoscale fields (implemented by JBR)  
           if (kfl_mom_nudging) then 

              do inode =1, nnode
                 ipoin = lnods(inode, ielem)
                 elu_meso(inode) = u_meso(ipoin)
                 elv_meso(inode) = v_meso(ipoin)
              end do
           end if
           if (kfl_tem_nudging) then
              do inode =1, nnode
                 ipoin = lnods(inode, ielem)
                 elth_meso(inode) = th_meso(ipoin)
              end do
           end if        

           ! Compute the characteristic length chale     
           chale = elcod(2) -elcod(1)     

           ! Cartesian derivates and jacobian (only valid for linear elements)
           cartd(1) = - 1.0d0/chale  
           cartd(2) =   1.0d0/chale
           detjm   =    chale*0.5d0

           gauss_points: do igaus=1,ngaus
              wmatr=0.0d0
              wrhsi=0.0d0

              ! derivatives and jacobian high order
              if (kfl_order==2) then  ! second order
                 !Jacobian at gaus point
                 detjm=0.0d0  ! dX/de
                 do inode =1, nnode
                    detjm = detjm + deriv(inode, igaus)*elcod(inode)
                 end do
                 ! first order derivative
                 do inode =1, nnode
                    cartd(inode) =  deriv(inode,igaus)/detjm
                 end do

                 ! Second order derivatives
                 accum =0.0d0  ! d2x/de2
                 do inode =1, nnode
                    accum = accum + deri2(inode, igaus)*elcod(inode)
                 end do
                 ! d2Na/dx2
                 do inode =1, nnode
                    hessi(inode) = (deri2(inode, igaus) - deriv(inode, igaus)*accum/detjm)/detjm/detjm
                 end do
              end if
              dvolu=weigp(igaus)*detjm
              shave(1:nnode) = shape(1:nnode, igaus)

              ! Interpolation (Gauss point values)
              gpvel     = 0.0d0 !gauss point velocity     
              gpkey     = 0.0d0 !gauss point key
              gpeps     = 0.0d0 !gauss point epsilon
              gpcod     = 0.0d0 !gauss point coordinate
              grvel     = 0.0d0 !velocity gradient
              hevel     = 0.0d0 !velocity hessian
              grmut     = 0.0d0 !turbulent viscosity gradient
              gptem     = 0.0d0 !gauss point temperature
              grtem     = 0.0d0 !temperature gradient
              grkey     = 0.0d0 !tke gradient
              greps     = 0.0d0 !dissipation gradient
              !Mesoscale variables  (implemented by JBR)  
              gptem_adv = 0.0d0 !Temperature advection
              gpu_adv   = 0.0d0 !U vel component advection
              gpv_adv   = 0.0d0 !V vel component advection
              gpu_geo   = 0.0d0 !U vel geostrophic component
              gpv_geo   = 0.0d0 !V vel geostrophic component
              gpu_meso  = 0.0d0 !U vel mesoscale component for nudging
              gpv_meso  = 0.0d0 !V vel mesoscale component for nudging
              gpth_meso = 0.0d0 !Th mesoscale for nudging
              
              do inode =1,nnode      
                 do itime =1, 2 ! now and previous time step
                    gpvel(1,itime)= gpvel(1,itime) +  shave(inode)*elvel(inode,1,itime)
                    gpvel(2,itime)= gpvel(2,itime) +  shave(inode)*elvel(inode,2,itime)
                    gpkey (itime) = gpkey(itime)   +  shave(inode)*  elkey(inode,itime)
                    gpeps (itime) = gpeps(itime)   +  shave(inode)*  eleps(inode,itime)       
                    gptem (itime) = gptem(itime)   +  shave(inode)*  eltem(inode,itime)       
                 end do
                 gpkey(3) = gpkey(3) + shave(inode)*  elkey(inode,3)
                 gpeps(3) = gpeps(3) + shave(inode)*  eleps(inode,3)
                 gpcod   = gpcod     + shave(inode)*  elcod(inode)
                 grvel(1)= grvel(1)  + elvel(inode, 1 ,1)*cartd(inode)
                 grvel(2)= grvel(2)  + elvel(inode, 2 ,1)*cartd(inode)
                 grkey  =  grkey     + elkey(inode, 1)*cartd(inode)
                 greps =   greps     + eleps(inode, 1)*cartd(inode)
                 grtem   = grtem     + eltem(inode,1)*cartd(inode)
                 hevel(1)= hevel(1) + elvel(inode, 1 ,1)*hessi(inode)
                 hevel(2)= hevel(2) + elvel(inode, 2 ,1)*hessi(inode)
                 grmut = grmut + densi*cmu*elkey(inode,1)*elkey(inode,1)/eleps(inode,1)*cartd(inode)
                 ! Advections  (implemented by JBR)  
                 gptem_adv = gptem_adv + shave(inode)*eltempe_adv(inode)
                 gpu_adv   = gpu_adv   + shave(inode)*elu_adv(inode)
                 gpv_adv   = gpv_adv   + shave(inode)*elv_adv(inode)
                 ! Geostrophic Pressure gradient along vertical  (implemented by JBR)  
                 gpu_geo   = gpu_geo   + shave(inode)*elu_geo(inode)
                 gpv_geo   = gpv_geo   + shave(inode)*elv_geo(inode)
                 ! Mesoscale Wind and temperature variable along vertical  (implemented by JBR)  
                 gpu_meso   = gpu_meso   + shave(inode)*elu_meso(inode)
                 gpv_meso   = gpv_meso   + shave(inode)*elv_meso(inode)
                 gpth_meso  = gpth_meso  + shave(inode)*elth_meso(inode)
              end do

              !
              ! Coefficients
              !
              ! REALIZABLE MODEL
              if (kfl_model==2) then !Realizable model               
                 uaste = sqrt(grvel(1)*grvel(1)+grvel(2)*grvel(2))
                 W=0.0d0
                 As=1.5d0*sqrt(2.0d0)
                 cmu = 1.0d0/(A0+As*uaste*gpkey(3)/gpeps(3))    
                 grmut = 0.0d0
              end if

              ! this is for k -l model
              !lm =kar*(gpcod+rough)*l_max/(kar*(gpcod+rough)+l_max)
              !gpeps(1) = densi/lm*(cmu*gpkey(3)*gpkey(3))**(0.75)   
              !gpmut = lm*(cmu*gpkey(1)*gpkey(1))**0.25d0

              ! turbulence viscosity
              if (kfl_logva) then
                 gpmut = max(0.00001d0,densi*cmu*exp(2.0*gpkey(3)-gpeps(3)))
              else
                 gpmut = max (0.00001d0, densi*cmu*gpkey(3)*gpkey(3)/gpeps(3)) 
                 ! if thermally coupled, turbulent visco > 1.0 
              end if
              gpmut = min (gpmut, 1.0e30)

              ! mechanical Production
              prodm = gpmut*(grvel(1)*grvel(1)+grvel(2)*grvel(2))
              !           if (iunkn.ge.3) then
              if (kfl_thcou) then     ! thermal coupling (stationary)
                 !thermal production
                 prodt = - gravi/teref*gpmut*grtem/sigte
                 if (kfl_thmod.eq.1)  then ! modification of sigte and prodt
                    if (kfl_trtem) l_max= lenmy !reduntant (JBR)
                    if (prodt.lt.1.0e-10) then ! stable
                       sigte = 0.74d0 
                       prodt = - gravi/teref*gpmut*grtem/sigte
                    else                   ! unstable
                       lm = (cmu*gpkey(1)*gpkey(1))**(0.75d0) /gpeps(3) 
                       if (lm.lt.0.0d0) lm=l_max
                       alpha =1.0d0-(1.0d0+(C2-1.0d0)/(C2-C1) )*lm/l_max
                       richg = - prodt/(prodm+abs(alpha/sigte*prodt))
                       if (.not.kfl_canop) then  ! modify prandtl to satisfy MO distribution     
                          sigte = 0.74d0*(1.0d0-15.0d0*richg)**(-0.25d0)
                          prodt = - gravi/teref*gpmut*grtem/sigte
                          richg = - prodt/(prodm+abs(alpha/sigte*prodt))
                          sigte = 0.74d0*(1.0d0-15.0d0*richg)**(-0.25d0)
                          prodt = - gravi/teref*gpmut*grtem/sigte
                       end if
                    end if
                 end if
                 !                    prodt = - gravi/298.0*hflx0/rhocp    
                 ! modification in the model due to thermal coupling
                 !                    beta = 0.0d0
                 !                    produ =  produ - beta*(C2-C1)/C1*prodt
                 !                    produ = max(1.0d-14,produ)                       
                 !                   produ = prodm + prodt               
                 !                   prodm = prodm*(1.0+prodt/prodm)
                 !                   gpmut = gpmut*(1.0+prodt/prodm)
                 produ = max(0.0d0,prodm+ prodt)
                 prodt = produ - prodm
                 
              else  ! no  thermal coupling 
                 produ = prodm
              end if

              ! velocity module, for canopy flows
              modve = sqrt(gpvel(1,1)*gpvel(1,1) + gpvel(2,1)*gpvel(2,1))
              percn =0.0d0 ! stabiliztaion term when canopy model
              dracn =0.0d0

              ! if canopy, leaf area index distribution
              if (kfl_canop.and.gpcod.lt.heica) then 
                 if (kfl_candi==0) then ! uniform distribution (default)
                    LADEN =LAD
                 else if (kfl_candi==1) then   ! Canopy distribution from LALIC and MIHAILOVIC, 
                    !http://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%282004%29043%3C0641%3AAERDLD%3E2.0.CO%3B2
                    zm =0.7*heica
                    lam =2.0d0*LAD
                    if (gpcod.lt.zm) then
                       n =6.0
                    else                      
                       n =0.5
                    end if
                    LADEN = lam * (((heica -zm) /max(heica -gpcod,0.001))**n )*exp(n*(1.0d0-(heica -zm) /max(heica -gpcod,0.001)))
                 else if(kfl_candi==2) then       ! Lopes da costa one-dim
                    ladta(1) = 0.43d0
                    ladta(2) = 0.45d0
                    ladta(3) = 0.56d0
                    ladta(4) = 0.74d0
                    ladta(5) = 1.1d0
                    ladta(6) = 1.35d0
                    ladta(7) = 1.48d0
                    ladta(8) = 1.47d0
                    ladta(9) = 1.35d0
                    ladta(10) = 1.01d0
                    ladta(11) = 0.00d0
                    do index= 1,11 ! coordinates of the table
                       codta(index) = dfloat(index-1)/10.0d0*heica
                    end do

                    if (gpcod.gt.heica) then
                       LADEN= 0.0d0
                    else 
                       iz =1
                       jz =11
                       kz = 11/2            
                       do while ((jz-iz).gt.1)                 
                          if (gpcod.lt.codta(kz)) then
                             jz = kz                  
                          else
                             iz = kz
                          end if
                          kz = (iz+jz)/2
                       end do
                       facto = (gpcod- codta(iz))/(codta(jz)-codta(iz))
                       LADEN = ladta(iz) + facto*(ladta(jz)- ladta(iz))
                       LADEN = LADEN*LAI/heica
                    end if              
                 end if
              end if

              ! If prescribed pressure gradient from mesoscale (implemented by JBR)
              if (kfl_pressgr) then
                 ugeos(1)=gpu_geo
                 ugeos(2)=gpv_geo
              end if
              
              force =0.0d0
              react =0.0d0
            ! damping from Palma's paper  
            !  leng_damp = 1000.0
            !  a_damp = 0.0d0
            !  if ( (coord(npoin)- gpcod).lt.leng_damp) then
            !     a_damp= damping*((sin(pi/2.0*(leng_damp - (coord(npoin)- gpcod)/leng_damp)))**2.0)
!                 print *, a_damp
            !  else
            !     a_damp = 0.0d0
            !  end if
              a_damp = damping*min(gpcod/z_damping,1.0d0)
              select case (iunkn)
              case(1) !x- velocity
                 !      Diffusive term
                 diffu = gpmut
                 !      Reactive  term
                 react = 0.0d0
                 react = react + a_damp !(implemented by JBR)
                 !      Rhs term
                 prgra = 0.0
                 ! stationary damping??
                 force =  densi*fcori*(gpvel(2,1)-ugeos(2)) + prgra !+ a_damp*ugeos(1)
                 !nudging = time_weight*spatial_weight*innovation
                 if (kfl_mom_nudging.and.gpcod.gt.500.0d0) then !(implemented by JBR)
                    force = force + densi*(1.0d0/3600.0d0)*(gpcod/coord(npoin))*(gpu_meso-gpvel(1,1)) 
                 end if
                 pert = abs(fcori)*densi
                 if (kfl_canop.and.gpcod.lt.heica) then
                    dracn = cdcan*LADEN*modve  !drag force 
                    react = react +dracn      !perturbation term
                    percn = dracn ! stabilization term (into tau and perturbation functions)
                 end if

                 if (kfl_local) dtinv = (4.0d0*diffu/chale/chale + pert)*safet
                 force = force + dtinv*densi*gpvel(1,2)
                 react = react + dtinv*densi
                 force = force + densi*fcori*gpu_adv !(implemented by JBR)

              case(2) !y- velocity
                 !     Diffusive term
                 diffu = gpmut
                 !     Reactive  term
                 react = 0.0d0
                 react = react + a_damp !(implemented by JBR)
                 !      Rhs term
                 force =  -densi*fcori*(gpvel(1,1)-ugeos(1)) !+ a_damp*ugeos(1)
                 !nudging = time_weight*spatial_weight*innovation
                 if (kfl_mom_nudging.and.gpcod.gt.500.0d0) then !(implemented by JBR)
                    force = force + densi*(1.0d0/3600.0d0)*(gpcod/coord(npoin))*(gpv_meso-gpvel(2,1))
                 end if
                 pert = abs(fcori)*densi
                 if (kfl_canop.and.gpcod.lt.heica) then
                    dracn = cdcan*LADEN*modve
                    react = react + dracn
                    percn = dracn ! stabilization term (into tau and perturbation functions)
                    if (gpcod.lt.1.2*heica) ustar_can = sqrt( gpmut*sqrt(grvel(1)*grvel(1) + grvel(2)*grvel(2)))
                 end if
                 
                 if (kfl_local) dtinv = (4.0d0*diffu/chale/chale + pert)*safet
                 force = force + dtinv*densi*gpvel(2,2)
                 react = react + dtinv*densi
                 force = force + densi*fcori*gpv_adv !(implemented by JBR)
                 
              case(3) ! k-equation
                 if (kfl_logva) then
                    !     Diffusive term
                    diffu = gpmut/sigka

                    rfact = densi*densi*cmu/gpmut*exp(gpkey(1))
                    
                    S = (grvel(1)*grvel(1)+grvel(2)*grvel(2))

                    react = rfact        
                    pert =  rfact 
                    force = rfact*(gpkey(1)-1.0d0)  + densi*cmu*exp(gpkey(1)-gpeps(1))* S + diffu*grkey*grkey

                    if (kfl_local) dtinv = (4.0d0*diffu/chale/chale +pert)*safet

                    !     Temporal integration
                    react = react + dtinv*densi
                    force = force + dtinv*densi*gpkey(2)
                    
                 else
                    !             lm =kar*(gpcod+rough +dwall)*l_max/(kar*(gpcod+rough +dwall)+l_max)
                    !             gpeps(1) = ((cmu*gpkey(3)*gpkey(3))**(0.75d0))/lm
                    !     Diffusive term
                    diffu = gpmut/sigka
                    !     reactive term
                    !              react = 2.0d0*densi*densi*cmu*gpkey(1)/gpmut
                    !              gpeps(1) = densi/lm*(cmu**gpkey(3)*gpkey(3))**(1.5d0)
                    !              react = max(1.0d-14,2.0d0*densi*gpeps(1)/gpkey(3)/gpkey(3)*gpkey(1))

                    !              rfact = densi*densi*Cmu*gpkey(1)/gpmut
                    rfact = densi*gpkey(1)*gpeps(1)/gpkey(3)/gpkey(3)

                    if (rfact.gt.1.0d-14) then
                       react = linke *rfact        
                       pert =  rfact !densi*gpeps(1)/gpkey(3)
                       force = (linke-1.0d0)*densi*gpeps(1) + produ
                    else
                       react = 0.0d0
                       pert =  0.0d0 !  rfact !densi*gpeps(1)/gpkey(3)
                       force = produ
                    end if

                    !              (linke-1.0d0)* rfact*gpkey(1) +produ 
                    !              max(0.0d0,1.0d0*react*gpkey(1)) + produ
                    !              force = densi*gpeps(1) +produ
                    !              force = 0.5*react*gpkey(1) + produ
                    !              force =  produ -  densi*cmu**(0.75)*gpkey(3)**(1.5)/lm
                    !              force = 1.0d0*densi*gpeps(1)/gpkey(3)/gpkey(3)*gpkey(1)*gpkey(1) + produ
                    if (kfl_canop.and.gpcod.lt.heica) then                    
                       force = force + betap*cdcan*LADEN*modve*modve*modve
                       react = react + betad*cdcan*LADEN*modve*isone
                       pert = pert + betad*cdcan*LADEN*modve*isone
                    end if
                    !     Temporal integration
                    react = react + dtinv*densi
                    force = force + dtinv*densi*gpkey(2) + densi*epsam
              end if

              case(4)  ! eps - equation!
                 if (kfl_logva) then ! logarithmic variables
                    S = (grvel(1)*grvel(1)+grvel(2)*grvel(2))

                    ! mixing length, 
                    lm = (cmu**(0.75d0)) *exp(1.5d0*gpkey(1)- gpeps(3)) 
     
                    F= lm/l_max -1.0d0

                    DC=(C2-C1)*(F+1.0d0)   !*lm/l_max                   
                    
                    C1p = C1 + DC 
                    !     Diffusive term
                    diffu = gpmut/sigep
                    !  EVALUATION OF DIFFUSION, REACTION AND FORCE TERMS
                    !  LINEARIZATION OF REACTIVE TERM
                    rfact = densi*C2*exp(gpeps(1)-gpkey(1)) ! linearized term

                    react = rfact        
                    pert =  rfact  
                    force = rfact*(gpeps(1)-1.0d0)  +  densi*C1p*Cmu*exp(gpkey(1)-gpeps(3))*S + greps*greps*diffu
                  

                 else ! not logarithmic variables
                    if (kfl_model==2) then !Realizable model
                       eta =sqrt(grvel(1)*grvel(1)+grvel(2)*grvel(2))*gpkey(3)/gpeps(3)
                       C1p= cmuch*max(0.43d0,eta/(5.0d0/cmuch+eta))  

                       ! Mixing length, 
                       lm = (cmu*gpkey(1)*gpkey(1))**(0.75d0) /gpeps(3) 
                       if (lm.lt.0.0d0) lm=l_max                    
                       C1p = C1p + (C2*sqrt(Cmu0) - C1p)*lm/l_max  ! not cmu???               

                       ! Reactive term
                       react = 2.0d0*densi*C2*gpeps(1)/(gpkey(1)+ sqrt(gpvis*gpeps(3))) !-  densi*C1p*uaste
                       
                       ! Rhs term
                       force = 1.0d0*densi*C2*gpeps(1)/(gpkey(1)+ sqrt(gpvis*gpeps(3)))*gpeps(1) +  densi*C1p*uaste*gpeps(3)
                       pert = 0.0d0

                    else  
                       if (kfl_model ==0.or.kfl_thmod==1.or.kfl_trtem) then   
                          lm = (cmu*gpkey(1)*gpkey(1))**(0.75d0) /gpeps(3) 
                          if (lm.lt.0.0d0) lm=l_max
                          !                    if (lm.lt.1.0d-15.or.lm.gt.l_max)  lm = l_max                    
                          F= lm/l_max -1.0d0
                          !                   F= -(lm/l_max+1.0d0)*((1.0d0-lm/l_max)**3.0d0)
                          !                   F = (lm/l_max)**(10.0d0**(-lm/l_max))-1.0d0
                          DC=(C2-C1)*(F+1.0d0)   !*lm/l_max                   
                          !                   Cg= -(C1 +DC)*0.8
                          difkw= cmu*(-(1.0d0/sigka+1.0d0/sigep)*gpkey(3)/gpeps(3)*grkey*greps + &
                               2.0d0/sigka*grkey*grkey)
                          difkw = 0.0d0
                          produ = prodm
                          if (prodt.lt.0.0) then ! stable
                             alpha = 1.0d0 -lm/l_max
                          else                   !unstable
                             alpha =1.0d0-(1.0d0+(C2-1.0d0)/(C2-C1) )*lm/l_max
                          end if
                          Cg=1.0d0+(C1-C2)*alpha
                       else if (kfl_model == 1) then !RNG model
                          eta =sqrt(grvel(1)*grvel(1)+grvel(2)*grvel(2))*gpkey(1)/gpeps(3)
                          DC= -eta*(1.0d0-eta/4.377d0)/(1.0+0.012d0*eta*eta*eta)                 
                       end if
                       C1p = C1 + DC 

                       !  EVALUATION OF DIFFUSION, REACTION AND FORCE TERMS
                       !  LINEARIZATION OF REACTIVE TERM
                       rfact = densi*C2*gpeps(1)/gpkey(1) ! linearized term
                       if (rfact.gt.1.0d-14) then
                          react = linre *rfact        
                          pert =  rfact  !  densi*gpeps(1)/gpkey(3)
                          force = (linre-1.0d0)*rfact*gpeps(1) +  max(0.0d0,gpeps(1)/gpkey(1)*(C1p*produ+Cg*prodt)) + difkw

                       else
                          react = 0.0d0
                          pert =  0.0d0 
                          force =  max(0.0d0,gpeps(1)/gpkey(1)*(C1p*produ+Cg*prodt)) + difkw
                       end if


                       !                 react = max(0.0d0,2.0d0*densi*C2*gpeps(1)/gpkey(1))
                       !     Rhs term
                       !                 force = 0.5d0*react*gpeps(1) + max(0.0d0,C1p*gpeps(1)/gpkey(1)*produ)
                       !                 force = 1.0d0*densi*C2*gpeps(1)/gpkey(1)*gpeps(1) + C1p*gpeps(1)/gpkey(1)*produ
                       !                 pert = 0.5d0*react !densi*c2/gpkey(1)*gpeps(1)
                    end if
                    if (kfl_canop.and.gpcod.lt.heica) then 
                       C5 = C1 - C2 
                       react = react + C5*betad*cdcan*LADEN*modve
                       pert  = pert  + C5*betad*cdcan*LADEN*modve
                       force = force + C4*betap*cdcan*LADEN*modve*modve*modve*gpeps(3)/gpkey(3)


                       !                    force = force - C5*betad*cdcan*LADEN*modve*gpeps(3)                   

                    end if
                    !if (kfl_canop.and.gpcod.lt.heica) force = force + (C2-C1p)*12.0d0*sqrt(cmu)*cdcan*LAD*modve*gpeps(3)
               
                    !     Diffusive term
                    diffu = gpmut/sigep
                    if(gpkey(1).lt.1.0e-14.or.gpeps(1).lt.1.0e-18) &
                         write(*,'(a,x,10(e12.5,x))'), 'gpcod,produ, gpkey, gpeps,diffu', & 
                         gpcod, produ, gpkey(1), gpeps(1), diffu
                 end if
                 if (kfl_local) dtinv = (4.0d0*diffu/chale/chale +pert )*safet
                 !  Temporal integration
                 react = react + dtinv*densi
                 force = force + dtinv*densi*gpeps(2) + densi*C2*epsam*epsam/keyam
                 
              case(5)  !Temperature - equation
                 !      Diffusive term
                 diffu = gpmut/sigte
                 pert = 0.0d0
                 if (kfl_local) dtinv = (4.0d0*diffu/chale/chale + pert)*safet
                 !      Reactive  term
                 react = dtinv*densi
                 !      Rhs term
                 force = dtinv*densi*gptem(2)
                 force = force + densi*gptem_adv !(implemented by JBR)
                 ! Nudging = time_weight*spatial_weight*innovation
                 if (kfl_tem_nudging.and.gpcod.gt.100.0d0) then !(implemented by JBR)
                    force = force + densi*(1.0d0/3600.0d0)*(gpcod/coord(npoin))*(gpth_meso-gptem(2))  
                 end if
                 !
                 if (kfl_canop.and.gpcod.lt.heica) then ! source or sink term due to canopy
                    force = force - htflx_can*eta_can*LAI/heica*exp(-eta_can*LAI*(heica-gpcod)/heica)
                 end if
              end select

              !
              ! GALERKIN TERMS
              !
              !  Diffusion terms
              do inode =1, nnode
                 do jnode =1, nnode
                    wmatr(inode, jnode)= diffu*cartd(inode)*cartd(jnode)
                 end do
              end do
              ! Stabilization paremeter
              tau = 1.0d0/(4.0d0*diffu/chale/chale +pert +percn)

              ! tau = 0.3d0/densi/fcori        


              if (iunkn.lt.3) then  ! assemble of velocity equations
                 pert= tau*densi*fcori*((-1.0d0)**(iunkn))  ! perturbation test function
                 ! pert term needs to be 0.0 in GABLS3 case type
                 if (kfl_case.eq.3.or.kfl_trtem) pert= 0.0 !(implemented by JBR) damps coriolis oscilations in the free atmosphere  

                 index = 3- iunkn ! points to the perpendicular velocity
                 ! residual in the perpendicular velocity
!                 grmut =0.0
                 resid = densi*fcori*(ugeos(iunkn))*((-1.0d0)**(iunkn)) + densi*dtinv*(gpvel(index,1)-gpvel(index,2))  - gpmut*hevel(index)  + dracn*gpvel(index,1) - grmut*grvel(index)
                 
                 ! Reaction 
                 do inode =1, nnode
                    do jnode =1, nnode
                       wmatr (inode, jnode)= wmatr(inode, jnode)  +  densi*shave(inode)*shave(jnode)*abs(pert*fcori)  + react*shave(inode)*shave(jnode)*(1.0d0-tau*percn)  !*(1.0d0-tau*pert)
                    end do
                 end do

                 ! RHS term
                 do inode =1, nnode
                    wrhsi(inode) = force*shave(inode)*(1.0d0-tau*percn)  + pert*resid*shave(inode)
                 end do

              else   ! for turbulence and temper equations               
                 ! pert =0.0d0
                 ! Reaction 
                 do inode =1, nnode
                    do jnode =1, nnode
                       wmatr (inode, jnode)= wmatr(inode, jnode)  +  react*shave(inode)*shave(jnode)*(1.0d0-tau*pert)
                    end do
                    ! RHS term
                    wrhsi(inode) = force*shave(inode)*(1.0d0-tau*pert) 
                 end do

              end if

              if (abs(ctime - 3600.0d0*35.0d0).lt.0.001.and.iunkn.eq.5 ) then ! calculates divergence of heat flux
                 do inode =1, nnode
                    gpdhf(inode)  =  gpdhf(inode) + diffu*grtem*cartd(inode)*dvolu
                 end do
              end if

              ! Add Gauss point contributions
              elrhs=elrhs+dvolu*wrhsi
              elmat=elmat+dvolu*wmatr
 
           end do gauss_points

           ! Matrix Assembly and RHS     
           !        write(*,*) 'assembly'
!!$           if (abs(ctime - 3600.0d0*35.0d0).lt.0.001.and.iunkn.eq.5 ) then ! calculates divergence of heat flux
!!$              do inode =1, nnode ! heat flux divergence
!!$                 ipoin = lnods(inode, ielem)
!!$                 divhf(ipoin)= divhf(ipoin) + gpdhf(inode)  
!!$              end do
!!$           end if
           do inode =1, nnode
              ipoin = lnods(inode, ielem)
              do jnode =1, nnode   
                 jpoin = lnods(jnode, ielem)
                 izdom = ia(ipoin)
                 do while (ja(izdom).ne.jpoin) 
                    izdom= izdom +1
                 end do
                 amatr(izdom) = amatr(izdom) + elmat(inode, jnode)
              end do
              rhsid(ipoin)= rhsid(ipoin) + elrhs(inode)                           
           end do

        end do elements
        
        if (iunkn.eq.5) then
           do ipoin =1, npoin
              rhsid(ipoin) = rhsid(ipoin) +  divhf(ipoin)
           end do
        end if
        
        ! Boundary assembly (only if wall law)
        if(kfl_bouco_vel.gt.0.and.(iunkn.lt.3.or.iunkn.eq.5)) then
           if (kfl_bouco_vel.eq.1.and.iunkn.lt.3) then  !ustar is fix (only valid for playing)
              vnorm= sqrt(veloc(1,1,1)*veloc(1,1,1)+veloc(1,2,1)*veloc(1,2,1))
              tract = densi*ustar*ustar*veloc(1,iunkn,1)/vnorm
              ustar2= ustar
              B = densi*ustar*ustar2/vnorm/vnorm
              ui = veloc(1,iunkn,1)
              dtrac =   B*(vnorm+ui*ui/vnorm)
              tract =   B*ui*ui*ui/vnorm
           else if(kfl_bouco_vel.eq.2) then ! wall b.c.
              !
              !*** IF THERMAL MODEL CALCULATES MONIN OBUKHOV CORRECTIONS FOR WALL MODEL
              !
              if (kfl_thmod.gt.0) then  ! thermal coupling 
                 eta = dwall/lmoni  ! dimensionless height parameter
                 if (lmoni.lt.-1.0d-6.and..not.kfl_canop) then     ! unstable 
                    phi_m = (1.0 -mo_gamm1*eta)**(-0.25)
                    phi_h = mo_prand*(1.0 -mo_gamm2*eta)**(-0.5)
                    logfu = log(1.0d0+dwall/rough)+log((8.0d0*phi_m**4.0)/((phi_m+1.0)*(phi_m+1.0)*(phi_m*phi_m+1.0)) &
                         - 0.5d0*pi + 2.0*atan(1.0d0/phi_m))
                    logft = log(1.0+ dwall/rough) -2.0*log(0.5*(1.0+ mo_prand/phi_h))
                    phi_e = 1.0d0 - eta
                 else if (lmoni.gt.1.0d-6.and..not.kfl_canop) then  ! stable
                    phi_m = 1.0 + mo_beta1*eta
                    phi_h = mo_prand +mo_beta2*eta
                    logfu = log(1.0+ dwall/rough) + phi_m -1.0d0 
                    logft = log(1.0+ dwall/rough)  + phi_h/mo_prand -1.0d0
                    phi_e = phi_m - eta                    
                 else ! neutral
                    logft =  log(1.0+ dwall/rough)
                    logfu =  log(1.0+ dwall/rough)
                    phi_e= 1.0
                 end if
              else ! not thermal model
                 logfu=log(1.0+ dwall/rough)
              end if
              if (iunkn.lt.3) then
                 vnorm= sqrt(veloc(1,1,1)*veloc(1,1,1)+veloc(1,2,1)*veloc(1,2,1))
                 if (.not.kfl_abl2) then ! ustar only in terms of velocity
                    ustar = kar*vnorm/logfu
                    ustar2= ustar
                    B = densi*ustar*ustar2/vnorm/vnorm
                    ui = veloc(1,iunkn,1)
                    dtrac =   B*(vnorm+ui*ui/vnorm)
                    tract =   B*ui*ui*ui/vnorm
                 else  !ABL2
                    ustar = kar*vnorm/logfu
                    ustar2 = cmu0**0.25*sqrt(keyva(1, 1))*(phi_m/phi_e)**0.25
                    if (kfl_logva)    ustar2 = cmu0**0.25*exp(0.5d0*keyva(1, 1))
                      !       dtrac = densi*ustar*ustar2/vnorm
                    dtrac = densi*ustar*ustar2/vnorm !*0.5d0
                    tract = 0.0d0 !- dtrac* veloc(1,iunkn,2)
                 end if
              else if (iunkn.eq.5) then ! wall law for temper
                 ustar2 = cmu0**0.25*sqrt(keyva(1, 1))*(phi_m/phi_e)**0.25
                 if (.not.kfl_abl2) ustar2=ustar              
                 hconv = densi*ustar2*kar/(logft*mo_prand)              
                 tract = hconv*tewal ! right h side coeff
                 dtrac = hconv       ! matrix assembly
                 ! for adiabatic uncoment
!                 tract = 0.0 ! htflx_can*exp(-0.6*6.0*(heica-coord(1))/heica) !   0.0 ! -60.0/1000
!                 dtrac = 0.0
              end if
           end if
           !  non linearized
           !        rhsid(1) = rhsid(1) -tract
           ! linearized
           !       tract = rho*ustar*ustar*u/un        
           !       ustar*ustar = A*un*un 
           !       tract = B*un*u => B=A*rho=rho*ustar*ustar/(un*un)
           !       dtrac = B(un*un+u*u)/un
           !       tract =  tract(ui) + dtrac*(u-ui)
           !       tract \approx  + dtrac*u -B*ui*ui*ui/un
           !       tract \approx    dtrat*u -rho*ustar*ustar*(ui/un)**3
           
           rhsid(1) = rhsid(1) + tract
           amatr(1) = amatr(1) + dtrac
        end if
        
           !
           ! Prescribe Dirichlet boundary conditions over the matrix
           !
        wvalu = -10000.0
        tvalu = -10000.0

        !DIRICHLET

        if(iunkn.lt.3) then        ! wall and top velocity values
           if(kfl_bouco_vel.eq.0)  wvalu = 0.0d0
           if (abs(fcori).lt.1.0d-9.or.kfl_topco_vel.eq.1) &
                tvalu = ugeos(iunkn) ! prescribe if not geostrophic pressure (implemented by JBR)
           ! If canopy  symmetry b.c. on top
!           if (kfl_canop) tvalu = -10000 
        elseif (iunkn.eq.3) then   ! wall and top kinet values
           if (.not.kfl_abl2) then
              wvalu = ustar*ustar/sqrt(cmu0)
              if (kfl_logva)  wvalu = log(wvalu)
           end if
!           if (abs(fcori).lt.1.0d-9) tvalu = wvalu
         ! impose tke on bound top when imposing velocity
           if (kfl_topco_vel.eq.1) & 
                tvalu = keyam
        elseif (iunkn.eq.4) then   ! wall and top epTsilon values
           if (kfl_thcou) then
              eta = dwall/lmoni
              if (lmoni.lt.-1.0d-6) then !unstable
                 phi_m   = (1.0 -mo_gamm1*eta)**(-0.25)
                 phi_e = 1.0d0 - eta
              else if (lmoni.gt.1.0d-6) then ! stable
                 phi_m   = 1.0 + mo_beta1*eta
                 phi_e = phi_m- eta
              end if
           end if
           if (.not.kfl_abl2) then
              wvalu = ustar*ustar*ustar/(kar*(dwall+rough))*phi_e
              if (kfl_logva)    wvalu = log(wvalu)
              !           if (abs(fcori).lt.1.0d-9) tvalu = ustar*ustar*ustar/(kar*(dwall+rough+ length))
              !        wvalu = cmu0**0.75*abs(keyva(1,1))**1.5/(kar*(dwall+rough))
              !        print *, 'rough', rough, dwall
              !        wvalu = ustar*ustar*ustar*(l_max+kar*rough)/(kar*rough*l_max)
              !        tvalu = cmu**0.75*((vegeo*vegeo*1.0e-6)**(1.50d0))/l_max
           else  !ABL2
              ustar2= cmu0**0.25*sqrt(keyva(1, 1))*(phi_m/phi_e)**0.25
              wvalu = ustar2*ustar2*ustar2/(kar*(dwall+rough))*phi_e

              if (kfl_logva)    then
                 ustar2 = cmu0**0.25*exp(keyva(1, 1)*0.5d0)
                 wvalu = log(ustar2*ustar2*ustar2/(kar*(dwall+rough)))
              end if
           endif
           if (kfl_topco_vel.eq.1) &
                tvalu = epsam
        else if (iunkn.eq.5) then   !wall temper value
           if (kfl_bouco_vel.eq.0)  wvalu = tewal
!           if (kfl_case.eq.3) then
!              wvalu = tewal
!           end if           
!           tvalu = tetop
        end if

        inico= 2
        finco= 1
        if (wvalu.gt.-1.0) then 
           inico = 1
           finco = 1
           bvalu(1)= wvalu   
        end if
        if (tvalu.gt.-1.0) then
           finco = 2 
           bvalu(2)= tvalu
        end if
    
        do inode = inico, finco
           ipoin = 1                    ! bottom dirichlet cond
           if (inode==2) ipoin =npoin   ! top dirichlet cond
           do izdom = ia(ipoin), ia(ipoin+1) -1
              jpoin = ja(izdom)
              if (jpoin.eq.ipoin) then
                 adiag = amatr(izdom)
              else
                 amatr(izdom) =0.0
              end if
           end do
           rhsid(ipoin) = bvalu(inode)*adiag
        end do
        !  write(*,*) 'going to solve'
        !
        ! Solve the algebraic system.
        !
        !  write (*,*) 'going to solve...'
        ! converts csr matrix to tri diagonal matrix
        avect(1) = 0.0d0    ! subdiagonal
        diago(1) = amatr(1) ! diagonal
        cvect(1) = amatr(2) ! supdiagonal
        do ipoin = 2, npoin -1
           avect(ipoin) = amatr (ia (ipoin))
           diago(ipoin) = amatr (ia (ipoin) +1)
           cvect(ipoin) = amatr (ia (ipoin) +2)      
        end do
        avect(npoin) = amatr (ia (npoin))
        diago(npoin) = amatr (ia (npoin) +1)
        cvect(npoin) = 0.0d0
        !       write(*,*) 'calling solver'
        call solvtr(avect, diago, cvect,rhsid, unkno, npoin)
        !          call solvpa(ia, ja, amatr, rhsid, unkno, npoin, 1, &
        !                       1, .false., 1)


        !
        !   Check convergence and updates unknown
        !
        tnume = 0.0d0
        tdeno = 0.0d0
        select case (iunkn)
        case (1) !ux velocity
           do ipoin =1, npoin
              tnume= tnume + (veloc(ipoin, 1, 1) - unkno(ipoin))*(veloc(ipoin, 1, 1) - unkno(ipoin))
              tdeno= tdeno + unkno(ipoin)*unkno(ipoin)
              veloc(ipoin, 1, 1) = unkno(ipoin)
           end do
        case (2) ! uy velocity
           do ipoin =1, npoin
              tnume= tnume + (veloc(ipoin, 2, 1) - unkno(ipoin))*(veloc(ipoin, 2, 1) - unkno(ipoin))
              tdeno= tdeno + unkno(ipoin)*unkno(ipoin)
              veloc(ipoin, 2, 1) = unkno(ipoin)
           end do
        case (3)
           do ipoin =1, npoin
              tnume= tnume + (keyva(ipoin, 1) - unkno(ipoin))*(keyva(ipoin, 1) - unkno(ipoin))
              tdeno= tdeno + unkno(ipoin)*unkno(ipoin)
              keyva(ipoin, 1) = (0.7d0*unkno(ipoin)+0.3d0*keyva(ipoin,1))
              !           lm =kar*(coord(ipoin)+rough)*l_max/(kar*(coord(ipoin)+rough)+l_max)
              !           epsil(ipoin,1) =  ((cmu*keyva(ipoin,1)*keyva(ipoin,1))**(0.75d0))/lm
              if(kfl_thmod.eq.1)  keyva(ipoin, 1)= max(keyva(ipoin, 1), keyam) !.and.coord(ipoin).gt.10.0)
           end do
        case (4)
           do ipoin =1, npoin
              tnume= tnume + (epsil(ipoin, 1) - unkno(ipoin))*(epsil(ipoin, 1) - unkno(ipoin))
              tdeno= tdeno + unkno(ipoin)*unkno(ipoin)
              epsil(ipoin, 1) =(0.7d0*unkno(ipoin)+0.3d0*epsil(ipoin,1))
              if(kfl_thmod.eq.1) then ! limits to k ans eps
                 epsil(ipoin, 1)= max(epsil(ipoin, 1), epsam)             
              end if
           end do
        case (5)    ! transient temper equation
           do ipoin =1, npoin
              tnume= tnume + (tempe(ipoin, 1) - unkno(ipoin))*(tempe(ipoin, 1) - unkno(ipoin))
              tdeno= tdeno + unkno(ipoin)*unkno(ipoin)
              tempe(ipoin, 1) =unkno(ipoin)
           end do
        end select

        if (tdeno.gt.1.0d-14) then
           error = sqrt(tnume/tdeno)*100.0d0
        else 
           error =0.0d0
        end if
        write (lun_conve(iunkn), 101) istep, iiter, itera, ctime, error, sqrt(tdeno/npoin)

        ! Update inner iteration counter
        itera = itera +1 
     end do !end convergence 

     if (error.lt.toler(iunkn)) kfl_goite=0 ! do not continue iterating
     ! relaxation for epsil
     if (iunkn.eq.4) epsil(1:npoin,1)=0.5d0*(epsil_aux(1:npoin)+epsil(1:npoin,1))

     ! actualization of ustar value when using Dirichlet cond. for velocity
     if (iunkn.eq.2.and.kfl_bouco_vel.eq.0) then
        gpmut = cmu*keyva(1,1)*keyva(1,1)/epsil(1,1)
        ipoin =2
        vnorm= sqrt(veloc(ipoin,1,1)*veloc(ipoin,1,1)+veloc(ipoin,2,1)*veloc(ipoin,2,1))
        ustar= sqrt(gpmut*vnorm/coord(ipoin))
        !     write (*,*) 'ustar=', usta
     end if
     if (.not.kfl_canop) ustar_can = ustar
     
     if (kfl_canop.and.iunkn.eq.1.and.flctr) then ! control of mean velocity throug pressure gradient
        !  perform the integration of velocity
        umean (1) = 0.0d0
        do ielem = 1, nelem 
           umean (1) = umean (1) + 0.5d0*(veloc(ielem+1,1,1)+veloc(ielem, 1, 1))*(coord(ielem+1)-coord(ielem))
        end do
        umean (1) = umean (1)/ (coord(npoin) -coord(1))
        if (istep ==1) umean(2) = umean(1)
        prgra = prgra + dtinv * ( umtar -  2.0d0*umean(1) + umean(2) )
        print *, 'pr grad', prgra, 'umean', umean(1)
!        prgra = prgra + dtinv* ( umtar - umean(1))
        ! actualizes u before ...
        umean (2)  = umean(1)
     end if
     ! calculates Mellor-Yosida length for transient thermal flows, after TKE loops
     if (iunkn.eq.3.and.(kfl_thmod.eq.1)) then
        tnume = 0.0d0
        tdeno = 0.0d0
        do ielem =1, nelem
           ipoin = ielem
           jpoin = ielem + 1
           chale = coord(jpoin)-coord(ipoin)
           sqkey = 0.50d0*(sqrt(abs(keyva(ipoin,1)))+sqrt(abs(keyva(jpoin,1))))
           gpcod = coord(ipoin) +0.50d0*chale
           tnume = tnume + gpcod*sqkey*chale        
           tdeno = tdeno + sqkey*chale     
        end do
        !     Mellor-Yosida length
        lenmy = 0.075*tnume/tdeno
     end if
     cmu = cmu0
     !
     ! Formats
     !
101  format(4x,i9,2x,i9,2x,i9,5(2x,e14.6))

   end subroutine solite
