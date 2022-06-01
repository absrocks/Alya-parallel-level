subroutine Turnof
  ! This routines writes the obtained results to a file 

  use def_master
  implicit none
  real(8) :: viscot, tempa, tempz, richf, grvel(2), prodt, prodm, phi_e, phi_m, thvel, thtem
  real(8) :: thvte, thmut, zbylm,tstar,  theps, gpnut, mixle, stres, deltz, grtem, htflx
  real(8) :: WW, uaste, A0, As, cmu_p
  integer :: ipoin

  richf = 0.0d0
  tempa = 0.0d0
  tempz = tempa
  print *, 'lmoni', lmoni
  open(lun_postp, FILE='plotresults',status='unknown')
  write(lun_postp,'(a)') '#1: coord          2:velocx          3:velocy          4: modvel         5:key             6:eps             7: mutur          8: mixlen         9:temp            10: stres         11:heatflx       12:rich_flux      13: prodm         14: prodt         15: theor vel   16: theor.temper  17: theor.virt temp 18: th.eps     19:th. mut'
  if (kfl_logva) then
     keyva(1:npoin,1) = exp(keyva(1:npoin,1))
     epsil(1:npoin,1) = exp(epsil(1:npoin,1))
  end if
  if (kfl_thmod.eq.1) then ! evaluate MO-length and wall heat flux
     gpnut =  cmu*keyva(1,1)*keyva(1,1)/epsil(1,1)
     ! qwall positive in stable atm, when floor is cold (heat flux going out)
     hflx0 =  rhocp*gpnut*(tempe(2,1)- tempe(1,1))/(coord(2)-coord(1))/sigte
     lmoni =  rhocp*teref*ustar*ustar*ustar/(kar*gravi*hflx0)
  end if
  ! realizable model
  if (kfl_model ==2)   A0=(1.0d0-3.0d0*sqrt(0.5d0*cmu0))/cmu0
  do ipoin = 1, npoin
     viscot = densi*cmu*keyva(ipoin,1)*keyva(ipoin,1)/epsil(ipoin,1)

     ! temper distribution          

     if (ipoin.eq.1.and.kfl_thmod.ne.0) tempz =  hflx0*sigte/rhocp/viscot*coord(1)
     if (ipoin.gt.1.and.kfl_thmod.ne.0.and.coord(ipoin).lt.ztsbl)  tempz = tempa + &
          hflx0*sigte/rhocp/viscot*(1.0d0-coord(ipoin)/ztsbl)*(coord(ipoin)-coord(ipoin-1))
     tempa = tempz
     if (kfl_thmod.eq.1) tempz = tempe(ipoin,1) 
     ! production distribution
     prodm =0.0d0
     prodt =0.0d0
     if (ipoin.lt.npoin) then
        deltz = coord(ipoin+1)-coord(ipoin)
        grvel(1) =(veloc(ipoin+1,1,1)-veloc(ipoin,1,1)) /deltz
        grvel(2) =(veloc(ipoin+1,2,1)-veloc(ipoin,2,1)) /deltz
        grtem =(tempe(ipoin+1,1)-tempe(ipoin,1)) /deltz
     else
        deltz = coord(npoin)-coord(npoin-1)
        grvel(1) =(veloc(npoin,1,1)-veloc(npoin-1,1,1)) /deltz
        grvel(2) =(veloc(npoin,2,1)-veloc(npoin-1,2,1)) /deltz
        grtem =(tempe(npoin,1)-tempe(npoin-1,1)) /deltz
     end if
     if (kfl_model==2) then ! realizable
        uaste = sqrt(grvel(1)*grvel(1)+grvel(2)*grvel(2))
        WW=0.0d0
        As=1.5d0*sqrt(2.0d0)
        cmu = 1.0d0/(A0+As*uaste*keyva(ipoin,1)/epsil(ipoin,1))
     end if
     viscot = densi*cmu*keyva(ipoin,1)*keyva(ipoin,1)/epsil(ipoin,1)
     prodm  =viscot*(grvel(1)*grvel(1)+grvel(2)*grvel(2)) 
     stres = viscot*sqrt(grvel(1)*grvel(1)+grvel(2)*grvel(2))
     htflx = rhocp*viscot*grtem/sigte
     !         if (ipoin.gt.1) then
     !            grvel(1) =(veloc(ipoin,1,1)-veloc(ipoin-1,1,1)) /(coord(ipoin)-coord(ipoin-1))
     !            grvel(2) =(veloc(ipoin,2,1)-veloc(ipoin-1,2,1)) /(coord(ipoin)-coord(ipoin-1))           
     !         end if
     !         thermal production
     richf=0.0d0
     if (kfl_thmod.eq.1) then
        prodt  = -gravi/teref*viscot*grtem/sigte
     end if
     ! richardson flux

     if (kfl_thcou.and.prodm.gt.1.0e-8) then
        Richf  = prodt/prodm
     else
        Richf  = 0.0
     end if

     ! Monin obukhov, theoretical solution (without coriolis)
     if (lmoni.gt.1.0d-4) then  ! stable
        zbyLm = (coord(ipoin))/lmoni
        phi_m = 1.0d0+ beta_mo*zbyLm
        phi_e = phi_m -(coord(ipoin)+rough)/lmoni
        tstar = hflx0/(rhocp*ustar)
        !            theoretical velocity
        thvel = ustar/kar*(log(1.0d0 + coord(ipoin)/rough) + phi_m -1.0d0)
        !            theoretical virtual temperature variation respect to the ground            
        thvte = tstar/kar*(log(1.0d0 + coord(ipoin)/rough) + phi_m -1.0d0) 
        !            theoretical physical temperature
        thtem = thvte  - 9.75d-3*coord(ipoin)

     else if (lmoni.lt.-1.0d-4) then ! unstable
        zbyLm = (coord(ipoin)+rough)/lmoni
        phi_m = (1.0d0- alpha_mo*zbyLm)**expon_mo
        phi_e = 1.0d0 -(coord(ipoin)+rough)/lmoni
        tstar = hflx0/(rhocp*ustar) !OJO, ver signos
        !            theoretical velocity
        thvel = ustar/kar*(log(1.0d0 + coord(ipoin)/rough) + log(8.0d0*phi_m**4.0d0/((phi_m+1.0d0)*(phi_m+1.0d0)*(phi_m*phi_m+1.0d0)) ) -pi/2.0d0 +2.0d0*atan(1.0d0/phi_m)) 
        !            theoretical virtual temperature variation respect to the ground            
        thvte = tstar/kar*(log(1.0d0 + coord(ipoin)/rough)-2.0d0*log(0.5d0*(1.0d0+1.0d0/(phi_m*phi_m)))) 
        !            theoretical physical temperature
        thtem = thvte  - 9.75d-3*coord(ipoin)
     else ! neutral
        phi_m = 1.0d0
        phi_e= 1.0d0
        thvel = ustar/kar*(log(1.0d0 + coord(ipoin)/rough))
        thvte = 0.0d0
        thtem = thvte  - 9.75d-3*coord(ipoin)          

     end if
     !        theroretical epsilon 
     theps = ustar*ustar*ustar/(kar*(coord(ipoin)+rough))*phi_e
     !        theroretical turbulence viscosity
     thmut = densi*kar*ustar*(coord(ipoin)+rough)/phi_m

     mixle  = ((cmu*keyva(ipoin,1)*keyva(ipoin,1))**(0.75d0))/epsil(ipoin,1)

     write (lun_postp, '(21(e17.10, x))') coord(ipoin),veloc(ipoin,1,1),veloc(ipoin,2,1), &
          sqrt(veloc(ipoin,1,1)*veloc(ipoin,1,1)+ veloc(ipoin,2,1)*veloc(ipoin,2,1)),  &
          keyva(ipoin,1), epsil(ipoin,1), viscot, &
          mixle, tempz ,stres, htflx, &
          Richf, prodm, prodt, thvel, thtem, thvte, theps, thmut

     if (kfl_canop.and.coord(ipoin).lt.heica) ustar_can = sqrt(stres)
  end do

  write (*,*) 'ustar=', ustar
  if (kfl_canop) write (*,*) 'ustar_can=', ustar_can


  ! Close files and deallocate structures 
  close(lun_postp)
  close(lun_conve(1))
  close(lun_conve(2))
  close(lun_conve(3))
  close(lun_conve(4))
  if (kfl_thmod.eq.1) then ! For transient problems
     close(lun_conve(5)) ! tem.cvg
     close(lun_globa)    ! resglobal
     do icuts =1, ncuts
        close(lun_cutre(icuts)) 
     end do
#ifdef _NETCDF_      ! This allows to compile without NETCDF library
     if (kfl_netCDF) then
        call netCDF_deallocate
     end if
#endif
  end if
  !
  ! Deallocte structures
  !
  deallocate (lnods)
  deallocate (unkno)
  deallocate (rhsid)
  deallocate (coord)
  deallocate (avect)
  deallocate (diago)
  deallocate (cvect)
  deallocate (veloc)
  deallocate (keyva)
  deallocate (epsil)
  deallocate (shape)
!  deallocate (amatr)
  deallocate (ja)
  deallocate (ia)

end subroutine Turnof
