   subroutine Postpr
  ! This routines writes 
  ! 1) The obtained time step profiles to a file 
  ! 2) Global variables after each time step
     ! 3) Point evolution variables after each time step
     ! And results to a given lenght
     ! set file names 
     ! profiles
     ! cuts
     ! global variables
     
     use def_master
     implicit none
     integer :: kpoin
     real(8) :: grvel(2)
     real(8) :: tstar,  gpnut, facto, lenmo, htpbl, stre0
     real(8) :: velcu(2,mcuts), keycu(mcuts), epscu(mcuts), temcu(mcuts)
     real(8) :: cmu_p! cmu postrprocess
     real(8) :: vel,mut,mixle,Richf,prodm,prodt,htflx,stres, grtem
     character*20 :: label, file
     !
     ! Write global evolution variables
     !
     gpnut =  cmu*keyva(1,1)*keyva(1,1)/epsil(1,1)
     ! qwall positive in stable atm, when floor is cold (heat flux going out)
     hflx0 =  rhocp*gpnut*(tempe(2,1)- tempe(1,1))/(coord(2)-coord(1))/sigte
     tstar =  hflx0/(rhocp*ustar_can)
     if (abs(tstar).lt.1.0e-8) tstar = 1.0e-8
     lenmo = ustar_can*ustar_can*teref/(kar*gravi*tstar)
     stre0 = ustar_can*ustar_can*densi
! looks for Height top boundary layer
     do ipoin = 2, npoin
        gpnut = densi*cmu*keyva(ipoin,1)*keyva(ipoin,1)/epsil(ipoin,1)       
        grvel(1) =(veloc(ipoin,1,1)-veloc(ipoin-1,1,1)) /(coord(ipoin)-coord(ipoin-1))
        grvel(2) =(veloc(ipoin,2,1)-veloc(ipoin-1,2,1)) /(coord(ipoin)-coord(ipoin-1))              
        stres = gpnut*sqrt(grvel(1)*grvel(1)+grvel(2)*grvel(2))
        kpoin = ipoin
        if (stres.lt.0.05*stre0.and.(coord(ipoin).gt.1.2*heica)) exit
     end do
     ! height top bound layer
     htpbl =coord(kpoin)
     if (istep.eq.0) then
        write(lun_globa,'(a,16a17)')  '#', '1:time', '2:ustar', '3:Htflx', '4:tstar', '5:tewal', '6:L_MY', '7:L_MO', '8:H_PBL', '9:uxgeo', '10:uygeo(2)', '11:mod_ugeo', '12: uxtop', '13:uy_top', '14:mod_utop'
     else
        write (lun_globa,'(20(e17.10,2x))')  ctime , ustar_can, hflx0, tstar, tewal, lenmy, lenmo, htpbl,ugeos(1),ugeos(2),sqrt(ugeos(1)*ugeos(1)+ugeos(2)*ugeos(2)),veloc(npoin,1,1),veloc(npoin,2,1),sqrt(veloc(npoin,1,1)*veloc(npoin,1,1) + veloc(npoin,2,1)*veloc(npoin,2,1))
     end if
          !
     ! Write cut evolution
     !
     do icuts =1, ncuts      
        ! values at coord icuts
        ipoin = ielec(icuts)
        jpoin = ipoin +1
        facto = (cutpr(icuts)-coord(ipoin))/(coord(jpoin)-coord(ipoin))
        ! interpolate values
        velcu(1:2,icuts) = veloc(ipoin,1:2,1) + facto*(veloc(jpoin,1:2,1)-veloc(ipoin,1:2,1))    
        keycu(icuts) = keyva(ipoin,1) + facto*(keyva(jpoin,1)-keyva(ipoin,1))
        epscu(icuts) = epsil(ipoin,1) + facto*(epsil(jpoin,1)-epsil(ipoin,1))
        temcu(icuts) = tempe(ipoin,1) + facto*(tempe(jpoin,1)-tempe(ipoin,1))
        gpnut = cmu*keycu(icuts)*keycu(icuts)/epscu(icuts)
        grtem = (tempe(jpoin,1)-tempe(ipoin,1))/(coord(jpoin)-coord(ipoin))
        htflx = -rhocp*gpnut/sigte*grtem
   
        write(lun_cutre(icuts),'(10(e17.10,2x))') ctime, velcu(1,icuts), velcu(2,icuts), &
             keycu(icuts), epscu(icuts), temcu(icuts), htflx
     end do

     !
     !*** Writes inital profiles
     !
     if (istep.eq.0) then
        file = 'Profiles.00'
        open(lun_ini, FILE=file, status='unknown' )
        write(lun_ini,'(a,6a17)')  '#', '1:coord', '2:U', '3:V', '4:key', '5:eps', '6:temp'
        do ipoin = 1, npoin
           write (lun_ini, '(6(e17.10, x))') coord(ipoin),veloc(ipoin,1,1),veloc(ipoin,2,1), &
                keyva(ipoin,1), epsil(ipoin,1), tempe(ipoin,1)
        end do
        close (lun_ini)
     end if
     
     !
     ! Write profiles at some time steps
     !
     if (mod(istep, stepr)==0) then
        print *,'Writing profiles at ',ctime/3600,'hours'
        write(label,'(f10.2)') ctime/3600.0 - mod(ctime/3600,1.0)*0.4 ! time in hours         
        file = 'plotprofile_'//trim(adjustl(label))  ! plotprof_hh.hh
        open(lun_timre, FILE=file,status='unknown')
        write(lun_timre,'(a)') '#1: coord          2:velocx          3:velocy          &
             4: modvel          5:key            6:eps             7: mutur         &
             8: mixlen         9:temp            10:rich_flux      11: prodm        &
12: prodt         13:htflx          14:stress         15: cmu    16:u/u_star   17:heat/ustar'
        
        cmu_p = cmu ! cmu for postprocess (needed for realizable model)
        do ipoin = 1, npoin
           !
           vel = 0.0d0
           call Postpr_ipoin(vel,mut,mixle,Richf,prodm,prodt,htflx,stres)
           !
           write (lun_timre, '(21(e17.10, x))') coord(ipoin),veloc(ipoin,1,1),veloc(ipoin,2,1), &
                vel, keyva(ipoin,1), epsil(ipoin,1), mut, mixle, tempe(ipoin,1), &
                Richf, prodm, prodt,htflx, stres, cmu_p, vel/ustar_can, htflx/ustar_can
           !
        end do
        close (lun_timre)
     end if

     return
   end subroutine Postpr

   
   subroutine Postpr_ipoin(vel,gpnut,mixle,Richf,prodm,prodt,htflx,stres)
     !This subroutine calculates the output variables at each node
     
     use def_master
     implicit none
     real(8),intent(out) :: vel, gpnut,mixle,Richf,prodm,prodt,htflx,stres
     real(8)             :: grvel(2), grtem, deltz,cmu_p
     real(8)             :: W, uaste, A0, As

     cmu_p = cmu
     ! realizable model
     if (kfl_model ==2)   A0=(1.0d0-3.0d0*sqrt(0.5d0*cmu0))/cmu0
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
        W=0.0d0
        As=1.5d0*sqrt(2.0d0)
        cmu_p = 1.0d0/(A0+As*uaste*keyva(ipoin,1)/epsil(ipoin,1))
     end if
     gpnut = cmu_p*keyva(ipoin,1)*keyva(ipoin,1)/epsil(ipoin,1)         
     prodm  =densi*gpnut*(grvel(1)*grvel(1)+grvel(2)*grvel(2)) 
     sigte =0.74d0
     if (grtem.lt.0.0d0) then ! unstable
        ! richf  = gravi/teref*grtem/(sigte*(grvel(1)*grvel(1)+grvel(2)*grvel(2)))    !JBR
        richf  = gravi/teref*grtem/(grvel(1)*grvel(1)+grvel(2)*grvel(2))
        if(.not.kfl_canop)  then
           sigte = 0.74d0*(1.0d0-15.d0*richf)**(-0.25d0)
           ! richf  = richf*0.74d0/sigte
           ! sigte = 0.74d0*(1.0d0-15.d0*richf)**(-0.25d0)
        end if
     end if
     ! thermal production
     prodt  = -densi*gravi/teref*gpnut*grtem/sigte
     ! richardson flux
     if (prodm.gt.1.0e-10 ) then
        richf  = prodt/prodm
     else 
        richf =1e6
     end if
     htflx = -rhocp*gpnut/sigte*grtem
     stres = densi*gpnut*sqrt(grvel(1)*grvel(1)+grvel(2)*grvel(2))
     mixle  = ((cmu_p*keyva(ipoin,1)*keyva(ipoin,1))**(0.75d0))/epsil(ipoin,1)
     vel   = sqrt(veloc(ipoin,1,1)*veloc(ipoin,1,1)+ veloc(ipoin,2,1)*veloc(ipoin,2,1))

   end subroutine Postpr_ipoin
