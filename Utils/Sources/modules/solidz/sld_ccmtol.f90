subroutine sld_ccmtol
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_ccmtol
  ! NAME 
  !    sld_ccmtol
  ! DESCRIPTION

  ! USES

  ! USED BY
  !    sld_bouope
  !***
  !
  !volst(1) -> V_old
  !volst(2) -> V_new
  !volst(3) -> deltaV_old
  !volst(4) -> deltaV_new
  !-----------------------------
  !volst(5) -> V start diastole 
  !volst(6) -> V end   diastole 
  !volst(7) -> V start ejection 
  !volst(8) -> V end   ejection 
  !volst(9) -> V end   iso. relax 
  !volst(10)-> V end   rapid filling 
  !volst(11)->    
  !
  !=============================
  !
  !dltap_sld(1) -> deltaP_old
  !dltap_sld(2) -> deltaP_new
  !dltap_sld(3) -> P_old  (and unsmooth P)
  !dltap_sld(4) -> P_new
  !
  !dltap_sld(11) -> P smooth  
  !dltap_sld(12) -> dP/dt during iso relax - will be used for rapid filling  
  !


  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_solidz
  implicit none


  real(rp)    :: baloc(ndime,ndime)
  real(rp)    :: tract(ndime),ptota,bidon,big,small,reverse,init,sidei
  real(rp)    :: rcste,ccste,rcst2,dfldt,cinve,cporr,xfunc,flowq,pfins,volin,volou,volim
  real(rp)    :: voluf,pfinf,tfinf,pfinr,mfilr,pfirf,alpha,paort
  integer(ip) :: ivolu



  tfinf = parcc_sld(1,1)  !time end fillig (For the first cycle)
  pfinf = parcc_sld(2,1)  !Pressure end filling
  pfins = parcc_sld(3,1)  !Pressure at end diastole (end of isovolumetric contraction)
  pfinr = parcc_sld(4,1)  !Pressure end of isovolum relax
  pfirf = parcc_sld(5,1)  !Pressure end rapid filling initial to zero for the firts cycle !?????? REVOIR

  ccste = parcc_sld(6,1)  !C (windkessel)
  rcste = parcc_sld(7,1)  !R1 (windkessel)
  rcst2 = parcc_sld(8,1)  !R2 (windkessel)

  call runend('sld_ccmtol: ojooooooo, el volumen estaba mal calculado!! revisar esta subru')



  do ivolu=1,1!mcavi_sld !swap all the cavities defined   !INSTEAD just calculate one pressure for the LV and scale it for the RV

     !
     ! PARAMETERS
     !
     !devrais Vient de la courbe "classique" de pression
     ! pfins=-113000.0_rp !Pressure at end systole (end of isovolumetric contraction)

     !  pfinf=-10000.0_rp !Pressure end filling 
     !  tfinf=0.4_rp      !time end fillig (For the first cycle)
     
     !mfili= pfinf/tfinf

     ! pfinr=-25000.0_rp  !Pressure end of isovolum relax

     !  pfirf=0.0_rp !Pressure end rapid filling initial to zero for the firts cycle !?????? REVOIR

     !
     ! DETERMINATION OF THE PHASE (kfase_sld)
     !


     if (kfase_sld==9) then  !Enter filling for the fisrt time in the first cycle 

        kfase_sld=0
        timst_sld(9)=0.0_rp !time end of rapid filling initial to zero for the firts cycle  
        !pfirf=0.0_rp          !Pressure end rapid filling initial to zero for the firts cycle !?????? REVOIR
        !POUR LE MOMENT ON L'IBITIALISE AU DEBUT
        ! Utiliser iwave_sld???


     else if (kfase_sld==0 .AND. dltap_sld(3,ivolu)<pfinf) then     !start of Isovolum contraction 

        kfase_sld=1
        volst(6,ivolu)=volst(2,ivolu)    !keep the volume at the end of Diastole                                  
        dltap_sld(6,ivolu)=dltap_sld(4,ivolu)
        timst_sld(6)=cutim

        !fix the time at which depolarization wave starts : 
        iwave_sld=iwave_sld+1
        aptim(iwave_sld) = cutim  


        !Initialization of the volume for smoothing:                     
        volst(1,ivolu)=volst(2,ivolu)
        volst(11,ivolu)=volst(2,ivolu) 
        dltap_sld(11,ivolu)=dltap_sld(4,ivolu)!300.0_rp

     else if (kfase_sld==1 .AND. dltap_sld(3,ivolu)<pfins) then    !End of isovolu contraction, beginning of ejection

        kfase_sld=2 

        volst(7,ivolu)=volst(2,ivolu) 
        dltap_sld(7,ivolu)=dltap_sld(4,ivolu)   
        timst_sld(7)=cutim 

        ifase_sld(2)=0       !reset the flag to detect the first entry to the isovolum phases     

     else if (kfase_sld==3 .AND. dltap_sld(3,ivolu)>pfinr) then         !End of Isovolum relaxation

        kfase_sld=4

        volst(9,ivolu)=volst(1,ivolu)  !Volume end ejection
        dltap_sld(9,ivolu)=dltap_sld(4,ivolu)  !Pressure end ejection
        timst_sld(9)=cutim

        dltap_sld(12,ivolu)=(dltap_sld(9,ivolu)-dltap_sld(8,ivolu))/(timst_sld(9)-timst_sld(8)) !slope that will be used for rapid filling   
        ifase_sld(2)=0       !reset the flag to detect the first entry to the isovolum phases  

     else if (kfase_sld==4 .AND. dltap_sld(3,ivolu)>pfirf) then         !End of Rapid Filling 
        !stop
        kfase_sld=0
        !pfirf=-4000.0_rp !Pressure end rapid filling  ????????????????? REVOIR

        volst(10,ivolu)=volst(1,ivolu) !Volume end rapid filling
        dltap_sld(9,ivolu)=dltap_sld(4,ivolu)  !Pressure end rapid filling
        timst_sld(9)=cutim         !time rapid filling          

        !stop !If only 1 cycle is simulated

     end if


     ! * * * * * * * * * * 
     !FILLING (ramp function) 
     ! * * * * * * * * * * 

     ! La version originale (valeur du fichier txt):
     if (kfase_sld==0) then

        !gppre = cutim*pfinf/tfinf !linear interpolation 
        ptota_sld(ivolu) = (cutim-timst_sld(9))*(pfinf+pfirf)/tfinf !linear interpolation 

        !Phase ctrl
        if (ifase_sld(1)==0) then  !it's the first iteration of the filling (no volume available for the 1st cycle...)
           ifase_sld(1)=1
        else if (ifase_sld(1)==1) then
           ifase_sld(1)=2
           volst(5,ivolu)=volst(2,ivolu) !keep the volume at the beginning of Diastole 
           dltap_sld(5,ivolu)=dltap_sld(4,ivolu)
           timst_sld(5)=cutim                              
        end if

        !
        dltap_sld(4,ivolu)=ptota_sld(ivolu) !P_new
        dltap_sld(2,ivolu)=dltap_sld(4,ivolu)-dltap_sld(3,ivolu) !deltaP_new


        !update 
        dltap_sld(3,ivolu)=dltap_sld(4,ivolu) 
        dltap_sld(1,ivolu)=dltap_sld(2,ivolu)


     end if !kfase=0

     ! * * * * * * * * * * 
     !ISOVOLUMETRIC (better PDI control theory sould be implemented)
     ! * * * * * * * * * * 

     if (kfase_sld==1 .or. kfase_sld==3) then


        !write(*,*) 'presion entree ',ptota_sld(ivolu)

        !
        !PARAMETERS
        !

        !parametre bidon pour le phase isovolumetrique
        big=1.0_rp
        small=0.0_rp
        reverse=1.0_rp

        !define the volume at which the isovolum phase has to stay either we are in...
        if (kfase_sld==1) then
           volim =volst(6,ivolu)  !...iso contraction (volume of end filling)
           init=-100.0_rp!-1500.0_rp!1500.0_rp
        else if (kfase_sld==3) then  
           volim=volst(8,ivolu)   !...iso relaxation (volume of end ejection)                       
           init=-100.0_rp
        end if



        !Replace "volst(2)" with a projected value 
        !voluf = volst(2,ivolu)+volst(4,ivolu) + 0.5_rp*(volst(4,ivolu)-volst(3,ivolu))/dtinv_sld
        voluf=volst(2,ivolu)


        !Exponential smoothing
        !alpha=0.9
        !voluf=alpha*volst(2,ivolu) + (1.0-alpha)*volst(11,ivolu)
        !volst(11,ivolu)=voluf

        
       ! ** May have to change according to orientation of normal ** 
       !spheroid is TOO BIG (last volume larger than volume at End Filling (or end end ejection))
        if (abs(voluf)>=abs(volim)) then  
       ! if (voluf<volim) then

!write(*,*) 'voluf<volim '
!write(*,*) 'volim target= ',volim
!write(*,*) 'volim target= ',voluf

           !check if we passed on the other side of V_ini
           sidev_sld(2,ivolu)=1.0_rp
           sidei=sidev_sld(2,ivolu)*sidev_sld(1,ivolu)

           if (sidei<0.0_rp) then

              ptota_sld(ivolu)=ptota_sld(ivolu)-init
              bidon=1.0_rp

           else if (abs(volst(2,ivolu))<=abs(volst(1,ivolu))) then ! Too big but Vnew < Vold, do nothing

              ptota_sld(ivolu)=ptota_sld(ivolu)
              bidon=3.0_rp
                                                     
           else if (abs(volst(2,ivolu))>abs(volst(1,ivolu))) then! Too big and Vnew > Vold, decrease pressure  

              ptota_sld(ivolu)=ptota_sld(ivolu)-init
              bidon=4.0_rp

           end if

        else if (abs(voluf)<abs(volim)) then !spheroid is TOO SMALL (last volume smaller than volume at End Filling (or end end ejection))
        !else if (voluf>=volim) then
        
!write(*,*) 'voluf>=volim '
!write(*,*) 'volim target= ',volim
!write(*,*) 'volim target= ',voluf

           !check if we passed on the other side of V_ini
           sidev_sld(2,ivolu)=-1.0_rp
           sidei=sidev_sld(2,ivolu)*sidev_sld(1,ivolu)

           if (sidei<0.0_rp) then

              ptota_sld(ivolu)=ptota_sld(ivolu)-init
              bidon=5.0_rp

           else if (abs(volst(2,ivolu))>abs(volst(1,ivolu))) then ! Too small but Vnew > Vold, do nothing

              ptota_sld(ivolu)=ptota_sld(ivolu)
              bidon=3.0_rp
                                                     
           else ! Too small and Vnew < Vold, increase pressure  

              ptota_sld(ivolu)=ptota_sld(ivolu)+init
              bidon=4.0_rp

           end if !petite boucle

        end if !cases - grande boucle de if

        dltap_sld(4,ivolu)=ptota_sld(ivolu)                    !P_new
        dltap_sld(2,ivolu)=dltap_sld(4,ivolu)-dltap_sld(3,ivolu) !deltaP_new


        !* * 
        !
        !"SMOOTH" PRESSURE
        !
        !* *   ...seems to be useless

        !ptota = ptota_sld(ivolu)+dltap_sld(2,ivolu) + 0.5_rp*(dltap_sld(2,ivolu)-dltap_sld(1,ivolu))/dtinv_sld

        !Expo smooting
        alpha=0.9_rp
        !ptota = alpha*dltap_sld(3,ivolu) + (1.0-alpha)*dltap_sld(11,ivolu)                       
        !ptota = alpha*ptota_sld(ivolu) + (1.0-alpha)*dltap_sld(11,ivolu)
        
        !NOUVELLE IDEE
        dltap_sld(11,ivolu)=(dltap_sld(11,ivolu)+dltap_sld(4,ivolu))/2.0_rp
        volst(11,ivolu)=(volst(11,ivolu)+volst(2,ivolu))/2.0_rp
       
        !No smoothing 
        ptota=ptota_sld(ivolu)
        !somooth nouvelle idee 
        !ptota= dltap_sld(3,ivolu) + 0.6_rp*((volst(2,ivolu)-volst(1,ivolu)) * (dltap_sld(11,ivolu)/volst(11,ivolu)))
   

        

        !dltap_sld(11,ivolu)=ptota !keep smooth  
        dltap_sld(3,ivolu)= ptota_sld(ivolu) !keep unsmooth

        ptota_sld(ivolu)=ptota


        !update variables  
        volst(1,ivolu)=volst(2,ivolu)   !V
        dltap_sld(1,ivolu)=dltap_sld(2,ivolu)  !dP
        volst(3,ivolu)=volst(4,ivolu)  !dV
        sidev_sld(1,ivolu)=sidev_sld(2,ivolu)  !side
        !dltap_sld(3,ivolu)=ptota_sld(ivolu) !P_old

!write(*,*) 'ptota_sld(1) ',ptota_sld(1)
     end if !kfase

     ! * * * * * * * * * * 
     !EJECTION (windkessel)
     ! * * * * * * * * * * 

      if (kfase_sld==2) then 

        if (rcst2 == 0.0_rp) then !Windkessel 2 parameters

           flowq = volst(4,ivolu)*dtinv_sld                ! I=dV/dt
           xfunc = (flowq/ccste) - (dltap_sld(3,ivolu)/(rcste*ccste)) !Attention au signe: +/-(flowq/ccste)???
           !xfunc = (flowq/ccste) - (dltap_sld(3,ivolu)/(rcste*ccste)) !OJO OJO TEST

        else !Windkessel 3 parameters

           cinve = 1.0_rp/ccste  ! 1/C 
           cporr = cinve*(1.0_rp+(rcste/rcst2)) !(1/C)*(1 + (R1/R2))

           flowq = volst(4,ivolu)*dtinv_sld                ! I=dV/dt
           dfldt = (volst(4,ivolu)-volst(3,ivolu))*dtinv_sld ! dI/dt

           xfunc = cporr*flowq + rcste*dfldt - (dltap_sld(3,ivolu)/(rcst2*ccste))                                  

        end if ! windkessel 2 or 3 param 

        !Euler: p_n+1 = p_n + deltaT * f(t_n,P_n) 
        dltap_sld(4,ivolu) =  dltap_sld(3,ivolu) + (1.0_rp/dtinv_sld)*xfunc
        ptota_sld(ivolu)=dltap_sld(4,ivolu)

        !write(*,*) 'ptota_sld(ivolu) ',ptota_sld(ivolu)

        !Volume balance (detect the end of ejection)
        volin=volst(6,ivolu)-volst(5,ivolu)
        volin=0.98_rp*volin  !tolerance on volum in
        volou=volst(2,ivolu)-volst(7,ivolu)


        if (abs(volou) >= abs(volin)) then 
           kfase_sld=3
           !write(*,*) 'FIN WINKESSEL'
           !stop
           volst(8,ivolu)=volst(1,ivolu)  !Volume end ejection
           dltap_sld(8,ivolu)=dltap_sld(4,ivolu)  !Pressure end ejection
           timst_sld(8)=cutim

           sidev_sld=0.0_rp     !reset      
           dltap_sld(1,ivolu)=dltap_sld(4,ivolu)-dltap_sld(3,ivolu) !define initial deltaP 

           ptota_sld(ivolu) = dltap_sld(3,ivolu)                                                                     


        end if

        !update variables 
        dltap_sld(3,ivolu)=dltap_sld(4,ivolu) 
        ! volst(11,ivolu)=volst(1,ivolu) !Initialization of the volume for smoothing

     end if !kfase_sld==2


     ! * * * * * * * * * * 
     ! TEST: USE WINKKESSEL ALWAYS TO CALCUTE P_AO 
     ! * * * * * * * * * * 

     !parameters  
     !ccste=0.00006_rp  !C
     !rcste=5330.0_rp!7500.0_rp !R1   !sensible
     !rcst2=15000.0_rp !R2 !insensible


     !     cinve = 1.0_rp/ccste  ! 1/C 
     !     cporr = cinve*(1.0_rp+(rcste/rcst2)) !(1/C)*(1 + (R1/R2))

     !Euler: p_n+1 = p_n + deltaT * f(t_n,P_n)
     !     flowq = volst(4,ivolu)*dtinv_sld                ! I=dV/dt
     !     if (volst(4,ivolu)>0.0_rp) flowq=0.0_rp
     !     dfldt = (volst(4,ivolu)-volst(3,ivolu))*dtinv_sld ! dI/dt
     !     if (volst(4,ivolu)>0.0_rp) dfldt=0.0_rp

     !     xfunc = cporr*flowq + rcste*dfldt - (dltap_sld(3,ivolu)/(rcst2*ccste))                                  
     !     dltap_sld(4,ivolu) =  dltap_sld(3,ivolu) + (1.0_rp/dtinv_sld)*xfunc
     !     paort =dltap_sld(4,ivolu)

     !     write(6699,*) cutim,paort,flowq
     !     !write(6699,*) cutim,volst(4,ivolu)
     !     flush(6699)

     !   if (kfase_sld==2) then  
     !     !update variables 
     !     dltap_sld(3,ivolu)=dltap_sld(4,ivolu) 
     !     volst(11,ivolu)=volst(1,ivolu) !Initialization of the volume for smoothing
     !    end if

     ! * * * * * * * * * * 
     ! FIN TEST
     ! * * * * * * * * * * 


     ! * * * * * * * * * * 
     !RAPID FILLING (just after isovolume relaxation, for cycle .ne. 1)
     ! * * * * * * * * * * 


     if (kfase_sld==4) then 


        !ptota_sld(ivolu) = dltap_sld(3,ivolu) + dltap_sld(12,ivolu)*(1.0_rp/dtinv_sld) !linear interpolation 
        ptota_sld(ivolu) = dltap_sld(3,ivolu) + (dltap_sld(12,ivolu)/2.5_rp)*(1.0_rp/dtinv_sld)


        !update variables  
        dltap_sld(3,ivolu)=ptota_sld(ivolu) !P_old


     end if !kfase_sld==4

     !
     ! Broadcast values to the slaves
     !  
     call pararr('BCT',0_ip,1_ip,ptota_sld(ivolu))

  end do !ivolu

101 format (9(e13.6,' ')) 
end subroutine sld_ccmtol
