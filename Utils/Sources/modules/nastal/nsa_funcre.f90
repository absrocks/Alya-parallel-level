subroutine nsa_funcre(&
     itask,idofn,ndime,venew,bvess,kfixn,palaw,npara,ifuge,kfixi,rtini,rtifi,timev)
  !------------------------------------------------------------------------
  !
  ! This subroutine computes transient boundary conditions
  !
  ! CAVEAT: THIS IS DIFFERENT THAN KERNEL'S FUNCRE!!! 
  !
  !------------------------------------------------------------------------
  use def_kintyp

  use      def_nastal

  implicit none
  integer(ip), intent(in)  :: npara,ifuge,kfixi,ndime,itask,idofn
  real(rp),    intent(in)  :: palaw(npara),timev
  integer(ip)              :: ipara,iftbc,kfixn(ndime+2)
  real(rp)                 :: timea,timeb,funca,funcb,zerom,timec
  real(rp)                 :: timei,timef,rtini,rtifi,tirel,tista,tiend,tirep
  real(rp)                 :: venew,bvess(ndime+2),xx_funcre,bvref(ndime+2)
  integer(ip)              :: ifunc,idata
  real(rp)                 :: t1,t2,p1,p2,m,b

  zerom=epsilon(1.0_rp)

  tista= palaw(11)
  tiend= palaw(12)
  tirep= palaw(13)
  bvref =palaw(15:15+ndime+2)

  if (itask == 0_ip) then                    ! initial fixities
     
     if (kfixi == 1 .or. kfixi==3) then      ! valve-like behavior:
        kfixn(1:ndime)= 1                    !   close the valve (viscous)
        kfixn(ndime+1)= 0                    
     else if (kfixi == 2 .or. kfixi==4) then               
        kfixn(      1)= 2                    !   close the valve (inviscid)
        kfixn(2:ndime)= 0
        kfixn(ndime+1)= 0                    
     end if
     
  else if (itask == 1_ip) then               ! check fixities 
     if ((timev+zerom) > rtini) then   
        if (timev < rtifi) then                    ! within the time lapse
           if (kfixi == 1 .or. kfixi==2) then      ! p-valve-like behavior:
              kfixn(1:ndime)= 0                    !   open the valve
              kfixn(ndime+1)= 2                    
              bvess(ndime+1)= bvref(ndime+1)
           else if (kfixi == 3 .or. kfixi==4) then      ! v-valve-like behavior:
              kfixn(1:ndime)= 1                         !   open the valve
              bvess(1:ndime)= bvref(1:ndime)
           end if
        else                                       ! outside the time lapse
           if (kfixi == 1) then                    ! valve-like behavior:
              kfixn(1:ndime)= 1                    !   p-valve: close the valve (viscous)
              kfixn(ndime+1)= 0                    
              bvess(1:ndime)= 0.0_rp
           else if (kfixi == 2) then               
              kfixn(      1)= 2                    !   p-valve: close the valve (inviscid)
              kfixn(2:ndime)= 0
              kfixn(ndime+1)= 0                               
              bvess(1:ndime)= 0.0_rp
           else if (kfixi == 3) then               
              kfixn(1:ndime)= 1                    !   v-valve: close the valve (viscous)
              kfixn(ndime+1)= 0                    
              bvess(1:ndime)= 0.0_rp
           else if (kfixi == 4) then               
              kfixn(      1)= 2                    !   v-valve: close the valve (inviscid)
              kfixn(2:ndime)= 0
              kfixn(ndime+1)= 0                               
              bvess(1:ndime)= 0.0_rp
           end if
        end if
     end if


  else if (itask == 2_ip) then

     venew=bvess(idofn)
     
     xx_funcre= 1.0_rp
     
     if(ifuge==0) then
        !
        ! No time dependence 
        !
        xx_funcre=1.0_rp
        
     else if(ifuge==1) then
        !
        ! Polynomial evolution
        !
        iftbc= 0
        if ((timev+zerom) > rtini) then   
           if (timev < rtifi) then  ! within the time lapse
              tirel= timev - rtini
              iftbc= 1
           end if
        end if
     
        if (iftbc == 0) return

        xx_funcre= palaw(1) &
             + palaw(2) * tirel &
             + palaw(3) * tirel * tirel &
             + palaw(4) * tirel * tirel * tirel

     else if(ifuge==2) then
        !
        ! Periodic evolution
        !
        if(palaw(1)-zerom<=timev.and.timev<=palaw(2)+zerom) then  !Normal evolution
           xx_funcre=palaw(3)*cos(palaw(4)*timev+palaw(5))+palaw(6)
        else if (timev>palaw(2)+zerom) then
           timea=palaw(2)
           xx_funcre=palaw(3)*cos(palaw(4)*timea+palaw(5))+palaw(6)
        else if (timev<palaw(1)-zerom) then
           timea=palaw(1)
           xx_funcre=palaw(3)*cos(palaw(4)*timea+palaw(5))+palaw(6)
        end if
     else if(ifuge==3) then
        !
        ! Discrete evolution
        !
        timei=palaw(1)
        timef=palaw(2)
        ifunc=int(palaw(3))
        
        if (timei<timev .and. timev<timef) then
           do idata= 1,npara-1
              if (timev >= tload_nsa(ifunc)%a(ndime+1,idata) .and. timev < tload_nsa(ifunc)%a(ndime+1,idata+1) ) then
                    t1= tload_nsa(ifunc)%a(ndime+1,idata)  !start window time
                    t2= tload_nsa(ifunc)%a(ndime+1,idata+1) !end window time
                    p1= tload_nsa(ifunc)%a(idofn,idata  ) !start value
                    p2= tload_nsa(ifunc)%a(idofn,idata+1) !end value
                    m=(p2-p1)/(t2-t1)
                    b=p1-m*t1
                    xx_funcre=m*timev+b                              ! return value
              end if
           end do
        end if
        
     else if(ifuge==4) then
        !
        ! Special function to change boundary values
        !
        xx_funcre=palaw(2)
        
     else if(ifuge==5) then
        !
        ! Marek Prymon's function
        !
        if(timev<palaw(1)) then
           xx_funcre = palaw(2)
        else if(timev<palaw(3)) then
           xx_funcre = palaw(4)-palaw(5)*timev
        else
           xx_funcre = palaw(6)
        end if
        
     end if

     !
     !  Compute new value
     !     
     venew= bvess(idofn)*xx_funcre

  else if (itask == 10_ip) then                    ! update time lapses
     
     if ((timev+zerom) > rtini) then   
        if (timev < rtifi) then                    ! within the time lapse
           ! ok, do nothing
        else                                       ! outside the time lapse
           rtini= rtifi + tirep
           rtifi= rtini + (tiend - tista)        
        end if
     end if
     
  end if
  
end subroutine nsa_funcre
