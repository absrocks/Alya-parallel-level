!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_comapp.f90
!> @author  Various
!> @brief   Compute of stimulus current
!> @details Setup the appfi_exm (applied current field) to the EP problem\n
!! The parameters are being read at exm_reaphy.f90 \n
!> @} 
!!-----------------------------------------------------------------------
subroutine exm_comapp

  use      def_master
  use      def_domain
  use      def_exmedi
  implicit none

  integer(ip) :: ipoin,istim,iboun,pblty,pnodb,inodb
  real(rp)    :: dicen,xicen,yicen,zicen,batim
  integer(ip) :: istim_compute,kmodel,nauxi
  real(rp)    :: sqrea(15000)
  integer(ip) :: kstim(15000)

  if (INOTMASTER) then

     batim = cutim  
     !  ltime = cutim -dtime                        !!! <<<---  Iapp is that of the CURRENT time step
     appfi_exm = 0.0_rp

     if (modst_exm < 0_ip) then
        !
        ! Stimuli in a field
        !
        do ipoin= 1,npoin
           if (xfiel(-modst_exm) % a(2,ipoin,1) .ge. 0.0_rp) then
              if (batim .ge. xfiel(-modst_exm) % a(2,ipoin,1)) then              
                 appfi_exm(ipoin) = xfiel(-modst_exm) % a(1,ipoin,1)
                 if (batim .gt. (xfiel (-modst_exm) % a (2,ipoin,1) + xfiel (-modst_exm) % a (3,ipoin,1))) then !Correction by Francesc
                 !if (batim .gt. xfiel(-modst_exm) % a(3,ipoin)) then
                    appfi_exm(ipoin) = 0.0_rp
                 end if
              end if
           end if
        end do

     else

        if (aploo_exm(1) > 0.0_rp) then
           if (batim > aploo_exm(3)) then
              !
              ! reset loop
              !
              kstim= 0_ip
              aptim= 0.0_rp
              aploo_exm(3) = aploo_exm(1) + batim - (batim-aploo_exm(3))
              if (aploo_exm(2) > 0.0_rp) then
                 !
                 ! set loop time
                 !
                 aptim = aptim + aploo_exm(2)
              end if
           end if
        end if

        if (kfl_appty_exm == 1_ip) then
           ! decay time lapse considered, good for TT, Ohara, etc.

           do istim= 1, nstim_exm

              istim_compute = 0_ip
              if (kfl_ptrig_exm == 1_ip) then
                 if (epres > aptim(istim)) then
                    istim_compute = 1_ip 
                    kfl_ptrig_exm = -1_ip
                    aptim(istim)  = batim  ! start the chrono
                 end if
              else if (kfl_ptrig_exm == -1_ip) then
                 if (batim <= (aplap_exm(istim)+ aptim(istim))) then
                    istim_compute = 1_ip           
                 end if
              else if (kfl_ptrig_exm == 0_ip) then
                 if (batim >= aptim(istim) .and. batim <= (aplap_exm(istim)+ aptim(istim))) then
                    istim_compute = 1_ip          
                 end if
              end if

              if (istim_compute == 1_ip) then 
                   
                 qneto_exm = 0.0_rp       
                 sqrea(istim)= aprea_exm(istim)*aprea_exm(istim)    
                 do ipoin= 1,npoin
                    xicen = coord(1,ipoin)-apcen_exm(1,istim)
                    yicen = coord(2,ipoin)-apcen_exm(2,istim)
                    zicen = 0.0_rp
                    if (ndime == 3) zicen = coord(ndime,ipoin)-apcen_exm(ndime,istim)

                    dicen = xicen*xicen + yicen*yicen + zicen*zicen

                    if (dicen <= sqrea(istim)) then
                       lapno_exm(ipoin) = kfl_appli_exm
                       appfi_exm(ipoin) = apval_exm(istim)                          
                       !!!!! write(666,*) ipoin, '  -80.0  0.000  0.002  '
                    end if
                 end do
                 !else if (batim < aptim(istim) .and. batim > (aplap_exm(istim)+ aptim(istim))) then
                 !appfi_exm(1:npoin,3)= 0._rp           
                 !appfi_exm(1:npoin,1)= 0._rp          

              end if
           end do

           if (nstim_exm < 0_ip) then
              !
              !
              !
              do iboun= 1,nboun
                 pblty = ltypb(iboun) 
                 pnodb = nnode(pblty)

                 if (lbset(iboun)== nstis_exm) then

                    do inodb= 1,pnodb
                       ipoin= lnodb(inodb,iboun)

                       if (batim >= aptim(istim) .and. batim <= (aplap_exm(istim)+ aptim(istim))) then
                          lapno_exm(ipoin) = kfl_appli_exm
                          appfi_exm(ipoin) = apval_exm(1)
                       end if

                    end do
                 end if
              end do

           end if

        else   if (kfl_appty_exm == 2_ip) then
           ! decay time lapse not considered, flash-like stimuly, good for FHN, but then, the result depends on the time step!

           istim_compute = 0_ip
           kstim=0_ip
           if (kfl_appli_exm == 200_ip) then
              nauxi= nstim_exm        
              if (nstim_exm < 0) nauxi=1 

              do istim= 1, nauxi
                 if (kfl_ptrig_exm == 1_ip) then
                    if (epres > aptim(istim)) then
                       istim_compute = 1_ip  
                       kfl_ptrig_exm = -1_ip
                       aptim(istim)  = batim  ! start the chrono
                    end if
                 else if (kfl_ptrig_exm == -1_ip) then
                    !              if (batim <= (aplap_exm(istim)+ aptim(istim))) then
                    !                 istim_compute = 1           
                    !              end if
                    !  no lapse considered, just a flash 
                    kfl_ptrig_exm = 0_ip
                 else if (kfl_ptrig_exm == 0_ip) then
                    if (aptim(istim) >= 0_ip) then
                       if (batim >= aptim(istim) ) then
                          kstim(istim) = istim
                          aptim(istim) = -1.0_rp
                          istim_compute = istim_compute + 1                    ! count the starting stimuli
                       end if
                    end if
                 end if
              end do

              if (istim_compute == 0_ip) return

              if (nstim_exm < 0_ip) then
                 !
                 !
                 !
                 do iboun= 1,nboun
                    pblty = ltypb(iboun) 
                    pnodb = nnode(pblty)

                    if (lbset(iboun)== nstis_exm) then

                       do inodb= 1,pnodb
                          ipoin= lnodb(inodb,iboun)

                          elmag(ipoin,1)= apval_exm(1)      
                          elmag(ipoin,2)= apval_exm(1)      
                          elmag(ipoin,3)= apval_exm(1)

                       end do
                    end if
                 end do

                 return

              end if

              ! Starting potential
              do istim= 1,nstim_exm
                 sqrea(istim)= 0
                 if (kstim(istim) > 0) then
                    sqrea(istim)= aprea_exm(kstim(istim))*aprea_exm(kstim(istim))
                 end if
              end do

              do ipoin= 1,npoin

                 kmodel= kfl_cellmod(nodemat(ipoin))           
                 ! when NO_MODEL, leave elmag equal to zero
!!!           if (kmodel > 0) then

                 do istim= 1,nstim_exm

                    if (kstim(istim) > 0) then

                       xicen = coord(1,ipoin)-apcen_exm(1,kstim(istim))
                       yicen = coord(2,ipoin)-apcen_exm(2,kstim(istim))
                       zicen = 0.0_rp
                       if (ndime == 3) zicen = coord(ndime,ipoin)-apcen_exm(ndime,kstim(istim))

                       dicen = xicen*xicen + yicen*yicen + zicen*zicen

                       if (dicen <= sqrea(istim)) then
                          elmag(ipoin,1)= apval_exm(kstim(istim))
                          elmag(ipoin,2)= apval_exm(kstim(istim))
                          elmag(ipoin,3)= apval_exm(kstim(istim))
                          appfi_exm(ipoin)= apval_exm(kstim(istim))
                       end if

                    end if
                 end do

                 !!           end if
              end do
           end if

        end if

     end if

  end if

end subroutine exm_comapp
