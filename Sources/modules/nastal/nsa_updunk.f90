!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_updunk.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Update unknowns
!> @details Update unknowns
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_updunk(itask)
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal

  implicit none
  integer(ip), intent(in) :: itask !> Who is calling 
  integer(ip)             :: ipoin,kpoin,itotv,idime,idofn,ievat,itime,jtinn
  real (rp)               :: varia, densitot, dt
  character               :: fnp1*3, fnp*72
  integer(ip)             :: irestart,itime_scheme,itime_last
  
  if(INOTMASTER) then

     select case (itask)
     case(1)   
        
        do ipoin = 1,npoin
           veloc(1:ndime,ipoin,ITER_AUX) = veloc(1:ndime,ipoin,TIME_N) 
           umome(1:ndime,ipoin,ITER_AUX) = umome(1:ndime,ipoin,TIME_N) 
           press(        ipoin,ITER_AUX) = press(        ipoin,TIME_N)  
           densi(ipoin,ITER_AUX) = densi(ipoin,TIME_N)
           tempe(ipoin,ITER_AUX) = tempe(ipoin,TIME_N)
           energ(ipoin,ITER_AUX) = energ(ipoin,TIME_N)
           visco(ipoin,ITER_AUX) = visco(ipoin,TIME_N)
        end do

     case(2)
        !
        ! Assign u(n,i,j) <-- u(n,i-1,j)
        !
        do ipoin = 1,npoin
           veloc(1:ndime,ipoin,ITER_K) = veloc(1:ndime,ipoin,ITER_AUX)
           press(ipoin,ITER_K)         = press(ipoin,ITER_AUX) 
           veloc(1:ndime,ipoin,ITER_K) = veloc(1:ndime,ipoin,ITER_AUX) 
           umome(1:ndime,ipoin,ITER_K) = umome(1:ndime,ipoin,ITER_AUX)  
           densi(ipoin,ITER_K)         = densi(ipoin,ITER_AUX) 
           tempe(ipoin,ITER_K)         = tempe(ipoin,ITER_AUX) 
           energ(ipoin,ITER_K)         = energ(ipoin,ITER_AUX) 
           visco(ipoin,ITER_K)         = visco(ipoin,ITER_AUX) 
        end do

!!!!!        call nsa_setvar(three,one)

     case(3)
        !
        ! Assign u(n,i,j-1) <-- Solver unknown, update conservative variables in inner iteration
        ! according to the time scheme (INCREMENTAL, see also case (9))
        do ipoin = 1,npoin

           do idime= 1,ndime
              ievat = (ipoin-1)*ndofn_nsa + idime
              umome(idime,ipoin,ITER_K) = unkno(ievat)
           end do
           densi(ipoin,ITER_K) =  unkno(ievat+1) 
           energ(ipoin,ITER_K) =  unkno(ievat+2)
        end do

     case(4)
        !
        ! Assign u(n,i,j-1) <-- u(n,i,j)
        !
        do ipoin = 1,npoin
           veloc(1:ndime,ipoin,ITER_AUX) = veloc(1:ndime,ipoin,ITER_K) 
           umome(1:ndime,ipoin,ITER_AUX) = umome(1:ndime,ipoin,ITER_K) 
           press(        ipoin,ITER_AUX) = press(        ipoin,ITER_K)  
           densi(        ipoin,ITER_AUX) = densi(ipoin,ITER_K) 
           tempe(        ipoin,ITER_AUX) = tempe(ipoin,ITER_K) 
           energ(        ipoin,ITER_AUX) = energ(ipoin,ITER_K) 
        end do

!!!!!!        write(6022,*) cutim,densi(390,ITER_AUX),densi(360,ITER_AUX)

     case(5)
        !
        ! Obtain u(n,*,*) for the Crank-Nicholson method and assign
        ! u(n-1,*,*) <-- u(n,*,*)
        !     

        ! This first loop is only active for BDF schemes through the value of tiaor
        do itime=2+kfl_tiaor_nsa,4,-1 
           do ipoin = 1,npoin
              veloc(1:ndime,ipoin,itime) = veloc(1:ndime,ipoin,itime-1) 
              umome(1:ndime,ipoin,itime) = umome(1:ndime,ipoin,itime-1) 
              press(        ipoin,itime) = press(        ipoin,itime-1)  
              densi(ipoin,itime) = densi(ipoin,itime-1) 
              tempe(ipoin,itime) = tempe(ipoin,itime-1) 
              energ(ipoin,itime) = energ(ipoin,itime-1) 
           end do           
        end do
     
        do ipoin = 1,npoin
           veloc(1:ndime,ipoin,TIME_N) = veloc(1:ndime,ipoin,ITER_AUX) 
           umome(1:ndime,ipoin,TIME_N) = umome(1:ndime,ipoin,ITER_AUX) 
           press(        ipoin,TIME_N) = press(        ipoin,ITER_AUX)  
           densi(ipoin,TIME_N) = densi(ipoin,ITER_AUX) 
           tempe(ipoin,TIME_N) = tempe(ipoin,ITER_AUX) 
           energ(ipoin,TIME_N) = energ(ipoin,ITER_AUX) 
        end do


     case(7)
        !
        ! Assign unkno <-- (rho,U,E)(TIME_N) or (p,U,T)(TIME_N) or (rho,U,T)(TIME_N) 
        !        
        do ipoin = 1,npoin
           ievat=(ipoin-1)*ndofn_nsa
           if (kfl_unkse_nsa == 1) then        ! primitive set (p,u,T)
              do idime=1,ndime
                 ievat=ievat+1
                 unkno(ievat)=umome(idime,ipoin,TIME_N) / densi(ipoin,TIME_N)
              end do
              unkno(ievat+1)=press(ipoin,TIME_N)   
              unkno(ievat+2)=tempe(ipoin,TIME_N)
           else if (kfl_unkse_nsa == 0) then   ! conservative set (rho,U,E)
              do idime=1,ndime
                 ievat=ievat+1
                 unkno(ievat)=umome(idime,ipoin,TIME_N) 
              end do
              
              unkno(ievat+1)=densi(ipoin,TIME_N)   
              unkno(ievat+2)=energ(ipoin,TIME_N)
              
           else if (kfl_unkse_nsa >= 10) then  ! conservative/heat set (rho,U,T)
              
              densitot = densi(ipoin,TIME_N) + rekee_nsa(ndime+1,ipoin)
              
              if(kfl_ncons_nsa == 0) then
                 !
                 ! Conservative set
                 !
                 do idime=1,ndime
                    ievat=ievat+1
                    unkno(ievat)=umome(idime,ipoin,TIME_N) 
                 end do
              else
                 !
                 ! Non-conservative set
                 !
                 do idime=1,ndime
                    ievat=ievat+1
                    unkno(ievat)=umome(idime,ipoin,TIME_N) / densitot
                 end do
              end if
              unkno(ievat+1)=densi(ipoin,TIME_N)
              unkno(ievat+2)=tempe(ipoin,TIME_N)
           end if
        end do
        

     case(8)
        !
        ! Assign unkno <-- (rho,U,E)(ITER_K) or (p,U,T)(ITER_K) or (rho,U,T)(ITER_K) 
        !        
        do ipoin = 1,npoin

           ievat=(ipoin-1)*ndofn_nsa
           do idime=1,ndime
              ievat=ievat+1
              dunkn_nsa(ievat)=0.0_rp
           end do
           dunkn_nsa(ievat+1)=0.0_rp
           dunkn_nsa(ievat+2)=0.0_rp

           ievat=(ipoin-1)*ndofn_nsa
           if (kfl_unkse_nsa == 1) then        ! primitive set (p,u,T)
              do idime=1,ndime
                 ievat=ievat+1
                 unkno(ievat)=umome(idime,ipoin,ITER_K) / densi(ipoin,ITER_K)
              end do
              unkno(ievat+1)=press(ipoin,ITER_K)   
              unkno(ievat+2)=tempe(ipoin,ITER_K)
           else if (kfl_unkse_nsa == 0) then   ! conservative set (rho,U,E)
              do idime=1,ndime
                 ievat=ievat+1
                 unkno(ievat)=umome(idime,ipoin,ITER_K) 
              end do
              
              unkno(ievat+1)=densi(ipoin,ITER_K)   
              unkno(ievat+2)=energ(ipoin,ITER_K)
              
           else if (kfl_unkse_nsa >= 10) then  ! conservative/heat set (rho,U,T)
              
              densitot = densi(ipoin,ITER_K) + rekee_nsa(ndime+1,ipoin)
              
              if(kfl_ncons_nsa == 0) then
                 !
                 ! Conservative set
                 !
                 do idime=1,ndime
                    ievat=ievat+1
                    unkno(ievat)=umome(idime,ipoin,ITER_K) 
                 end do
              else
                 !
                 ! Non-conservative set
                 !
                 do idime=1,ndime
                    ievat=ievat+1
                    unkno(ievat)=umome(idime,ipoin,ITER_K) / densitot
                 end do
              end if
              unkno(ievat+1)=densi(ipoin,ITER_K)
              unkno(ievat+2)=tempe(ipoin,ITER_K)
           end if
        end do
        

     case(10)

        ! Initialization of fields for possible coupling with other modules

        do ipoin = 1,npoin
           veloc(1:ndime,ipoin,ITER_K) = veloc(1:ndime,ipoin,TIME_N) 
           veloc(1:ndime,ipoin,ITER_AUX) = veloc(1:ndime,ipoin,TIME_N) 
           umome(1:ndime,ipoin,ITER_K) = umome(1:ndime,ipoin,TIME_N)
           umome(1:ndime,ipoin,ITER_AUX) = umome(1:ndime,ipoin,TIME_N) 

           press(ipoin,ITER_K) = press(ipoin,TIME_N)
           densi(ipoin,ITER_K) = densi(ipoin,TIME_N)
           tempe(ipoin,ITER_K) = tempe(ipoin,TIME_N)
           energ(ipoin,ITER_K) = energ(ipoin,TIME_N)
           visco(ipoin,ITER_K) = visco(ipoin,TIME_N)

           press(ipoin,ITER_AUX) = press(ipoin,TIME_N)
           densi(ipoin,ITER_AUX) = densi(ipoin,TIME_N)
           tempe(ipoin,ITER_AUX) = tempe(ipoin,TIME_N)
           energ(ipoin,ITER_AUX) = energ(ipoin,TIME_N)
           visco(ipoin,ITER_AUX) = visco(ipoin,TIME_N)
        end do

     case(11)
        !
        ! Correcting unkno when the delta form is used
        !

        !
        ! No pseudo time used: update from the last time
        !
        itime_last  = TIME_N
        itime_scheme= TIME_N  ! this is required to avoid forbidden memory acces of the global vectors
        if (kfl_tisch_nsa == 2) itime_scheme= TIME_N_MINUS_1        
           
        if (kfl_pseud_nsa == 1) then
           !
           ! Pseudo time used: update from the last iteration
           !
           itime_last  = ITER_K  ! this is required to avoid forbidden memory acces of the global vectors
           itime_scheme= ITER_K  ! this is required to avoid forbidden memory acces of the global vectors
           if (kfl_tisch_nsa == 2) then
              call runend('NSA_UPDUNK: CN SCHEMES AND PSEUDO-TIME NOT YET PROGRAMMED.')
           end if
        end if

        if (kfl_linea_nsa == 1) then
           ! Jacobi updates are from itime_last
           do ipoin = 1,npoin
              !         
              ! There is a "-" because of how the bdf parameters are defined
              !         
              do idime= 1,ndime
                 ievat = (ipoin-1)*ndofn_nsa + idime
!                 unkno(ievat)= umome(idime,ipoin,itime_last) +dunkn_nsa(ievat)
                 unkno(ievat)= &
                      (- pabdf_nsa(2) * umome(idime,ipoin,itime_last) &
                      - pabdf_nsa(3) * umome(idime,ipoin,itime_scheme))/pabdf_nsa(1)  &
                      + dunkn_nsa(ievat)
              end do

              unkno(ievat+1)= (- pabdf_nsa(2) * densi(ipoin,itime_last) &
                   - pabdf_nsa(3) * densi(ipoin,itime_scheme) )/pabdf_nsa(1)&
                   + dunkn_nsa(ievat+1) 
              unkno(ievat+2)= (- pabdf_nsa(2) * energ(ipoin,itime_last) &
                   - pabdf_nsa(3) * energ(ipoin,itime_scheme) )/pabdf_nsa(1)&
                   + dunkn_nsa(ievat+2)

!              unkno(ievat+1)= densi(ipoin,itime_last) + dunkn_nsa(ievat+1) 
!              unkno(ievat+2)= energ(ipoin,itime_last) + dunkn_nsa(ievat+2)
           end do
        else if (kfl_linea_nsa == 2) then
           ! Newton Raphson updates are from ITER_K, stored in unkno
           do ipoin = 1,npoin
              do idofn= 1,ndofn_nsa
                 ievat = (ipoin-1)*ndofn_nsa + idofn
                 unkno(ievat)= unkno(ievat) + dunkn_nsa(ievat)
              end do
           end do           
        end if

     case(12)
        !
        ! Assigning rhsmo_nsa to rhsmo 
        !
        do ipoin = 1,npoin
           do idime= 1,ndime
              ievat = (ipoin-1)*ndofn_nsa + idime
              rhsid(ievat)= resmo_nsa(ievat)
           end do
           rhsid(ievat+1)= resmo_nsa(ievat+1)
           rhsid(ievat+2)= resmo_nsa(ievat+2)
        end do
        
     end select

  end if

end subroutine nsa_updunk

