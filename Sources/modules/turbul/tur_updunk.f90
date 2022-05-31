!-----------------------------------------------------------------------
!> @addtogroup Turbul
!> @{
!> @file    tur_updunk.f90
!> @author  Guillaume Houzeaux
!> @date    02/11/2015
!> @brief   Solution updates
!> @details Solution updates:
!>          do time
!>             tur_begste (itask=1) ..................... (:,2) <= (:,3)
!>             do outer
!>                tur_begite (itask=2) .................. (:,1) <= (:,2)
!>                do inner
!>                   tur_endite (itask=3, inner loop) ... (:,1) <= UNKNO
!>                end do
!>                tur_endite (itask=4, outer loop) ...... (:,2) <= (:,1)
!>             end do
!>             tur_endste (itask=5) ..................... (:,3) <= (:,1)
!>          end do
!> @} 
!-----------------------------------------------------------------------

subroutine tur_updunk(itask)
  use def_parame
  use def_master
  use def_domain
  use def_turbul
  use mod_postpr
  use def_kermod
  use mod_ker_proper,         only : ker_updpro
  use mod_output_postprocess, only : output_postprocess_check_variable_postprocess
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,iturb,itotn,itime
  real(rp)                :: rela1

  if( INOTMASTER ) then

     select case (itask)

     case(1_ip) 
        !
        ! Assign f(n,0,*) <-- f(n-1,*,*), initial guess for outer iterations
        !     
        do ipoin=1,npoin
           do iturb=1,nturb_tur
              untur(iturb,ipoin,2) = untur(iturb,ipoin,nprev_tur)        
           end do
        end do

     case(2_ip)
        !
        ! Assign f(n,i,0) <-- f(n,i-1,*), initial guess for inner iterations
        !
        do ipoin=1,npoin  
           do iturb=1,nturb_tur
              untur(iturb,ipoin,1) = untur(iturb,ipoin,2)
           end do
        end do 


     case(3_ip)
        !
        ! Assign f(n,i,j-1) <-- f(n,i,j), update of the untur
        !
        if(iunkn_tur==nturb_tur) then
           if( output_postprocess_check_variable_postprocess(17_ip) ) then
              do ipoin=1,npoin
                 unold_tur(nturb_tur+1,ipoin) = turmu(ipoin)
              end do
           end if
        end if

        if(kfl_algor_tur==1) then
           do ipoin=1,npoin 
              untur(iunkn_tur,ipoin,1)=unkno(ipoin)
           end do
        else
           do ipoin=1,npoin 
              itotn=(ipoin-1)*nturb_tur
              do iunkn_tur=1,nturb_tur
                 itotn=itotn+1
                 untur(iunkn_tur,ipoin,1)=unkno(itotn)
              end do
           end do
        end if


     case(4_ip)
        !
        ! Assign f(n,i-1,*) <-- f(n,i,*)
        !   
        if(    output_postprocess_check_variable_postprocess(15_ip).and.&
             & output_postprocess_check_variable_postprocess(16_ip) ) then
           do ipoin=1,npoin
              do iturb=1,nturb_tur
                 unold_tur(iturb,ipoin) = untur(iturb,ipoin,2)       
              end do
           end do
        end if

        do ipoin=1,npoin
           do iturb=1,nturb_tur
              untur(iturb,ipoin,2) = untur(iturb,ipoin,1)            
           end do
        end do

     case(5_ip) ! endste
        !
        ! Obtain f(n,*,*) for the Crank-Nicolson method and assign
        ! f(n-1,*,*) <-- f(n,*,*)
        !     
        if( kfl_tisch_tur == 1 .and. kfl_tiacc_tur==2 ) then
           !
           ! Crank-Nicolson method 
           !
           untur(1:nturb_tur,1:npoin,1) = 2.0_rp * untur(1:nturb_tur,1:npoin,1) - untur(1:nturb_tur,1:npoin,3)
        else if( kfl_tisch_tur == 2 ) then
           !
           ! BDF scheme
           !
           do itime = 2+kfl_tiaor_tur,4,-1  ! ....., 5=4, 4=3   -- 3=1 is done later
              untur(1:nturb_tur,1:npoin,itime) = untur(1:nturb_tur,1:npoin,itime-1)
           end do
        end if
     
        untur(1:nturb_tur,1:npoin,3) = untur(1:nturb_tur,1:npoin,1)
        

     case(6_ip) 
        !
        ! Assign f(n,0,*) <-- f(n-1,*,*) when using restart file
        !
        do ipoin=1,npoin
           do iturb=1,nturb_tur
              untur(iturb,ipoin,1) = untur(iturb,ipoin,nprev_tur)
           end do
        end do

     case(7_ip)
        !
        ! Just after the solver, impose threshold according to clipping strategy
        !
        call tur_clippi()

     case(8_ip)
        !
        ! Relax UNKNO
        !
        if(relax_tur/=1.0_rp) then
           rela1=1.0_rp-relax_tur
           if(kfl_algor_tur==1) then
              do ipoin=1,npoin 
                 if(kfl_fixno_tur(1,ipoin,iunkn_tur)<=0) then
                    unkno(ipoin)=relax_tur*unkno(ipoin)&
                         & + rela1*untur(iunkn_tur,ipoin,1)   
                 end if
              end do
           else
              do ipoin=1,npoin
                 itotn=(ipoin-1)*nturb_tur
                 do iunkn_tur=1,nturb_tur
                    itotn=itotn+1
                    if(kfl_fixno_tur(1,ipoin,iunkn_tur)<=0) then
                       unkno(itotn)=relax_tur*unkno(itotn)&
                            & + rela1*untur(iunkn_tur,ipoin,1)
                    end if
                 end do
              end do
           end if
        end if

     case(9_ip)
        !
        ! Solver initial guess
        !
        if(kfl_algor_tur==1) then
           do ipoin=1,npoin
              unkno(ipoin)=untur(iunkn_tur,ipoin,1)
              if (kfl_fixno_tur(1,ipoin,iunkn_tur)>0) unkno(ipoin) = bvess_tur(1,ipoin,iunkn_tur) ! Force initial guess to satisfy bcs
           end do
        else
           do ipoin=1,npoin
              itotn=(ipoin-1)*nturb_tur
              do iunkn_tur=1,nturb_tur
                 itotn=itotn+1
                 unkno(itotn)=untur(iunkn_tur,ipoin,1)
                 if (kfl_fixno_tur(1,ipoin,iunkn_tur)>0) unkno(itotn) = bvess_tur(1,ipoin,iunkn_tur) ! Force initial guess to satisfy bcs
              end do
           end do
        end if

     case(10_ip)
        !
        ! Initial guess
        !
        do ipoin=1,npoin
           do iturb=1,nturb_tur
              untur(iturb,ipoin,1)=untur(iturb,ipoin,nprev_tur)
           end do
        end do

     case(11_ip)
        !
        ! Put turbulent viscosity to zero (in case we have a delay)
        !
        do ipoin=1,npoin
           turmu(ipoin)=0.0_rp
        end do
     case(12_ip) 
        !
        ! updates unknown at the end of inner iteration (for each unknown)
        !
        if (iunkn_tur<=2.and.niter_tur.gt.1) then
     
           rela1=relax_tur

           do ipoin=1,npoin
              untur(iunkn_tur,ipoin,2) = rela1*untur(iunkn_tur,ipoin,1) +(1.0_rp-rela1)*untur(iunkn_tur,ipoin,2)
              untur(iunkn_tur,ipoin,1) = untur(iunkn_tur,ipoin,2)
              !   initial guess for epsilon, keeping the same nut
              !                 untur(2,ipoin,1) = (1.0_rp-rela1)*untur(2,ipoin,2)+rela1*untur(2,ipoin,2)*untur(1,ipoin,1)*untur(1,ipoin,1)/(untur(1,ipoin,2)*untur(1,ipoin,2) )
           end do
           call tur_updeda() ! update eddy at gauss point
         
           call ker_updpro(ITASK_ENDINN)           

        else
           do ipoin=1,npoin
              untur(iunkn_tur,ipoin,2) = untur(iunkn_tur,ipoin,1)
           end do
        end if
       
     end select

  end if

end subroutine tur_updunk
