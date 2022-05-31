!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_updbcs.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Update boundary conditions
!> @details Update boundary conditions
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_updbcs(itask)
  use      def_master
  use      def_domain
  use      def_nastal
  use      def_parame
  use      def_kermod
  use      mod_ker_space_time_function

  implicit none
  integer(ip), intent(in) :: itask !< When is it call the subroutine
  integer(ip)             :: ipoin,kpoin,ifunc,idofn
  real(rp)                :: dinew(ndofn_nsa),diold(ndofn_nsa),unfun(ndofn_nsa),unold

  if (itask == 0 .or. itask==1) then

     !  
     ! Before a time step
     !     

!     if (INOTSLAVE) then
!        if (kfl_tisch_nsa == 2) call runend('NSA_UPDBCS: NOT PROGRAMMED YET FOR BDF SCHEMES')
!     end if

     if( INOTMASTER .and. kfl_conbc_nsa==0) then
        !
        ! Space/Time functions
        !
        if( number_space_time_function > 0 ) then
           !
           ! From kernel, defined through functions
           !           
           do ipoin = 1,npoin
              if( kfl_funno_nsa(ipoin) < 0 ) then 
                 ifunc = -kfl_funno_nsa(ipoin)            
                 call ker_space_time_function(&
                      ifunc,coord(1,ipoin),coord(2,ipoin),coord(ndime,ipoin),cutim,unfun(1:ndofn_nsa))
                 do idofn = 1,ndofn_nsa
                    unfun(idofn) = unfun(idofn) * bvess_nsa(idofn,ipoin,2)
                 end do
                 if( kfl_timei_nsa /= 0 .and. kfl_tiacc_nsa == 2 .and. kfl_tisch_nsa == 1 ) then
                    do idofn = 1,ndime
                       unold = veloc(idofn,ipoin,ncomp_nsa)
                       bvess_nsa(idofn,ipoin,1) = 0.50_rp*(unfun(idofn)+unold)
                    end do
                    unold = densi(ipoin,ncomp_nsa)
                    bvess_nsa(ndime+1,ipoin,1) = 0.50_rp*(unfun(idofn)+unold)
                    unold = tempe(ipoin,ncomp_nsa)
                    bvess_nsa(ndime+2,ipoin,1) = 0.50_rp*(unfun(idofn)+unold)
                 else
                    do idofn = 1,ndofn_nsa
                       bvess_nsa(idofn,ipoin,1) = unfun(idofn)
                    end do
                 end if
              end if
           end do
        end if

        do ipoin = 1,npoin

           if(kfl_fixno_nsa(1,ipoin)==9) then
              !
              ! U(t)=U(0) - w x r
              !

              ! unprogrammed yet!

           else
              !
              ! Function: U(t)=f(t)*U(0)
              !
              ifunc= kfl_funno_nsa(ipoin)
             
              if( ifunc > 0 ) then

                 if (itask == 0) then   ! This is to update the values at ittim = 0
                     bvess_nsa(1:ndofn_nsa,ipoin,2) = bvess_nsa(1:ndofn_nsa,ipoin,1) 
                 endif

                 dinew(1:ndofn_nsa) = 0.0_rp
                 diold(1:ndofn_nsa) = bvess_nsa(1:ndofn_nsa,ipoin,2)

                 call nsa_funbou( 2_ip, ipoin,diold,dinew)

                 bvess_nsa(1:ndofn_nsa,ipoin,1) = dinew(1:ndofn_nsa)

              end if
             
           end if

        end do
        
        !   boundary elements: unprogrammed yet!
        
     end if
     
  else if (itask == 2) then
     !
     ! Before a global iteration
     !  
     if( INOTMASTER ) then
     end if

  else if (itask == 3) then
     !
     ! Before an inner iteration
     ! 
     if( INOTMASTER ) then
     end if

  end if

end subroutine nsa_updbcs
