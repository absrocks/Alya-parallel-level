!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_output.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Output 
!> @details Output 
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_output()
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal
  use      mod_postpr
  use      mod_matrix, only : matrix_output_gid_format

  implicit none
  integer(ip)             :: ivarp,ivari,ivort

  ivort = 0
  kfl_crens_nsa= 0
  do ivarp = 1,nvarp
     ivari = ivarp
     call posdef(11_ip,ivari)
     if( ivort == 0 .and. ( ivari == 4 .or. ivari == 29 .or. ivari == 31 .or. ivari == 32 ) ) then
        call vortic(-1_ip)
        ivort = 1
     end if
     call nsa_outvar(ivari)
  end do

  if( postp(1) % npp_stepi(55) /= 0 ) then
     if( mod(ittim,  postp(1) % npp_stepi(55) ) == 0 ) then ! AVVEL frequency
        avtim_nsa = cutim  ! Update reference time for time-averaging
     endif
  endif
 
  if( ivort == 1 ) call vortic(3_ip)

  if( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! Calculations on sets
     !
     call nsa_outset()

     !
     ! Postprocess on witness points
     !
     call nsa_outwit()            
     

!   else if( ittyp == ITASK_ENDRUN ) then

!      if( INOTMASTER ) then
!         call pspltm(&
!              npoin, npoin,ndime+2,0_ip,c_dom,r_dom,amatr,&
!              trim(title)//': '//namod(modul),0_ip,18.0_rp,'cm',&
!              0_ip,0_ip,2_ip,90+kfl_paral)
!     else
!      call matrix_output_gid_format(npoin,ndime+2_ip,r_dom,c_dom,amatr,momod(modul) % solve(1) % invpr_gs)
!     end if
!  end if

  !...special stuff:     
  !
  !   if (ivari==36 .and. ndime==2) then
  !      if (kfl_pro2d_nsa == 1) call nsa_outpro(two)
  !      call nsa_outpro(three)
     end if

  if (kfl_chkpo_nsa(2) < 0) kfl_chkpo_nsa(2) = - kfl_chkpo_nsa(2)
 

end subroutine nsa_output
