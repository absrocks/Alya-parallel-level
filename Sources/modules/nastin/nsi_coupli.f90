!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_coupli.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   Nastal coupling 
!> @details Nastal coupling
!> @} 
!-----------------------------------------------------------------------
subroutine nsi_coupli(itask)
  use def_master
  use def_domain
  use def_elmtyp
  use def_nastin
  use def_coupli,      only : FIXED_UNKNOWN
  use mod_nsi_commdom, only : nsi_commdom_lm2_code_i, nsi_commdom_lm2_code_j
  use mod_local_basis, only : local_basis_global_to_local
  use mod_local_basis, only : local_basis_local_to_global
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ipoin,idime,idofn,incnt
  integer(ip)             :: ipoin_fluid,ipoin_solid

  select case ( itask )

  case ( ITASK_CONCOU )
     call nsi_commdom_lm2_code_i(-1_ip, -1_ip)
     call nsi_commdom_lm2_code_j(-1_ip, -1_ip)

  case ( ITASK_BEGSTE )
     call nsi_commdom_lm2_code_i(-1_ip, -1_ip)
     call nsi_commdom_lm2_code_j(-1_ip, -1_ip)


  case ( ITASK_INIUNK )
     call nsi_commdom_lm2_code_i(-1_ip, -1_ip)
     call nsi_commdom_lm2_code_j(-1_ip, -1_ip)

     !-------------------------------------------------------------------
     !
     ! INIUNK
     !
     !-------------------------------------------------------------------

  case ( ITASK_BEGITE )

     !-------------------------------------------------------------------
     !
     ! BEGITE
     !
     !-------------------------------------------------------------------

     if( INOTMASTER ) then       
        !
        ! Coupling with ALEFOR: ALE + boundary condition
        !
        if( coupling('NASTIN','ALEFOR') >= 1 ) then
           if( kfl_conbc_nsi /= 1 ) then
              kfl_advec = 2
              call local_basis_global_to_local(kfl_fixrs_nsi,velom)           ! Global to local     
              do ipoin = 1,npoin
                 if( kfl_funno_nsi(ipoin) == 0 ) then ! to not overwrite the bc from a function
                    do idime = 1,ndime
                       if( kfl_fixno_ale(idime,ipoin) == 1 .or. kfl_fixno_ale(idime,ipoin) == FIXED_UNKNOWN ) then
                         bvess_nsi(idime,ipoin,1) = bvess_nsi(idime,ipoin,2) + velom(idime,ipoin)
                       end if
                    end do
                 end if
              end do
              call local_basis_local_to_global(kfl_fixrs_nsi,velom)           ! Local to global
           end if
        end if
     end if

  case ( ITASK_ENDITE )

     !-------------------------------------------------------------------
     !
     ! ENDITE
     !
     !-------------------------------------------------------------------
     !
     ! Coupling with Immbou: compute force
     !
     if( coupling('IMMBOU','NASTIN') >= 1 ) then
        call nsi_immbou()
     end if
     !
     ! Coupling with Alefor: compute force for rigid body
     !
     if( ( coupling('ALEFOR','NASTIN') >= 1 ) .and. ( nrbod > 0_ip ) ) then
        call nsi_rbobou()
     end if

  case ( ITASK_MATRIX )

     !-------------------------------------------------------------------
     !
     ! MATRIX: After assembly
     !
     !-------------------------------------------------------------------

     if( coupling('NASTIN','IMMBOU') == 1 ) then
        !
        ! Coupling with Immbou: mass matrix conservation variables
        !                
        call nsi_massma()
        if( INOTMASTER ) then
           !
           ! Coupling with Immbou: impose force FORCF: interpolate b.c.
           !     
           call nsi_embedd(&
                amatr(poauu_nsi),amatr(poaup_nsi),amatr(poapu_nsi),amatr(poapp_nsi),&
                lapla_nsi,rhsid,rhsid(ndbgs_nsi+1),unkno,unkno(ndbgs_nsi+1))           
        end if
     end if
     !
     ! Coupling with Partis: take off momentum from momentum equations
     !    
     call nsi_partis()

  end select

!===============================================================| contains |===!
contains
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_nsi_set_displacement( dummy, sld_displacement )
  use      def_master, only: bvess_ale
  implicit none
  real(rp), intent(out)   :: dummy(:,:)
  real(rp), intent(inout) :: sld_displacement(:,:)
  !
  if(INOTMASTER) then
      bvess_ale(1:ndime,1:,1) = sld_displacement(1:ndime,1:)
  endif
  !
  dummy(1:ndime,1:) = 0.0_rp
  sld_displacement(1:ndime,1:) = 0.0_rp
  !
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_nsi_get_force( force_in, dummy, i_fixbo)
  use def_domain, only: ndime,r_dom,c_dom,r_sol,c_sol
  use def_domain, only: kfl_codno
  use def_nastin, only: intfo_nsi
  use def_master, only: veloc, press
  implicit none
  real(rp), intent(out) :: force_in(:,:)
  real(rp), intent(out) :: dummy(:,:)
  integer(ip), intent(in) :: i_fixbo 

  real(rp) :: foref(ndime)
  integer(ip) :: jzdom, jpoin, ipoin, idime, jdime, izdom

  if ( INOTMASTER ) then
        force_in(1:ndime,1:npoin) = 0_rp
        !
        do ipoin = 1_ip,npoin
          !
          foref(1:ndime) = 0.0_rp
          if(kfl_codno(1,ipoin) == i_fixbo) then
              foref(1:ndime) = intfo_nsi(ipoin) % bu(1:ndime)
              !
              jzdom = 0
              do izdom = r_dom(ipoin),r_dom(ipoin+1) - 1
                jpoin  = c_dom(izdom)
                jzdom  = jzdom + 1
                !
                foref(1:ndime) = foref(1:ndime) - intfo_nsi(ipoin) % Aup(1:ndime,jzdom) * press(jpoin,1)
                foref(1:ndime) = foref(1:ndime) &
                                      - matmul( intfo_nsi(ipoin) % Auu(1:ndime,1:ndime,jzdom), &
                                                veloc(1:ndime,jpoin,1) )
                !
                !do idime = 1,ndime
                !  foref(idime) = foref(idime) - dot_product( intfo_nsi(ipoin) % Auu(1:ndime,idime,jzdom), &
                !                                 veloc(1:ndime,jpoin,1) )
                !enddo
                !
              end do
              !
          endif
          !
          force_in(1:ndime,ipoin) = force_in(1:ndime,ipoin) + foref(1:ndime)
          !
        end do
        !
        call rhsmod(ndime, force_in(1:ndime,1:npoin) )
        !
  endif
  !
  dummy(1:ndime,1:npoin) = 0.0_rp
  ! 
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

!===============================================================| contains |===!
end subroutine nsi_coupli
