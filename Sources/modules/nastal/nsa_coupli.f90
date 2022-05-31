!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    nsa_coupli.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   This routine manages nastal coupling with IMMBOU, ALEFOR, SOLIDZ, TEMPER...
!> @details This routine manages nastal coupling with IMMBOU, ALEFOR, SOLIDZ, TEMPER...
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_coupli(itask)
  use def_master
  use def_elmtyp
  use def_domain
  use def_nastal
  use mod_nsa_commdom, only: nsa_commdom_lmc_code_i, nsa_commdom_lmc_code_j
  implicit none
  integer(ip), intent(in) :: itask   !< who is calling nsa_coupli
  integer(ip)             :: ipoin,kpoin,idime
  integer(ip)             :: ipoin_fluid,ipoin_solid
  real(rp)                :: verot(3)

  select case ( itask )

  case ( ITASK_CONCOU )
     call nsa_commdom_lmc_code_i(-1_ip, -1_ip)
     call nsa_commdom_lmc_code_j(-1_ip, -1_ip)

  case ( ITASK_BEGSTE )
     call nsa_commdom_lmc_code_i(-1_ip, -1_ip)
     call nsa_commdom_lmc_code_j(-1_ip, -1_ip)


  case ( ITASK_INIUNK )
     call nsa_commdom_lmc_code_i(-1_ip, -1_ip)
     call nsa_commdom_lmc_code_j(-1_ip, -1_ip)
     !
     ! Initialize
     !
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
        if( coupling('NASTAL','ALEFOR') >= 1 ) then
           kfl_advec = 2 
           do ipoin = 1,npoin
              verot(1:ndime)= velom(1:ndime,ipoin)
              call nsa_rotunk(-1,ipoin,ipoin,verot) ! global to local
              do idime = 1,ndime
                 if( kfl_fixno_ale(idime,ipoin) == 1 ) then
                    bvess_nsa(idime,ipoin,1) = bvess_nsa(idime,ipoin,2) + verot(idime)
                 end if
              end do
           end do
        end if

     end if
     ! recalculates initial level set & hydrostatic pressure on the deformed mesh 
!!     if(( kfl_colev_nsi /= 0 ).and.( coupling('NASTIN','ALEFOR') >= 1 )) call nsi_rechyd()


  case ( ITASK_ENDITE )

     !-------------------------------------------------------------------
     !
     ! ENDITE
     !
     !-------------------------------------------------------------------
     !
     ! Coupling with Immbou: compute force
     !
!     if( coupling('IMMBOU','NASTAL') >= 1 ) then
!        call nsa_immbou()   ! TO BE DEFINED
!     end if
     !
     ! Coupling with Alefor: compute force for rigid body
     !
!     if( ( coupling('ALEFOR','NASTIN') >= 1 ) .and. ( nrbod > 0_ip ) ) then
!        call nsa_rbobou()   ! TO BE DEFINED
!     end if
     !
     ! Coupling with solidz: compute force
     !
     if( coupling('SOLIDZ','NASTAL') >= 1 ) then
!        call nsa_solidz()   ! TO BE DEFINED
     end if

  case ( ITASK_MATRIX )

     !-------------------------------------------------------------------
     !
     ! MATRIX: After assembly
     !
     !-------------------------------------------------------------------

     if( coupling('NASTAL','IMMBOU') == 1 ) then      ! TO BE DEFINED
        !
        ! Coupling with Immbou: mass matrix conservation variables
        !                
!!        call nsi_massma()
!!        if( INOTMASTER ) then
!!           !
!!           ! Coupling with Immbou: impose force FORCF: interpolate b.c.
!!           !     
!!           call nsi_embedd(&
!!                amatr(poauu_nsi),amatr(poaup_nsi),amatr(poapu_nsi),amatr(poapp_nsi),&
!!                lapla_nsi,rhsid,rhsid(ndbgs_nsi+1),unkno,unkno(ndbgs_nsi+1))           
!!        end if
     end if

     if( coupling('NASTAL','PARTIS') == 1 .and. INOTMASTER ) then        ! TO BE DEFINED
        !
        ! Coupling with Partis: take off momentum from momentum equations
        !    
!!        call nsi_partis()
     end if


  end select
end subroutine nsa_coupli
