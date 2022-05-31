!-----------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_roback.f90
!> @author  Mariano Vazquez
!> @date    11/12/2013
!> @brief   Rotate back boundary conditions
!> @details Rotate back boundary conditions
!> @} 
!-----------------------------------------------------------------------
subroutine nsa_roback
  use      def_master
  use      def_domain
  use      def_parame
  use      def_kermod
  use      mod_ker_proper
  use      def_nastal
  implicit none

  integer(ip)     :: ipoin,kpoin,idime,idofn,iffix,nofix,iffit,ievat,jevat,itask,ibopo,&
       kfl_nofix(ndofn_nsa),dummy_matrix(ndofn_nsa,ndofn_nsa)

  if (INOTMASTER) then
     
     do ipoin = 1,npoin
        iffix=0
        nofix=0
        do idofn=1,ndofn_nsa
           kfl_nofix(idofn)= kfl_fixno_nsa(idofn,ipoin)
           if (kfl_nofix(idofn) .ne. 0) then
              iffix=1
              nofix=nofix + 1
           end if
        end do
        iffit=kfl_nofix(ndofn_nsa)            ! energy / temperature prescription
        !
        !when all dof's are fixed, conservative quantities are used as boundary conditions
        !
        if (nofix == ndofn_nsa) iffit= 0      

        if (kfl_nofix(1).eq.2) then
           ievat=1
        end if

        ! check if there is a boundary condition on this node
        if (iffix == 1) then
           ibopo = lpoty(ipoin)
           if( ibopo > 0 ) then
              ievat = (ipoin-1)*ndofn_nsa + 1
              jevat = (ipoin-1)*ndofn_nsa + ndofn_nsa
              call nsa_rotsys(2_ip,&
                   1_ip,1_ip,ndofn_nsa,ndofn_nsa,dummy_matrix,dunkn_nsa(ievat:jevat),&
                   jacrot_du_dq_nsa(1,1,ibopo),jacrot_dq_du_nsa(1,1,ibopo),kfl_linea_nsa,kfl_timet_nsa)
           end if
           
        end if
     end do
     
  end if

end subroutine nsa_roback
