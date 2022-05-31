!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_upcell.f90
!> @author  Mariano Vazquez
!> @date    28/11/2016
!> @brief   Update cell model variables
!> @details Update cell model variables
!> @} 
!-----------------------------------------------------------------------
subroutine exm_upcell(itask)
  use      def_parame
  use      def_master
  use      def_domain
  use      def_exmedi

  implicit none
  integer(ip), intent(in) :: itask !> case to update
  integer(ip)             :: ipoin,itime,kmodel_ipoin

  if( INOTMASTER ) then
     
     select case ( itask )
 
     case( ITASK_BEGSTE )
        !
        ! (:,2) <= (:,1)
        !
        do ipoin = 1,npoin
           kmodel_ipoin = kfl_cellmod(nodemat(ipoin))
           if (kmodel_ipoin==CELL_FITZHUGH_EXMEDI) then
              refhn_exm(ipoin,2)= refhn_exm(ipoin,1)
           else if (kmodel_ipoin==CELL_TT2006_EXMEDI) then !TenTuscher 2006 Heterogeneo
              !do idofn=1,nicel_exm   
              !   vicel_exm(idofn,ipoin,ITER_AUX) = vicel_exm(idofn,ipoin,TIME_N) 
              !end do
           end if
        end do
     
     case ( ITASK_ENDITE )
        !
        ! (:,2) <= (:,1)
        !        
        do ipoin = 1,npoin
           kmodel_ipoin = kfl_cellmod(nodemat(ipoin))
           if (kmodel_ipoin==CELL_FITZHUGH_EXMEDI) then   ! fhn
              refhn_exm(ipoin,2) =refhn_exm(ipoin,1)
           end if
        end do

     case( ITASK_ENDSTE )
        !
        ! (:,3) <= (:,1)
        ! (:,4) <= (:,3)...
        !
        if( kfl_tiacc_exm == 2 ) then
           do itime=2+kfl_tiacc_exm,4,-1
              do ipoin=1,npoin                 
                 kmodel_ipoin = kfl_cellmod(nodemat(ipoin))
                 if (kmodel_ipoin==CELL_FITZHUGH_EXMEDI) then   ! fhn
                    refhn_exm(ipoin,itime)=refhn_exm(ipoin,itime-1)                 
                 end if
              end do
           end do
        end if
        
        do ipoin = 1,npoin
           kmodel_ipoin = kfl_cellmod(nodemat(ipoin))
           if (kmodel_ipoin==CELL_FITZHUGH_EXMEDI) then   ! fhn
              refhn_exm(ipoin,TIME_N)=refhn_exm(ipoin,ITER_K)
           end if
        end do
        
     end select

  end if

end subroutine exm_upcell
