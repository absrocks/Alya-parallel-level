subroutine ada_coamsh(itask)
!-----------------------------------------------------------------------
!****f* adapti/ada_divmsh
! NAME 
!    ada_divmsh
! DESCRIPTION
!    This routine performs adaptivity mesh division 
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_elmtyp

  use      def_adapti
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: &
       ielem,iechi,ieori,imark,imori

 
  if (kfl_redgr_ada >= 1) then

     if (itask == ITASK_ENDSTE) then
        
        call ada_livinf('Red coarsening -> Proceed.')

        do ieori=1,neori_ada
           if (lsori_ada(ieori)>0) then
              ! redened element: check if it should be coarsened
              imori= 0
              do iechi= 1,nechi_ada
                 ielem= lfaso_ada(iechi,ieori)
                 imark = 1           
                 call ada_markel(ielem,imark)                 
                 imori= imori + imark
              end do              
           end if

           if (imori == 0) then
              ! coarse this cell: children's error estimator below refinement threshold
              
           end if

        end do
        
        call ada_livinf('Red coarsening -> Done.')

     end if

  end if


end subroutine ada_coamsh
