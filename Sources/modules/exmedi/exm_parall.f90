subroutine exm_parall(itask)
!-----------------------------------------------------------------------
!****f* Exmedi/exm_parall
! NAME
!    exm_parall
! DESCRIPTION
!    This routine is a bridge to Parall service  
! USED BY
!    Exmedi
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_exmedi
  implicit none
  integer(ip), intent(in) :: itask


  call exm_chkpar

  if (weparal) then

     select case (itask)
        
     case(1)
        call exm_sendat(one)
     case(2)
        call exm_sendat(two)
     case(3)
        if(islave) then
           party =  3   ! per-node
           pardi =  1   ! vector (one dimension)
           parki =  2   ! real vector
           pard1 =  1
           parr1 => rhsid(1:npoin)
           call par_slexch()
        end if
     case(4)
        if(islave) then
           party =  3   ! per-node
           pardi =  1   ! vector (one dimension)
           parki =  2   ! real vector
           pard1 =  1
           parr1 => rhsid(1:npoin)
           call par_slexch()
        end if
     case(5)
        if(islave .or. imaster) then
           call vocabu(NPOIN_REAL_2DIM,ndime,0_ip)
           parr2 => fiber(1:ndime,1:npoin)
           call par_mygather()
           call vocabu(NPOIN_REAL_3DIM,ndime,ndime)
           parr3 => cedif_exm(1:ndime,1:ndime,1:npoin)
           call par_mygather()
           call vocabu(NELEM_INTE_1DIM,0_ip,0_ip)
           pari1 => kgrfi_exm(1:nelem)
           call par_mygather()
        end if         

     end select
     
  end if

end subroutine exm_parall
