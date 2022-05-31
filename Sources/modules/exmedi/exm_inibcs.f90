subroutine exm_inibcs()
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_inibcs
  ! NAME 
  !    exm_inibcs
  ! DESCRIPTION
  !    This routine reads the boundary conditions... 
  !
  ! USES
  !    exm_bcntoe
  !    ecoute
  !    memchk
  !    runend
  ! USED BY
  !    exm_turnon
  !***
  !-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk

  use      def_exmedi

  implicit none

  if( ISLAVE ) then

     !-------------------------------------------------------------
     !
     ! Allocate memory
     !
     !-------------------------------------------------------------

     if( kfl_exboc_exm == 0 ) then
        call exm_membcs(zero)           
     else
        call exm_membcs(one)
     end if

     if( kfl_exboc_exm == 1 ) then

        !-------------------------------------------------------------
        !
        ! Node codes
        !
        !-------------------------------------------------------------

        if( kfl_icodn > 0 ) then

           iffun      =  0
           ifloc     =  0
           ifbop     =  0
           kfl_fixno => kfl_fixno_exm
           bvess     => bvess_exm
           tncod     => tncod_exm
           call reacod(10_ip)

        end if

        !-------------------------------------------------------------
        !
        ! Boundary codes
        !
        !-------------------------------------------------------------

        if( kfl_icodb > 0 ) then

           kfl_fixbo => kfl_fixbo_exm
           bvnat     => bvnat_exm(:,:,1)
           tbcod     => tbcod_exm
           call reacod(20_ip)

        end if

     end if

  end if

end subroutine exm_inibcs
