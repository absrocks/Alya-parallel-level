subroutine ada_prject
!-----------------------------------------------------------------------
!****f* adapti/ada_prject
! NAME 
!    ada_prject
! DESCRIPTION
!    This routine projects variables from the coarse to the refined mesh
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_inpout
  use      mod_postpr
  use      def_adapti
  use      mod_iofile
  implicit none
  character(5)    :: wopad(2)

  if (kfl_inter_ada == 0) return

  if (kfl_inter_ada == 2) then
     call runend('ADA_INTERP: Elsest interpolation not programmed yet.')
  else if (kfl_inter_ada == 1) then

     call ada_memall(5)

     if (kfl_chpoi_ada > 0) then        

        wopad(1)='START'        
        call ada_chkpoi(wopad)
        wopad(1)='PRESS'        
        call ada_interp(wopad)
        call ada_chkpoi(wopad)
        wopad(1)='VELOC'        
        call ada_interp(wopad)
        call ada_chkpoi(wopad)
        wopad(1)='UMOME'        
        call ada_interp(wopad)
        call ada_chkpoi(wopad)
        wopad(1)='DENSI'        
        call ada_interp(wopad)
        call ada_chkpoi(wopad)
        wopad(1)='ENERG'        
        call ada_interp(wopad)
        call ada_chkpoi(wopad)
        wopad(1)='TEMPE'        
        call ada_interp(wopad)
        call ada_chkpoi(wopad)
        wopad(1)='VISCO'        
        call ada_interp(wopad)
        call ada_chkpoi(wopad)

     end if
  
  end if


end subroutine ada_prject
