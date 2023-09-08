subroutine exm_begite
!-----------------------------------------------------------------------
!****f* Exmedi/exm_begite
! NAME 
!    exm_begite
! DESCRIPTION
!    This routine starts the internal iteration 
! USES
!    exm_inisol
!    exm_updunk
! USED BY
!    exm_doiter
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain

  use      def_exmedi
  use mod_messages, only : livinf

  implicit none
  integer (ip) :: imate,kmodel_maxvalue
!
! Initializations
!
  kfl_goite_exm = 1
  itinn(modul)  = 0
  if(miinn_exm==0) kfl_goite_exm=0
!  if(itcou==1) call exm_tistep()
  call livinf(15_ip,' ',modul)
  !
  ! Obtain the initial guess for inner iterations: y(2) <- y(1)
  !

  !
  ! check if any TT-like model is present
  !
  kmodel_maxvalue= 0
  do imate= 1,nmate
     kmodel_maxvalue = max(kfl_cellmod(imate),kmodel_maxvalue)
  end do
  
  call exm_updunk(ITERK_EQ_ITERAUX)    ! u(,ITER_K) <-- u(,ITER_AUX)

 ! if (kmodel_maxvalue > 1) then

  !   if (kmodel_maxvalue == CELL_TT2006_EXMEDI) then  !TenTuscher 2006 Heterogeneo
        !call exm_ceihet  !!! compute CURRENTS  TT Het
  !   else if (kmodel_maxvalue ==  CELL_OHARA_EXMEDI) then  !OHara-Rudy 2011
        !call exm_oharaf  !!! compute CURRENTS O'Hara-Rudy
   !  end if

    ! call exm_upcell(ITERK_EQ_ITERAUX)  !!update the currents for the next step

  !end if


  
end subroutine exm_begite
