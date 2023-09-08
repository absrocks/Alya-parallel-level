subroutine iniste(itask)
  !-----------------------------------------------------------------------
  !****f* master/iniste
  ! NAME
  !    iniste
  ! DESCRIPTION
  !    This routine initialize variables at the beginning of a time step
  !    ITASK=1 ... Before modules have computed their time steps
  !    ITASK=2 ... After modules have computed their time steps
  ! USES
  ! USED BY
  !    Begste
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_postpr, only : ncoun_pos
  use def_coupli, only : kfl_gozon
  use def_domain, only : coord, npoin, lnods, nelem, ltype

#ifdef CATA
  use tcp
#endif

  implicit none
  integer(ip), intent(in) :: itask

  select case ( itask )

  case ( 1_ip )
     !
     ! Before modules have computed their time steps
     !
     if( kfl_timco == 1 ) dtinv = 0.0_rp

  case ( 2_ip )
     !
     ! After modules have computed their time steps 
     !
     kfl_goblk = 1
     kfl_gocou = 1
     kfl_gozon = 1
     itcou     = 1
     iblok     = 1
     iitrf     = 1 
     if( micou(iblok) == 0 ) kfl_gocou = 0
     !
     ! Postprocess counter
     !
     ncoun_pos = 0

     if ( kfl_reset < 1_ip ) then 
        ittim     = ittim + 1
        itti2     = itti2 + 1
     endif
     !
     ! call coprocessing 
     !
#ifdef CATA
     if (kfl_paral==1) call testcoprocessor(ittim,cutim,npoin,coord,nelem,lnods,ltype)
#endif
 
  end select

end subroutine iniste
