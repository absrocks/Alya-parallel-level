subroutine ibm_memall(itask)
  !-----------------------------------------------------------------------
  !****f* immbou/ibm_memall
  ! NAME
  !    ibm_memall
  ! DESCRIPTION
  !    Allocate/Deallocate the geometry arrays 
  !    ITASK=1 ... Allocate memory
  !    ITASK=2 ... Deallocate memory
  ! OUTPUT
  ! USED BY
  !    reageo
  !    sengeo
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_parame
  use def_master
  use def_domain
  use def_inpout
  use mod_memchk
  use def_immbou
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iimbo,ppoib,pboib,pelem,ppoin,nsize
  integer(4)              :: istat

  select case ( itask )

  case ( 0_ip )

  case ( 1_ip )
     !
     ! Immersed boundaries (IB) 
     !
     do iimbo = 1,nimbo

        if( imbou(iimbo) % kfl_typeb == 0 .or. imbou(iimbo) % kfl_typeb == -2 ) then
           !
           ! Embedded / body force type
           !
           ppoib = imbou(iimbo) % npoib
           pboib = imbou(iimbo) % nboib
           allocate(imbou(iimbo) % lnoib(mnoib,pboib),stat=istat)
           call memchk(zero,istat,memor_dom,'imbou(iimbo) % lnoib','ibm_memall',imbou(iimbo) % lnoib) ! LNOIB
           allocate(imbou(iimbo) % ltyib(pboib),stat=istat)
           call memchk(zero,istat,memor_dom,'imbou(iimbo) % ltyib','ibm_memall',imbou(iimbo) % ltyib) ! LTYIB
           allocate(imbou(iimbo) % cooib(ndime,ppoib),stat=istat)
           call memchk(zero,istat,memor_dom,'imbou(iimbo) % cooib','ibm_memall',imbou(iimbo) % cooib) ! COOIB
           allocate(imbou(iimbo) % cooin(ndime,ppoib),stat=istat)
           call memchk(zero,istat,memor_dom,'imbou(iimbo) % cooin','ibm_memall',imbou(iimbo) % cooin) ! COOIN
           allocate(imbou(iimbo) % cooi2(ndime,ppoib),stat=istat)
           call memchk(zero,istat,memor_dom,'imbou(iimbo) % cooi2','ibm_memall',imbou(iimbo) % cooi2) ! COOI2

        else if( imbou(iimbo) % kfl_typeb == -1 ) then
           !
           ! Volume type
           !
           ppoib = imbou(iimbo) % npoib
           pboib = imbou(iimbo) % nboib
           ppoin = imbou(iimbo) % npoin
           pelem = imbou(iimbo) % nelem
           allocate(imbou(iimbo) % lnoib(mnoib,pboib),stat=istat)
           call memchk(zero,istat,memor_dom,'imbou(iimbo) % lnoib','ibm_memall',imbou(iimbo) % lnoib) ! LNOIB
           allocate(imbou(iimbo) % ltyib(pboib),stat=istat)
           call memchk(zero,istat,memor_dom,'imbou(iimbo) % ltyib','ibm_memall',imbou(iimbo) % ltyib) ! LTYIB
           allocate(imbou(iimbo) % lnods(mnodi,pelem),stat=istat)
           call memchk(zero,istat,memor_dom,'imbou(iimbo) % lnods','ibm_memall',imbou(iimbo) % lnods) ! LNODS
           allocate(imbou(iimbo) % ltype(pelem),stat=istat)
           call memchk(zero,istat,memor_dom,'imbou(iimbo) % ltype','ibm_memall',imbou(iimbo) % ltype) ! LTYPE
           allocate(imbou(iimbo) % coord(ndime,ppoin),stat=istat)
           call memchk(zero,istat,memor_dom,'imbou(iimbo) % coord','ibm_memall',imbou(iimbo) % coord) ! COORD
           imbou(iimbo) % cooib => imbou(iimbo) % coord
           !imbou(iimbo) % cooin => imbou(iimbo) % coord
           !imbou(iimbo) % cooi2 => imbou(iimbo) % coord
           allocate(imbou(iimbo) % cooin(ndime,ppoib),stat=istat)
           call memchk(zero,istat,memor_dom,'imbou(iimbo) % cooin','ibm_memall',imbou(iimbo) % cooin) ! COOIN
           !allocate(imbou(iimbo) % cooib(ndime,ppoib),stat=istat)
           !call memchk(zero,istat,memor_dom,'imbou(iimbo) % cooib','ibm_memall',imbou(iimbo) % cooib) ! COOIB
           allocate(imbou(iimbo) % cooi2(ndime,ppoib),stat=istat)
           call memchk(zero,istat,memor_dom,'imbou(iimbo) % cooi2','ibm_memall',imbou(iimbo) % cooi2) ! COOI2

        end if
     end do

  case ( 2_ip )
     !
     ! IB type
     !
     if( nimbo > 0 ) then
        allocate(imbou(nimbo),stat=istat)
     end if

  case ( 3_ip )
     !
     ! Immersed boundaries (IB) of body fitted type
     !    
     iimbo = igene
     ppoib = imbou(iimbo) % npoib
     pboib = imbou(iimbo) % nboib
     allocate(imbou(iimbo) % cooin(ndime,ppoib),stat=istat)
     call memchk(zero,istat,memor_dom,'imbou(iimbo) % cooin','ibm_memall',imbou(iimbo) % cooin) ! COOIN
     allocate(imbou(iimbo) % cooib(ndime,ppoib),stat=istat)
     call memchk(zero,istat,memor_dom,'imbou(iimbo) % cooib','ibm_memall',imbou(iimbo) % cooib) ! COOIB
     allocate(imbou(iimbo) % cooi2(ndime,ppoib),stat=istat)
     call memchk(zero,istat,memor_dom,'imbou(iimbo) % cooi2','ibm_memall',imbou(iimbo) % cooi2) ! COOI2
     allocate(imbou(iimbo) % lnoib(mnoib,pboib),stat=istat)
     call memchk(zero,istat,memor_dom,'imbou(iimbo) % lnoib','ibm_memall',imbou(iimbo) % lnoib) ! LNOIB
     allocate(imbou(iimbo) % ltyib(pboib),stat=istat)
     call memchk(zero,istat,memor_dom,'imbou(iimbo) % ltyib','ibm_memall',imbou(iimbo) % ltyib) ! LTYIB
     allocate(imbou(iimbo) % lninv(ppoib),stat=istat)
     call memchk(zero,istat,memor_dom,'imbou(iimbo) % lninv','ibm_memall',imbou(iimbo) % lninv) ! LNINV
     allocate(imbou(iimbo) % lbinv(pboib),stat=istat)
     call memchk(zero,istat,memor_dom,'imbou(iimbo) % lbinv','ibm_memall',imbou(iimbo) % lbinv) ! LBINV

  case ( 4_ip )
     !
     ! Model type (IB): Int. and real parameters: IPARA and RPARA
     !  
     iimbo = igene
     if( imbou(iimbo) % kfl_model /= 0 ) then
        nsize = imbou(iimbo) % npari
        allocate(imbou(iimbo) % ipara(nsize),stat=istat)
        call memchk(zero,istat,memor_dom,'imbou(iimbo) % ipara','ibm_memall',imbou(iimbo) % ipara) ! IPARA
        nsize = imbou(iimbo) % nparr
        allocate(imbou(iimbo) % rpara(nsize),stat=istat)
        call memchk(zero,istat,memor_dom,'imbou(iimbo) % rpara','ibm_memall',imbou(iimbo) % rpara) ! RPARA        
     end if

  end select

end subroutine ibm_memall
