subroutine ale_mealrb(itask)
  !-----------------------------------------------------------------------
  !****f* alefor/ale_mealrb
  ! NAME
  !    ale_mealrb
  ! DESCRIPTION
  !    Similar to ibm_memall but for rigid body in ale
  !    For the moment only case 2 & 3 are ready
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
  use def_alefor
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iimbo,ppoib,pboib
  integer(4)              :: istat

  select case ( itask )

  case ( 0_ip )

  case ( 1_ip )
     !
     ! Immersed boundaries (IB) 
     !
     call runend('ale_mealrb: it should not enter here')  ! what ibm_memall has in case 1 is for non body fitted

  case ( 2_ip )
     !
     ! IB type
     !
     if( nrbod > 0 .and. kfl_rigid_ale ==1 ) then
        allocate(rbbou(nrbod),stat=istat)
     end if

  case ( 3_ip )
     !
     ! Immersed boundaries (IB) of body fitted type
     !
     if( kfl_rigid_ale ==1 ) then    
        iimbo = igene
        ppoib = rbbou(iimbo) % npoib
        pboib = rbbou(iimbo) % nboib
        allocate(rbbou(iimbo) % cooin(ndime,ppoib),stat=istat)
        call memchk(zero,istat,memor_dom,'rbbou(iimbo) % cooin','ale_mealrb',rbbou(iimbo) % cooin) ! COOIN
        allocate(rbbou(iimbo) % cooib(ndime,ppoib),stat=istat)
        call memchk(zero,istat,memor_dom,'rbbou(iimbo) % cooib','ale_mealrb',rbbou(iimbo) % cooib) ! COOIB
        !     allocate(rbbou(iimbo) % cooi2(ndime,ppoib),stat=istat)
        !     call memchk(zero,istat,memor_dom,'rbbou(iimbo) % cooi2','ale_mealrb',rbbou(iimbo) % cooi2) ! COOI2
        allocate(rbbou(iimbo) % lnoib(mnoib,pboib),stat=istat)
        call memchk(zero,istat,memor_dom,'rbbou(iimbo) % lnoib','ale_mealrb',rbbou(iimbo) % lnoib) ! LNOIB
        allocate(rbbou(iimbo) % ltyib(pboib),stat=istat)
        call memchk(zero,istat,memor_dom,'rbbou(iimbo) % ltyib','ale_mealrb',rbbou(iimbo) % ltyib) ! LTYIB
        allocate(rbbou(iimbo) % lninv(ppoib),stat=istat)
        call memchk(zero,istat,memor_dom,'rbbou(iimbo) % lninv','ale_mealrb',rbbou(iimbo) % lninv) ! LNINV
        allocate(rbbou(iimbo) % lbinv(pboib),stat=istat)
        call memchk(zero,istat,memor_dom,'rbbou(iimbo) % lbinv','ale_mealrb',rbbou(iimbo) % lbinv) ! LBINV
     end if

  end select

end subroutine ale_mealrb
