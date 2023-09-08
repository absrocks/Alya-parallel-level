subroutine symgra()
  !-----------------------------------------------------------------------
  !****f* Domain/symgra
  ! NAME
  !    symgra
  ! DESCRIPTION
  !    This routine generates the symmetric graph
  ! OUTPUT
  !    NZSYM ... Number of nonzero coefficients of the graph
  !    R_SYM ... Pointer to the array of rows r_dom(npoin+1) (r_dom(ipoin) = 
  !              coefficient of the graph where row ipoin starts)
  !    C_SYM ... Pointer to the array of columns c_dom(nzdom) (c_dom (izdom)
  !              = column of the izdom coefficient of mesh graph)
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use mod_memchk
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: izsym,izdom,ipoin,jpoin,izsta,izlen
  integer(4)  :: istat
  !
  ! Symmetric graph is needed if there exist a solver with symmetric graph
  !
  call parari('BCT',0_ip,1_ip,kfl_symgr)

  if( kfl_symgr == 1 ) then

     if( INOTMASTER ) then

        kfl_symgr = 2
        call livinf(62_ip,' ',0_ip)
        allocate(r_sym(npoin+1),stat=istat)
        call memchk(zero,istat,memor_dom,'R_SYM','symgra',r_sym)
        allocate(c_sym(nzsym),stat=istat)
        call memchk(zero,istat,memor_dom,'C_SYM','symgra',c_sym)
        izsym = 0
        izsta = 1
        do ipoin = 1,npoin
           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
              jpoin = c_dom(izdom)
              if( jpoin <= ipoin ) then
                 izsym        = izsym + 1
                 c_sym(izsym) = jpoin
              end if
           end do
           r_sym(ipoin) = izsta
           izlen        = izsym - izsta + 1
           call heapsorti1(2_ip,izlen,c_sym(izsta))
           izsta = izsym + 1
        end do
        r_sym(npoin+1) = izsym + 1

     else

        allocate(r_sym(1),stat=istat)
        call memchk(zero,istat,memor_dom,'R_SYM','symgra',r_sym)
        allocate(c_sym(1),stat=istat)
        call memchk(zero,istat,memor_dom,'C_SYM','symgra',c_sym)

     end if

  end if

end subroutine symgra
