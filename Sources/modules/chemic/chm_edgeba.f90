subroutine chm_edgeba(rhsid,unkno,amatr,pmatr)
  !-----------------------------------------------------------------------
  !****f* master/chm_edgeba
  ! NAME 
  !    chm_edgeba
  ! DESCRIPTION
  !    Diagonal solver:
  !    x^{i+1} = x^{i} + M^{-1} r
  ! USES
  !    memchk
  !    mediso
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_master, only       :  NPOIN_REAL_12DI,parr1,INOTMASTER,INOTSLAVE
  use def_master, only       :  dtinv,gesca
  use def_domain, only       :  npoin,c_dom,r_dom,c_sym,r_sym
  use def_domain, only       :  vmass
  use def_solver, only       :  solve_sol,memma,resal
  use mod_memchk
  use def_chemic
  implicit none
  real(rp),    intent(inout) :: unkno(*)
  real(rp),    intent(in)    :: amatr(1,1,*)
  real(rp),    intent(in)    :: pmatr(*)
  real(rp),    intent(inout) :: rhsid(*)
  integer(ip)                :: jnode,jpoin,itotn,jtotn,iters
  integer(ip)                :: idofn,ihang,ndofn,ipoin,izdom,jzdom,kzdom
  real(rp)                   :: fact1,xdiag,resin,resfi,dummr,dumm2
  real(rp)                   :: kij,kji,dij,dji

  call memgen(zero,npoin,zero)
  xdiag = 1.0_rp/dtinv_chm
  !
  ! IZDOM= edge i-j
  ! KZDOM= edge j-i
  !
  do ipoin = 1,npoin
     dummr = rhsid(ipoin)
     do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
        jpoin = c_dom(izdom)
        !
        ! Look for edge j-i
        !
        jzdom = r_dom(jpoin)
        do while( jzdom < r_dom(jpoin+1) )
           if( c_dom(jzdom) == ipoin ) then
              kzdom = jzdom
              jzdom = r_dom(jpoin+1)
           end if
           jzdom = jzdom + 1
        end do
        
        kij = -amatr(1,1,izdom)
        kji = -amatr(1,1,kzdom)
        dij = max(0.0_rp,-kji,-kij)
        dij = abs(kij-kji)/2.0_rp-(kij+kji)/2.0_rp

        if( ipoin /= jpoin ) then
           dummr = dummr + ( kij + dij ) * ( unkno(jpoin) - unkno(ipoin) )
        end if

     end do     

     gesca(ipoin) = xdiag / vmass(ipoin) * dummr 

  end do

  do ipoin = 1,npoin
     if( kfl_fixno_chm(1,ipoin) < 1 ) then
        unkno(ipoin) = unkno(ipoin) + gesca(ipoin)
     end if
  end do

  call memgen(two,npoin,zero)
  

end subroutine chm_edgeba
