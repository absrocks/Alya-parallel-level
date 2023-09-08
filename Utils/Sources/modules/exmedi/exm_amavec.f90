!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_amavec.f90
!> @author  Mariano Vazquez
!> @date    23/01/2017
!> @brief   Compute matrix*vector, where both are nodally defined 
!> @details Compute matrix*vector, where both are nodally defined 
!> @} 
!-----------------------------------------------------------------------
subroutine exm_amavec(amat_nodal,pp,qq,xfact_opt)

  !
  ! qq = amat_nodal pp
  !

  use def_kintyp
  use def_domain
  use def_exmedi
  use def_master
  use mod_parall, only : par_omp_npoin_chunk
  implicit none
  real(rp),    intent(in)            :: amat_nodal(*)
  real(rp),    intent(in)            :: pp(npoin)
  real(rp),    intent(out), target   :: qq(npoin)
  real(rp),    intent(in),  optional :: xfact_opt
  integer(ip)                        :: kk,ipoin,jpoin,izsym,izdom
  real(rp)                           :: raux,raux2,xfact


  if( present(xfact_opt) ) then
     xfact = xfact_opt
  else
     xfact = 1.0_rp
  end if

  if ( INOTMASTER ) then

     if( solve(1) % kfl_symme == 1 ) then
        do ipoin = 1, npoin
           raux  = 0.0_rp
           raux2 = pp(ipoin)
           do izsym= r_sym(ipoin), r_sym(ipoin+1)-2
              jpoin     = c_sym(izsym)
              raux      = raux + amat_nodal(izsym) * pp(jpoin)
              qq(jpoin) = qq(jpoin) + amat_nodal(izsym) * raux2
           end do
           kk        = r_sym(ipoin+1)-1
           qq(ipoin) = qq(ipoin) + raux + amat_nodal(kk) * pp(c_sym(kk))
        end do
     else
        do ipoin = 1,npoin
           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1    
              qq(ipoin) = qq(ipoin) + xfact * amat_nodal(izdom) * pp(c_dom(izdom))
           end do
        end do
     end if
  else
!!$     if( solve(1) % kfl_symme == 1 ) then
!!$        do ipoin = 1,npoin
!!$           raux  = 0.0_rp
!!$           raux2 = pp(ipoin)
!!$           do izsym= r_sym(ipoin), r_sym(ipoin+1)-2
!!$              jpoin     = c_sym(izsym)
!!$              raux      = raux + amat_nodal(izsym) * pp(jpoin)
!!$              qq(jpoin) = qq(jpoin) + amat_nodal(izsym) * raux2
!!$           end do
!!$           kk        = r_sym(ipoin+1)-1
!!$           qq(ipoin) = qq(ipoin) + raux + amat_nodal(kk) * pp(c_sym(kk))
!!$        end do
!!$     else
!!$        !$OMP PARALLEL DO SCHEDULE (DYNAMIC, par_omp_npoin_chunk) &
!!$        !$OMP DEFAULT  ( NONE )                                   &
!!$        !$OMP PRIVATE  ( ipoin, izdom, jpoin )                    &           
!!$        !$OMP SHARED   ( r_dom, c_dom, amat_nodal, pp, qq ,       &
!!$        !$OMP            npoin, xfact ) 
!!$        do ipoin = 1,npoin
!!$           do izdom = r_dom(ipoin),r_dom(ipoin+1)-1    
!!$              qq(ipoin) = qq(ipoin) + xfact * amat_nodal(izdom) * pp(c_dom(izdom))
!!$           end do
!!$        end do
!!$        !$OMP END PARALLEL DO
!!$     end if

  end if

end subroutine exm_amavec
