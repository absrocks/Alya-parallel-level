subroutine matspr(i_row,j_col,nzso2)
  !-----------------------------------------------------------------------
  !****f* solmum/matspr
  ! NAME
  !    matspr
  ! DESCRIPTION
  !    Need when ndofn_sol>1 to obtain the correct values of i_row and j_col
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only    :  ip
  use def_domain, only    :  npoin,nzsol,r_sol,c_sol
  use def_solver, only    :  solve_sol
  implicit none
  integer(4), intent(in)  :: nzso2
  integer(4), intent(out) :: i_row(nzso2),j_col(nzso2)
  integer(4)              :: irs,ire,ipoin,jpoin,idofn,jdofn,i,izdom
  integer(4)              :: ndofn_sol4

  if(ip==4) then
     i=0
     do ipoin=1,npoin
        irs=r_sol(ipoin)
        ire=r_sol(ipoin+1)-1
        do izdom=irs,ire
           jpoin=c_sol(izdom) 
           do idofn=1,solve_sol(1)%ndofn
              do jdofn=1,solve_sol(1)%ndofn
                 i=i+1
                 i_row(i)=(ipoin-1)*solve_sol(1)%ndofn+idofn
                 j_col(i)=(jpoin-1)*solve_sol(1)%ndofn+jdofn
              end do
           end do
        end do
     end do
  else
     i=0
     ndofn_sol4=int(solve_sol(1)%ndofn,4)
     do ipoin=1,npoin
        irs=int(r_sol(ipoin),4)
        ire=int(r_sol(ipoin+1),4)-1
        do izdom=irs,ire
           jpoin=int(c_sol(izdom),4)
           do idofn=1,ndofn_sol4
              do jdofn=1,ndofn_sol4
                 i=i+1
                 i_row(i)=(ipoin-1)*ndofn_sol4+idofn
                 j_col(i)=(jpoin-1)*ndofn_sol4+jdofn
              end do
           end do
        end do
     end do
  end if

end subroutine matspr
