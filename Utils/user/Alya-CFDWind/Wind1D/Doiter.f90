subroutine doiter
  !-----------------------------------------------------------------------
  !****f* doiter
  ! NAME 
  !    doiter
  ! DESCRIPTION
  !    This routine solves an iteration of the linearized
  !    equations.
  ! USES
  !   solite  !   
  ! USED BY
  !   Wind
  !***
  !-----------------------------------------------------------------------

  use      def_master
  implicit none
  integer(4) iterb, iunkn

  iterb =0  
  kfl_goite=1
  do while(kfl_goite==1.and.iterb.lt.maxit(1)) !Inner iteration
     iterb = iterb +1 
     do iunkn=1, 2   ! velocity block         
        call solite(iunkn) ! Solves a linear iteration for unkno iunkn
     end do

  end do
  iterb =0  
  kfl_goite=1
  do while(iterb.lt.maxit(3)) !Inner iteration
     iterb = iterb +1 
     do iunkn=3, 4   ! k-epsilon block         
        call solite(iunkn) ! Solves a linear iteration for unkno iunkn      
     end do
  end do
  
  ! If thermal transient, solve temper equation.
  if (kfl_thmod.eq.1)   call solite(5)

end subroutine doiter
