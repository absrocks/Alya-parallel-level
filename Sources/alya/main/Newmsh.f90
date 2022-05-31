subroutine Newmsh(order)
!-----------------------------------------------------------------------
!****f* master/Newmsh
! NAME
!    Newmsh
! DESCRIPTION
!    This routine defines a new mesh
! USES
!    updmsh
!    defmsh
!    promsh
! USED BY
!    Alya
!***
!-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use mod_moduls, only : moduls
  use mod_module_interface
  implicit none
  integer(ip), intent(in)  :: order

  if (kfl_paral==-1) then

     if(kfl_algor_msh==0) then
        return
     end if

     select case (order)

     case(zero)


     case(one)

        if(kfl_algor_msh==1) then

           if(iblok == 2) then    ! Only after the ALE module
           else
              return
           end if
           !
           ! Deallocate old mesh vectors
           !
           call updmsh(one)
           !
           ! Define a new mesh
           !
           !call defmsh
           !
           ! Project data from old to new mesh for each field
           !
           call Nastin(five)
           !
           ! Stop the SSL process (No iterations)
           !
           !iblok = iblok + 1
        else if(kfl_algor_msh==2) then
           !
           ! Save old mesh info and deallocate vectors
           !
           call updmsh(one)
           !
           ! Define a new mesh
           !
           !call defmsh
           !
           ! Project data from old to new mesh for each field
           !
           call moduls(ITASK_NEWMSH)
           !
           ! Deallocate old mesh vectors
           !
           call updmsh(two)
           !
           ! Check if the refinement loop has to be stopped or not.
           !
           iitrf = iitrf + 1
           if(iitrf>mitrf) iblok = iblok + 1        ! End refinement loop

        end if

     end  select

  end if

end subroutine Newmsh
