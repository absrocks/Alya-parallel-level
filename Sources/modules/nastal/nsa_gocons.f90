subroutine nsa_gocons
!-----------------------------------------------------------------------
!****f* Nastal/nsa_gocons
! NAME
!    nsa_gocons
! DESCRIPTION
!    Advance conservative set rho, U and E
! USES
!    nsa_...
! USED BY
!    nsa_itexco
!    nsa_itexin
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_nastal
  use      def_solver
  use      mod_nsa_newelmoperations
  use      mod_nsa_elmoperations
  implicit none
  real(rp)    :: time1,time2
  !
  ! Initialization (rhsid and amatr set to zero)
  !
  call inisol()
  !
  ! Construct the system matrix and right-hand-side
  !
  call cputim(time1)

    !-----------------------------------------------------------------------!
    ! LODI
    !-----------------------------------------------------------------------!
    if(lodi_nsa==1) then
      call nsa_lodi_normaltype()
      call nsa_lodi()
    endif
    !-----------------------------------------------------------------------!

  if (INOTMASTER) then

     if (kfl_unkse_nsa < 10 ) then

        if (kfl_timet_nsa == 2 .or. kfl_matri_nsa == 1) then

           ! implicit scheme for the conservative form
           ! or explicit but with rhs evaluated as matrix * vector

           if (kfl_mod_elmop_nsa == 0)  then
              call nsa_elcons
           else
              call nsa_elmoperations
!!              call nsa_newelmoperations
           end if

        else if (kfl_timet_nsa == 1) then

           if (kfl_mod_elmop_nsa == 0)  then
              call nsa_elconsxy
           else
              call nsa_elmoperations
!!              call nsa_newelmoperations
           end if

        end if

     else if (kfl_unkse_nsa < 20 ) then

        call nsa_elchea

     end if

     ! DMM the call for nsa_bouope should be here

  end if
  !
  ! Timing
  !
  call cputim(time2)
  cpu_modul(CPU_ASSEMBLY,modul) = cpu_modul(CPU_ASSEMBLY,modul) + time2 - time1
  !
  ! Advance explicitly or implicitly
  !

!!!  write(6,*) 'totoro',rhsid(1:4)

  call nsa_upcons

!  write(6,*) rhsid(1:10)
!  write(6,*) dunkn_nsa(1:10)
!  write(6,*) unkno(1:10)
!  stop


end subroutine nsa_gocons
