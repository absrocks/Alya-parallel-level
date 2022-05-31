module mod_bourgogne_pinotnoir

contains
  !DEC$ ATTRIBUTES FORCEINLINE :: bourgogne
  subroutine bourgogne(itask)

    use def_kintyp, only : ip,rp
    use def_master, only : INOTSLAVE,cpu_initi
    use def_parall
    implicit none
    integer(ip), intent(in) :: itask
    !CHARACTER(50)           :: s = "Alya-License.lic"
    CHARACTER(500)          :: licenseName
    INTEGER(4)              :: hourCounter = 50
    REAL(rp)                :: caca
    LOGICAL                 :: license_exist

#ifdef BOURGOGNE

    if( INOTSLAVE ) then
      !CALL get_command_argument(1, licenseName
      !licenseName = s

      CALL get_command_argument(2, licenseName)

      INQUIRE(FILE=licenseName, EXIST=license_exist)
      if(license_exist) then


        if( itask == 1 ) then
          call check(licenseName)

        else
          !Actualize License
          call cputim(caca)
          caca = ((caca - cpu_initi))
          call actualize(licenseName, caca)     
        end if
      else
        STOP "License FILE do not exist"
      end if
    end if

#endif

  end subroutine bourgogne

end module mod_bourgogne_pinotnoir
