module mod_module_interface

  use def_kintyp,          only : ip,rp
  implicit none
  public

contains

#ifdef CMAKE

  subroutine nastin(itask)
#ifdef NASTIN_MODULE
    NASTIN_MODULE mod_nastin
#endif
    integer(ip), intent(in) :: itask
#ifdef NASTIN_MODULE
    call nastin_main(itask)
#else
    call runend("NASTIN IS NOT AVAILABLE! PLEASE COMPILE THE NASTIN LIBRARY")
#endif
  end subroutine

  subroutine temper(itask)
#ifdef TEMPER_MODULE
    TEMPER_MODULE mod_temper
#endif
    integer(ip), intent(in) :: itask
#ifdef TEMPER_MODULE
    call temper_main(itask)
#else
    call runend("TEMPER IS NOT AVAILABLE! PLEASE COMPILE THE TEMPER LIBRARY")
#endif
  end subroutine
  
  subroutine partis(itask)
#ifdef PARTIS_MODULE
    PARTIS_MODULE mod_partis
#endif
    integer(ip), intent(in) :: itask
#ifdef PARTIS_MODULE
    call partis_main(itask)
#else
    call runend("PARTIS IS NOT AVAILABLE! PLEASE COMPILE THE PARTIS LIBRARY")
#endif
  end subroutine

  subroutine alefor(itask)
#ifdef ALEFOR_MODULE
    ALEFOR_MODULE mod_alefor
#endif
    integer(ip), intent(in) :: itask
#ifdef ALEFOR_MODULE
    call alefor_main(itask)
#else
    call runend("ALEFOR IS NOT AVAILABLE! PLEASE COMPILE THE ALEFOR LIBRARY")
#endif
  end subroutine

  subroutine solidz(itask)
#ifdef SOLIDZ_MODULE
    SOLIDZ_MODULE mod_solidz
#endif
    integer(ip), intent(in) :: itask
#ifdef SOLIDZ_MODULE
    call solidz_main(itask)
#else
    call runend("SOLIDZ IS NOT AVAILABLE! PLEASE COMPILE THE SOLIDZ LIBRARY")
#endif
  end subroutine
  
  subroutine exmedi(itask)
#ifdef EXMEDI_MODULE
    EXMEDI_MODULE mod_exmedi
#endif
    integer(ip), intent(in) :: itask
#ifdef EXMEDI_MODULE
    call exmedi_main(itask)
#else
    call runend("EXMEDI IS NOT AVAILABLE! PLEASE COMPILE THE EXMEDI LIBRARY")
#endif
  end subroutine

  subroutine chemic(itask)
#ifdef CHEMIC_MODULE
    CHEMIC_MODULE mod_chemic
#endif
    integer(ip), intent(in) :: itask
#ifdef CHEMIC_MODULE
    call chemic_main(itask)
#else
    call runend("CHEMIC IS NOT AVAILABLE! PLEASE COMPILE THE CHEMIC LIBRARY")
#endif
  end subroutine

  subroutine turbul(itask)
#ifdef TURBUL_MODULE
    TURBUL_MODULE mod_turbul
#endif
    integer(ip), intent(in) :: itask
#ifdef TURBUL_MODULE
    call turbul_main(itask)
#else
    call runend("TURBUL IS NOT AVAILABLE! PLEASE COMPILE THE TURBUL LIBRARY")
#endif
  end subroutine
  
#ifndef SOLFE2_MODULE

  subroutine solfe2(itask)
    integer(ip), intent(in) :: itask
  end subroutine

#endif

#ifndef INSITU_MODULE

  subroutine insitu(itask)
    integer(ip), intent(in) :: itask
  end subroutine

#endif

#ifndef NASTAL_MODULE

  subroutine nastal(itask)
    integer(ip), intent(in) :: itask
  end subroutine

#endif

#ifndef NEUTRO_MODULE

  subroutine neutro(itask)
    integer(ip), intent(in) :: itask
  end subroutine

#endif

#ifndef CASIMI_MODULE

  subroutine casimi(itask)
    integer(ip), intent(in) :: itask
  end subroutine

#endif

#ifndef HELMOZ_MODULE

  subroutine helmoz(itask)
    integer(ip), intent(in) :: itask
  end subroutine

#endif

#ifndef LEVELS_MODULE

  subroutine levels(itask)
    integer(ip), intent(in) :: itask
  end subroutine

#endif

#ifndef LATBOL_MODULE

  subroutine latbol(itask)
    integer(ip), intent(in) :: itask
  end subroutine

#endif

#endif

end module mod_module_interface
