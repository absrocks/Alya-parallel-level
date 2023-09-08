subroutine elsest_memerr(itask,vanam,vacal,istat)
  !-----------------------------------------------------------------------
  !****f* elsest_memerr
  ! NAME
  !    elsest_memerr
  ! DESCRIPTION
  !    This routine ends Alya when an error has been found
  !    allocating or deallocating memory.
  ! USES
  ! USED BY
  !    mod_elsest
  !***
  !-----------------------------------------------------------------------
  use def_elsest, only       : ip,elsest_intost
  implicit none
  integer(ip),   intent(in) :: itask,istat
  character*(*), intent(in) :: vanam,vacal

  if(itask==0) then
     call elsest_runend(&
          trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE ALLOCATED.'&
          //' RUN TIME ERROR: '//elsest_intost(istat))
  else if(itask==1) then
     call elsest_runend(trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE REALLOCATED')
  else
     call elsest_runend(trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE DEALLOCATED')
  end if

end subroutine elsest_memerr
