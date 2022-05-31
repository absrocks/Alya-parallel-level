! Created by  on 23/06/2020.

module mod_perf

  use def_kintyp,  only : ip, rp
  use mod_strings, only : lower_case
  use def_master

  implicit none

  private
  character(150) :: fil_perf
  character(150) :: pname
  character(150) :: phost
  character(150) :: pdate
  integer(ip) :: pcpus

  public :: init_perf, write_perf_header, write_perf_line

contains

  subroutine init_perf()
    use mod_iofile,         only : iofile_open_unit
    fil_perf = adjustl(trim(namda))//'-performance.csv' ! Unit 52
    call iofile_open_unit(lun_perf,fil_perf,'PERFORMANCE')
    pname = title
    call date_and_time(DATE=pdate)
    call get_environment_variable("HOST",phost)
    if (phost == "") then
      phost = "unknown"
    end if
    pcpus = npart+1

  end subroutine init_perf

  subroutine write_perf_header()
    write(lun_perf,'(a)') 'name, hostname, date, cores, module, function, type, value'
  end subroutine write_perf_header

  subroutine write_perf_line(pmodule, pfunction, ptype, pvalue)
    character(len=*), intent(in) :: pmodule
    character(len=*), intent(in) :: pfunction
    character(len=*), intent(in) :: ptype
    real(rp),         intent(in) :: pvalue
    ! 1rst line: metadata, 2nd line: value
    write(lun_perf, 1) trim(pname), trim(phost), trim(pdate), trim(intost(pcpus)), &
      trim(lower_case(pmodule)), trim(pfunction), trim(ptype), pvalue

    1 format(a, ', ', a, ', ', a, ', ', a, ', ' , a, ', ' , a, ', ', a, ', ', e15.8e3)
  end subroutine write_perf_line

end module mod_perf
