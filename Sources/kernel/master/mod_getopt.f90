 !> \file getopt.f90  Procedures for a getopt and getopt_long implementation to parse command-line parameters in Fortran
 !!
 !! \example getopt_example.f90
 !! \example getopt_long_example.f90
 
 
 !  Copyright (c) 2002-2017  Marc van der Sluys - marc.vandersluys.nl
 !   
 !  This file is part of the libSUFR package, 
 !  see: http://libsufr.sourceforge.net/
 !   
 !  This is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
 !  by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 !  
 !  This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 !  
 !  You should have received a copy of the GNU General Public License along with this code.  If not, see 
 !  <http://www.gnu.org/licenses/>.
 
 
 
 
 
 !***********************************************************************************************************************************
 !> \brief  Procedures for a getopt and getopt_long implementation to parse command-line parameters in Fortran
 !! 
 !! Getopt implementation in Fortran, similar but not identical to the GlibC implementation.  Two possibilities are available:
 !! - function getopt(): parses short options only
 !! - function getopt_long(): parses short and long options
 !!
 !! Both functions return the letter of the option (e.g. 'a' for -a) and set the variable optArg which contains the argument if
 !!   required and present.  In addition, getopt_long() sets the variable longOption which contains the full option found (e.g.
 !!   "-a" or "--all").  Both functions and variables can be accessed through the SUFR_getopt module.
 !!
 !! In addition, a quick help output can be generated (e.g. when no option, a help option or a non-existing option is given) using:
 !! - subroutine getopt_help():       prints a help list of all short options and their required arguments
 !! - subroutine getopt_long_help():  prints a help list of all short/long options, their required arguments and their descriptions
 !!
 !!
 !! \remarks
 !! This implementation of getopt was inspired by:
 !! - [1] GlibC getopt:  man 3 getopt  or  https://www.gnu.org/software/libc/manual/html_node/Getopt.html
 !! - [2] The getopt_long_module by Joe Krahn:  http://fortranwiki.org/fortran/show/getopt_long_module
 !! - [3] The getoptions module by Dominik Epple:  http://www.dominik-epple.de/getoptions/
 !! 
 !! The general idea comes from [1] and [2], while I followed [3] for the return values indication an issue (>!.).  Unlike [3],
 !!   I wanted both short and long options, allow the argument to be glued to the option (e.g. -ffile.txt) and get a warning if
 !!   a required argument is missing.  Unlike [2] I thought a non-OOP solution might be simpler.  In addition, I wanted to allow
 !!   an equal sign in e.g. --file=file.txt, and to provide a short description for each option, and to simplify the generation
 !!   of an explanatory list of options, which is provided through getopt_help() and getopt_long_help().
 
module mod_getopt

  use def_kintyp, only : ip
  implicit none

  integer(4), parameter, private :: longoptlen = 99        ! Maximum length of a long option (without '--')

  integer(ip)                    :: kfl_check_data_file 
  !> \brief The option's argument, if required and present
  character :: optarg*(999)
  !> \brief The short or long option found, including leading dash(es)
  character :: longoption*(longoptlen+2)
  !> \brief The current option count
  integer(4), save :: optcount = 0

  !> \brief Struct to define short and long options for getopt_long()
  type getopt_t
     !> \brief The short option (single character, without the leading dash)
     character :: short             = ''
     !> \brief The long option (without the leading dashes, max 99 characters long)
     character :: long*(longoptlen) = ''
     !> \brief Argument required? 0-no, 1-yes
     integer(4)   :: reqarg            = 0
     !> \brief A (short) description (recommended: <1 screen width; max 999 characters)
     character :: descr*(999)       = ''
  end type getopt_t

contains

  !*********************************************************************************************************************************
  !> \brief  Parse a command-line parameter and return short options and their arguments.  A warning is printed to stderr if an 
  !!            argument is required but not present.
  !!
  !! \param  optStr   String defining the allowed short options.  Characters followed by a colon have a required argument.
  !!                    Example: 'ab:c' for -a, -b <arg> and -c.
  !!
  !! \retval getopt  Either:
  !!                 - a '>' if no further command-line parameters were found
  !!                 - a '!' if an unidentified option was found
  !!                 - a '.' if a non-option argument was found
  !!                 - a single character identifying the short version of the option identified (without the leading dash)
  !!
  !! The following 'global' variable is set (accessible through the module SUFR_getopt):
  !! - optArg: the argument, if required and present

  function getopt(optStr)
    implicit none
    character, intent(in) :: optStr*(*) 
    integer(4) :: Narg, optStrI
    character :: getopt, option, arg*(999)
    logical :: found

    optcount = optcount+1

    ! Default values:
    getopt = ''
    optarg = ''

    narg = command_argument_count()
    if(optcount.gt.narg) then
       getopt = '>'
       return
    end if

    call get_command_argument(optcount, arg)

    if(arg(1:1).eq.'-') then                              ! Found a short option

       option = arg(2:2)                                  ! The short-option character
       found = .false.

       do optstri=1,len(optstr)                           ! Loop over all defined options for a match
          if(optstr(optstri:optstri).eq.option) then      ! Current option matches character in option string
             found = .true.

             if(optstr(optstri+1:optstri+1).eq.':') then  ! Option requires an argument
                if(len_trim(arg).gt.2) then               ! Argument is glued to option (no space)
                   optarg = trim(arg(3:))

                else                                      ! Next parameter should be an argument

                   optcount = optcount+1
                   call get_command_argument(optcount, optarg)
                   if(optcount.gt.narg .or. optarg.eq.'') write(0,'(A)') 'WARNING: option -'//option//' requires an argument'
                end if

             end if  ! Argument required

             exit
          end if  ! Match
       end do  ! optStrI


       if(found) then
          getopt = option
       else
          getopt = '!'
          optarg = arg
       end if

    else  ! no '-'

       getopt = '.'
       optarg = arg
    end if

  end function getopt
  !*********************************************************************************************************************************



  !*********************************************************************************************************************************
  !> \brief  Parse a command-line parameter and return short and/or long options and their arguments.  A warning is printed to 
  !!            stderr if an argument is required but not present or vice versa.
  !! 
  !! \param  longopts     Long-options used for getopt_long().  This is an array of getopt_t structs, each of which contains:
  !!                      - short:   the short option (single character, without the leading dash)
  !!                      - long:    the long option (without the leading dashes, max 99 characters long)
  !!                      - reqArg:  argument required? 0-no, 1-yes
  !!                      - descr:   a (short) description (recommended: <1 screen width; max 999 characters)
  !!                      An example entry would be getopt_t('f','file',1,'Input file name').
  !! 
  !! \retval getopt_long  Either:
  !!                      - a '>' if no further command-line parameters were found
  !!                      - a '!' if an unidentified option was found
  !!                      - a '.' if a non-option argument was found
  !!                      - a single character identifying the short version of the option identified (without the leading dash)
  !! 
  !! The following 'global' variables are set (accessible through the module SUFR_getopt):
  !! - longOption:  the short or long option found, including leading dash(es)
  !! - optArg:      the option's argument, if required and found
  !!

  function getopt_long(longopts)
    implicit none
    type(getopt_t), intent(in) :: longopts(:)
    integer(4) :: Narg, optI, pos,                                                  debug=0  ! 0-1
    character :: getopt_long, option, arg*(999), longOpt*(longoptlen)
    logical :: found, hasEql

    optcount = optcount+1

    ! Default values:
    getopt_long = ''
    optarg      = ''
    longoption  = ''

    ! Get the current command-line parameter:
    narg = command_argument_count()
    if(optcount.gt.narg) then
       getopt_long = '>'
       return
    end if

    call get_command_argument(optcount, arg)
    if(debug.ge.1) write(*,'(A,I0,A)') 'getopt_long():  option ', optcount, ': '//trim(arg)


    ! Check for long options, short options and arguments:
    if(arg(1:2).eq.'--') then                                           ! Found a long option
       longopt = trim(arg(3:))

       ! Allow for an argument connected through an equal sign, e.g. --file=file.txt
       haseql = .false.
       pos = scan(trim(longopt), '=')
       if(pos.gt.0) then  ! Separate option and argument
          optarg = trim(longopt(pos+1:))
          longopt = longopt(1:pos-1)
          haseql = .true.
       end if

       found = .false.
       do opti=1,size(longopts)
          if(longopts(opti)%long.eq.longopt) then                       ! Current option matches character in option string
             found = .true.
             longoption = '--'//trim(longopt)
             if(longopts(opti)%reqArg.gt.0 .and. .not.haseql) then      ! Option requires an argument, not glued using =

                optcount = optcount+1
                call get_command_argument(optcount, optarg)

                if(optcount.gt.narg .or. optarg.eq.'') write(0,'(A)') 'WARNING: option --'//option//' requires an argument'

             else if(longopts(opti)%reqArg.eq.0 .and. haseql) then
                write(0,'(A)') 'WARNING: option --'//option//' does not require an argument'
             end if

             exit
          end if
       end do

       if(found) then
          getopt_long = longopts(opti)%short
       else
          getopt_long = '!'
          optarg = arg
       end if


    else if(arg(1:1).eq.'-') then  ! Short option
       option = arg(2:2)

       found = .false.
       do opti=1,size(longopts)
          if(longopts(opti)%short.eq.option) then            ! Current option matches character in option string
             found = .true.
             longoption = '-'//option
             if(longopts(opti)%reqArg.gt.0) then             ! Option requires an argument
                if(len_trim(arg).gt.2) then                  ! Argument is glued to option (no space)
                   optarg = trim(arg(3:))

                else  ! Next parameter should be argument

                   optcount = optcount+1
                   call get_command_argument(optcount, optarg)

                   if(optcount.gt.narg .or. optarg.eq.'') write(0,'(A)') 'WARNING: option -'//option//' requires an argument'

                end if
             end if  ! Argument

             exit
          end if
       end do

       if(found) then
          getopt_long = option
       else
          getopt_long = '!'
          optarg = arg
       end if

    else  ! no '-'

       getopt_long = '.'
       optarg = arg
    end if

    if(debug.ge.1) write(*,'(2(A,I0))') 'optCount: ',optcount, ' -> ', optcount+1

  end function getopt_long
  !*********************************************************************************************************************************



  !*********************************************************************************************************************************
  !> \brief  Print a help list of all short options and their required arguments
  !!
  !! \param optStr  Short-options string used for getopt()

  subroutine getopt_help(optStr)
    implicit none
    character, intent(in) :: optStr*(*)
    character :: curChar
    integer(4) :: ichar4

    write(*,'(A)', advance='no') 'Available options: '
    do ichar4=1,len_trim(optstr)
       curchar = optstr(ichar4:ichar4)

       if(curchar==':') then
          write(*,'(A)', advance='no') ' <arg>'
       else
          if(ichar4>1) write(*,'(A)', advance='no') ','
          write(*,'(A)', advance='no') ' -'//curchar
       end if
    end do
    write(*,'(A)') '.'

  end subroutine getopt_help
  !*********************************************************************************************************************************


  !*********************************************************************************************************************************
  !> \brief  Print a help list of all short/long options, their required arguments and their descriptions
  !!
  !! \param longopts  Long-options struct used for getopt_long()
  !! \param lineBef   Number of lines to print before option list (default: 0)
  !! \param lineAft   Number of lines to print after option list (default: 0)

  subroutine getopt_long_help(longopts, lineBef, lineAft)
    implicit none
    type(getopt_t), intent(in) :: longopts(:)
    integer(4), intent(in), optional :: lineBef, lineAft
    type(getopt_t) :: curOpt
    integer(4) :: iLine, iOpt, nchar4, iSpc

    if(present(linebef)) then
       if(linebef.gt.0) then
          do iline=1,linebef
             write(*,*)
          end do
       end if
    end if

    write(*,'(A)') 'Available options:'
    do iopt=1,size(longopts)
       curopt = longopts(iopt)

       nchar4 = 0
       ! Print short option if defined:
       if(trim(curopt%short).ne.'') then
          write(*,'(A4)', advance='no') '  -'//curopt%short
          nchar4 = nchar4 + 4
       end if

       ! Print long option if defined:
       if(trim(curopt%long).ne.'') then
          write(*,'(A)', advance='no') '  --'//trim(curopt%long)
          nchar4 = nchar4 + len_trim(curopt%long) + 3  ! max: 99+3 = 102
       end if

       ! Print argument if required:
       if(curopt%reqArg.gt.0) then
          write(*,'(A7)', advance='no') '  <arg>'
          nchar4 = nchar4 + 7
       end if

       ! Justify description using spaces:
       do ispc=1,30-nchar4
          write(*,'(1x)', advance='no')
       end do

       ! Print description:
       write(*,'(5x,A)') trim(curopt%descr)
    end do

    if(present(lineaft)) then
       if(lineaft.gt.0) then
          do iline=1,lineaft
             write(*,*)
          end do
       end if
    end if

  end subroutine getopt_long_help
  !*********************************************************************************************************************************


end module mod_getopt
 !***********************************************************************************************************************************
