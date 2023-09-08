module mod_ecoute

  use def_kintyp, only : ip,rp
  use def_inpout
  use def_master, only : fil_pdata, kfl_split_plus,intost

  implicit none

  private

  public :: ecoute

contains

  subroutine ecoute(subna,STOP_END_OF_FILE, DO_NOT_READ_INCLUDE)

    !-----------------------------------------------------------------------
    !
    ! Reads a string and interprets it as words and parameters.
    !
    !   - Maximum number of words and parameters = maxwp.
    !   - Only the first five characters of each word are decoded.
    !   - The underline characters '_' are discarted.
    !   - Lower case letters are converted to upper case.
    !   - Each word or parameter must be separated by ' ', '=', ':' or ','
    !   - A "comment" begins with '$', '!', '/'  or '/' in any place.
    !   - A line that has a comment beginning with '/' or '/'
    !     continuates in the next line.
    !   - A line beginning with title is not decoded. title directive.
    !   - A line beginning with include is an include directive.
    !   - A line beginning with echo turns 'on' or 'off' the echo.
    !
    !-----------------------------------------------------------------------
    !
    ! Var
    !
    character(*), intent(in)           :: subna
    logical(lg),  intent(in), optional :: STOP_END_OF_FILE
    logical(lg),  intent(in), optional :: DO_NOT_READ_INCLUDE
    real(rp)                           :: digit
    integer(ip)                        :: first,firsp,i,last,lastp,ptrwo,npptr,nwptr
    integer(ip)                        :: leng,flag,resum
    logical(lg)                        :: newline=.false.
    logical(lg), save                  :: echo=.false.     ! default echo off. to change it use: echo on
    logical(lg)                        :: stop_end_of_file_opt
    logical(lg)                        :: do_not_read_include_opt
    !
    ! Options
    !
    stop_end_of_file_opt = .true.
    if( present(STOP_END_OF_FILE) ) then
       stop_end_of_file_opt = STOP_END_OF_FILE
    end if
    do_not_read_include_opt = .false.
    if( present(DO_NOT_READ_INCLUDE) ) then
       do_not_read_include_opt = DO_NOT_READ_INCLUDE
    end if
    !
    ! Data
    !
    if(lispa==0) then
       nunit = lisda                                  !  initial data file
       lispa = 1
    end if
    !
    ! Begin
    !
    ccard=' '
    nnwor=0                                           ! initialize.
    nnpar=0
    nwptr=0
    npptr=0
    resum=0
    do i=1,maxwp
       words(i)=' '
       param(i)=0.0_rp
    end do

    !
    ! Binary ecoute reading
    !
99  continue
    if (kfl_binin == 2) then
       call ecoute_bin(subna)
       ! End of binary file: close include file
       if (kfl_binin == 3) then
          ! Disable binary mode and close file
          kfl_binin = 0
          ! Roundabout to bug?
          newline=.false.
          last=0
          lastp=0
          go to 101
       end if

       ! Return with data
       return
    end if

100 continue
    do while(((nnwor==0).and.(nnpar==0)&              ! don't return without answer
         .or.newline).or.resum==1)                    ! continue reading if / or \

       if (resum==0) then
          newline=.false.                             ! initialize.
          last=0
          lastp=0
       end if
       firsp=1
       resum=0
       read(nunit,10,end=101,err=1) ccard             ! read a card
       if( subna(1:4) == 'DONT' ) return
       !     leng=lnbln1(ccard)                             ! calculate the length.
       leng=len_trim(ccard)                             ! calculate the length.

       decode_card: do while(last<leng)               ! decode all the card.
          first=last+1
          loop_first: do while(            &
               ccard(first:first)=='_'.or. &          ! jump null character (_)
               ccard(first:first)==' '.or. &          ! jump separators ( =:,)
               ccard(first:first)=='='.or. &
               ccard(first:first)==':'.or. &
               ccard(first:first)==','.or. &
               iachar(ccard(first:first)) == 9 .or. &         ! Tab character ASCII code is 9
               (kfl_split_plus==1 .and. ccard(first:first)=='+') )  !!FER for CHEMIC
             first=first+1
             if(first>leng) exit loop_first
          end do loop_first
          if(last==0) firsp=first                     ! save first to print card
          last=first
          loop_last: do while(             &
               ccard(last:last)/=' ' .and. &          ! look for separator ( =:,).
               ccard(last:last)/='=' .and. &
               ccard(last:last)/=':' .and. &
               ccard(last:last)/=',' .and. &
               ccard(last:last)/='$' .and. &          ! look for coment ($!/\).
               ccard(last:last)/='!' .and. &
               ccard(last:last)/='\' .and. &
               ccard(last:last)/='\' .and. &
               iachar(ccard(last:last)) /= 9  .and. &      ! Tab character ASCII code is 9     )
               (kfl_split_plus/=1 .or. ccard(last:last)/='+') ) !!FER for CHEMIC, use De Morgan´s laws to understand
             last=last+1
             if(last>leng) exit loop_last
          end do loop_last

          if(last<=251            .and.(&
               ccard(last:last)=='$'.or.&
               ccard(last:last)=='!')) leng=last-1
          if(last<=251            .and.(&
               ccard(last:last)=='\'.or.&
               ccard(last:last)=='\')) then
             leng=last-1                              ! deal with continuation
             newline=.true.                           ! set new line flag.
          end if
          last=last-1
          if(last>=first) then                        ! is it a word or a parameter
             lastp=last                               ! save last to print cardlogic
             call decod1(last-first+1_ip,&
                  ccard(first:last),flag,digit)
             wname= ccard(first:last)                 ! keep the useful ccard part in wname
             if(flag==0) then                         ! it is a parameter.
                nnpar=nnpar+1                         ! # of parameters
                npptr=npptr+1                         ! integer :: to next parameter
                if(npptr>maxwp) go to 4               ! error.
                if(nwptr>npptr) npptr=nwptr
                param(npptr)=digit
             else                                     ! it is a word.
                nnwor=nnwor+1                         ! # of words
                nwptr=nwptr+1                         ! integer :: to next word
                if(nwptr>maxwp) go to 5               ! error.
                if(npptr>=nwptr) nwptr=npptr+1
                ptrwo=1
                do while ((first<=last).and.(ptrwo<=5))
                   words(nwptr)(ptrwo:ptrwo)=ccard(first:first)
                   ptrwo=ptrwo+1
                   first=first+1
                   do while (ccard(first:first)=='_') ! jump null character
                      first=first+1
                   end do
                end do
                call upcase(words(nwptr))             ! convert to upper case.
             end if
          end if ! (last>=first)
          if((nnwor==1).and.(nnpar==0).and.&          ! deal with title or include
               ((words(1)=='TITLE').or.((words(1)=='INCLU').and.&
               (subna/='NOREAD')))) then
             if(echo.and.(subna/='noecho'))&
                  write(lisre,20) 'ecoute',ccard(firsp:leng)
             last=last+2
             do while(ccard(last:last)==' ')          ! remove blank spaces
                last=last+1
             end do
             if(leng<last) go to 6                    ! error
             ccard=ccard(last:leng)                   ! remove words(1) from ccard
             leng=leng-last+1
             if(words(1)=='TITLE') then               ! deal with titles
                firsp=1
                lastp=leng
             else if (.not.do_not_read_include_opt) then     ! deal with include directive
                !if(nunit==lisin) go to 3             ! error
                last=1                                ! remove tail comments
                do while((last<=leng).and.&           ! look for end (last=leng)
                     ccard(last:last)/=' ')           ! look for separator ( )
                   last=last+1
                end do
                ccard=ccard(1:last-1)                 ! remove tail comments
                if(nunit==lisin) then
                   nunit=lisi1
                else
                   nunit=lisin
                end if
                call opincl()
                lastp=0                               ! to ignore the echo
                nnwor=0                               ! forget all
                nwptr=0
                words(1)=' '
                if (kfl_binin==2) go to 99            ! go to ecoute_bin if binary mode
             else if (do_not_read_include_opt) then
                kfl_binin=0
             end if
             last=leng                                ! to break the do while
          else if(subna=='NOREAD') then
             newline=.false.
             nnwor=1
          end if
       end do decode_card

       if((words(1)=='ECHO') .and.&
            (nnpar==0).and.(nnwor==2)) then           ! deal with echo
          if(words(2)=='OFF') then
             echo=.false.
             if(subna/='noecho') write(lisre,20) 'ECOUTE','ECHO OFF'
          else
             echo=.true.
             if(subna/='noecho') write(lisre,20) 'ECOUTE','ECHO ON'
          endif
          nnwor=0                                     ! forget all
          nwptr=0
          do i=1,maxwp
             words(i)=' '
          end do
       else                                           ! print card
          if((echo).and.(firsp<=lastp).and.(subna/='noecho')) then
             if(newline) then
                lastp=lastp+2
                ccard(lastp-1:lastp)=' /'
             end if
             write(lisre,20) subna,ccard(firsp:lastp)
          end if
       end if

    end do ! while ((nnwor==0).and.(nnpar==0).or.newline)

    nwopa=max(npptr,nwptr)
    return
    !
    ! End of include
    !
101 continue
    if(nunit/=lisin.and.nunit/=lisi1) then
       if( endst == 1 .and. stop_end_of_file_opt ) then
          goto 2 ! error
       else
          words(1)='ENDFI'
          return
       end if
    end if
    close(unit=nunit,status='keep')
    if(nunit==lisi1) then
       nunit=lisin
    else
       nunit=lisda
    end if
    if(echo.and.(subna/='noecho'))&
         write(lisre,20) 'ECOUTE','END OF INCLUDE FILE'
    resum=1
    go to 100 ! for resume the error return to the same place.
    !
    ! Errors:
    !
1   call runend('ECOUTE: ERROR DETECTED WHEN READING')
2   call runend('ECOUTE: END OF FILE DETECTED IN SUBROUTINE '//trim(subna))
3   call runend('ECOUTE: ERROR: INCLUDE FROM INCLUDE')
4   call runend('ECOUTE: TOO MANY PARAM IN COMMAND  ')
5   call runend('ECOUTE: TOO MANY WORDS IN COMMAND  ')
6   call runend('ECOUTE: BLANK IS ILEGAL HERE       ')
7   call runend('ECOUTE: COULD NOT OPEN FILE: '//adjustl(trim(ccard)))
    !
    ! Format
    !
10  format(a250)
20  format(1x,a6,' <-- ',a)

  end subroutine ecoute


  subroutine ecoute_bin(subna)
    !-----------------------------------------------------------------------
    !
    ! THIS ROUTINE IS OBSOLETE
    !
    ! Reads a binary record/set and interprets it as parameters.
    !
    !   - Maximum number of parameters = maxwp.
    !   - ecoute is previously used to read INCLUDE statement
    !   - First read parses header of binary, following ecoute_bin calls
    !     read binary records which are composed of params(ndime)
    !   - A binary file can be composed of several sets (e.g. vel, pres ...)
    !   - Each set is composed of a header and a body of records (NPOINTS x NDIMENSIONS)
    !   - Data is stored ALWAYS in LITTLE ENDIAN format
    !
    !    /
    !   | * HEADER OF BINARY SET 1 *
    !   |
    !   | INT_FORMAT (INTEGER*1: 4-INTEGER*4/4BYTES, 8-INTEGER*8/8BYTES)
    !   | REAL_FORMAT (INTEGER*1: 0-NO_REALS, 2-INTEGER*2/2BYTES, 4-REAL*4/4BYTES, 8-REAL*8/8BYTES)
    !   | ID_PRESENT (INTEGER*1: 0-NO/1-YES)
    !   | NPOINTS (INTEGER*8)
    !   | NDIMENSIONS (INTEGER*8)
    ! S | LOWER_VALUE_DIM1 HIGHER_VALUE_DIM1 (REAL*8) <-- If REAL_FORMAT==0
    ! E | LOWER_VALUE_DIM2 HIGHER_VALUE_DIM2 (REAL*8) <-- As many as NDIMENSIONS
    ! T | ...
    !   | ...
    ! 1 |
    !   | * BODY OF BINARY SET 1 (EACH RECORD IS A LINE)                       *
    !   | * DATA FORMAT CAN BE FIXED OR FLOAT POINT:                           *
    !   | * {INT_FORMAT}[INT_FORMAT_1 ]..[INT_FORMAT_n ] <-- If REAL_FORMAT=-1 *
    !   | * Or                                                                 *
    !   | * {INT_FORMAT}[REAL_FORMAT_1]..[REAL_FORMAT_n] <-- Otherwise         *
    !   |
    !   | {id_1}[val_1]..[val_n]
    !   | ...
    !   | ...
    !   | {id_npoints}[val_1]..[val_n]
    !    \
    !
    !   ...
    !
    !    /
    ! S |
    ! E | HEADER OF BINARY SET N
    ! T |
    !   | BODY OF BINARY SET N
    ! N |
    !    \
    !
    ! Note about compression in INTEGER*2 format using unsigned shorts:
    ! Fortran does not allow unsigned shorts, then INTEGER*2 is converted
    ! to INTEGER*4 extending left-hand bits with 0s avoiding sign extension
    ! from signed shorts to signed integers and ensuring unsigned representation
    !
    !-----------------------------------------------------------------------
    !
    ! Var
    !
    character(*)         :: subna
    logical(lg), save    :: echo=.false.      ! default echo off. to change it use: echo on
    logical(lg), save    :: header = .false.  ! signals if header has been read

    ! Header vars
    integer(1), save      :: iform, rform      ! integer format and real format
    integer(1), save      :: idpre             ! id present as 1st element in record
    integer(8), save      :: npoin, ndime      ! number of records and dimension of each record
    real(8), save         :: lvalue(maxwp), hvalue(maxwp), rscale(maxwp) ! variable range compression vars for INT*2

    ! Record vars
    integer(8)            :: idime
    integer(4), save      :: ipoin4            ! id counter for ID_PRESENT=0
    integer(8), save      :: ipoin8            ! id counter for ID_PRESENT=0
    integer(4), parameter :: USHRT_MAX = 65536 ! Number of ranges in unsigned short
    integer(4)            :: parami4(maxwp)    ! INT*4  format      (iform=4 && rform=0)
    integer(8)            :: parami8(maxwp)    ! INT*8  format      (iform=8 && rform=0)
    integer(2)            :: paramr2(maxwp)    ! INT*2  compression (rform=2)
    real(4)               :: paramr4(maxwp)    ! REAL*4 format      (rform=4)
    real(8)               :: paramr8(maxwp)    ! REAL*8 format      (rform=8)
    !
    ! HEADER SECTION
    !
    ! Read header for initial read
    if (.not. header) then
       ! Precision and dimensions (NPOIN x NDIME)
       read(nunit,err=110,end=120) iform, rform, idpre, npoin, ndime

       ! Error checking
       ! Integer and real range
       if ((iform/=4).and.(iform/=8))      go to 111
       if ((rform/=0).and.(rform/=2)&
            .and.(rform/=4).and.(rform/=8)) go to 112

       ! Integer and real precision matching with Alya
       if (iform > ip)                     go to 113
       if (rform > rp)                     go to 114

       ! ID range
       if ((idpre/=0).and.(idpre/=1))      go to 115

       ! NPOINT and NDIMENSION range
       if (npoin<=0)                       go to 116
       if (ndime<=0)                       go to 117


       ! Initialize id poin counter
       if (idpre == 0) then
          ipoin8 = 0
       end if

       ! Compressed binary data in INTEGER*2 format
       if (rform == 2) then
          ! Read value range for each dimension
          do idime=1, ndime
             read(nunit,err=110,end=120) lvalue(idime), hvalue(idime)
             ! Check range
             if (lvalue(idime)>hvalue(idime)) go to 118

             ! Compute real scale
             rscale(idime) = (hvalue(idime) - lvalue(idime))/USHRT_MAX
          end do
       end if

       ! Set header as read
       header = .true.
       if (echo) write(*,*) 'ECOUTE_BIN: HEADER --->', 'INT_FORMAT:', iform, 'REAL_FORMAT:', rform,&
            'ID_PRESENT:', idpre, 'NPOIN:', npoin, 'NDIME:', ndime,&
            'LOWER_VALUE:', lvalue, 'HIGHER_VALUE:', hvalue, 'RSCALE:', rscale
    end if


    !
    ! ID SECTION
    !
    ! Read binary data with ipoin id
    ! param(1) is left for idpoin
    !
    ! {id_1}[val_1]..[val_n]
    if (idpre == 1) then

       ! Read id poin and convert to REAL_RP
       if (iform == 4) then
          read(nunit,err=110,end=120) ipoin4
          param(1) = REAL(ipoin4,rp)
       else
          read(nunit,err=110,end=120) ipoin8
          param(1) = REAL(ipoin8,rp)
       end if

       ! Read binary data without ipoin id
       !
       ! [val_1]..[val_n]
    else

       ! Advance id poin
       ipoin8 = ipoin8 + 1

       ! Convert ipoin to REAL_RP
       param(1) = REAL(ipoin8,rp)
    end if


    !
    ! VALUE SECTION
    !
    ! NO REALS in record
    if (rform == 0) then

       ! INTEGER*4 in record
       if (iform == 4) then
          read(nunit,err=110,end=120) (parami4(idime),idime=1,ndime)
          do idime=1, ndime
             param(idime+1) = REAL(parami4(idime),rp)
          end do

          ! INTEGER*8 in record
       else
          read(nunit,err=110,end=120) (parami8(idime),idime=1,ndime)
          do idime=1, ndime+1
             param(idime+1) = REAL(parami8(idime),rp)
          end do
       end if

       ! Reals in record
       ! INTEGER*2 format compression
    else if (rform == 2) then
       read(nunit,err=110,end=120) (paramr2(idime),idime=1,ndime)
       do idime=1, ndime
          call runend('MOD_ECOUTE: THIS IS TOTALLY OBSOLETE AND NOT STANDARD')
          !param(idime+1) = REAL(lvalue(idime) + (ZEXT(paramr2(idime)) * rscale(idime)),rp)
       end do

       ! REAL*4 format
    else if (rform == 4) then
       read(nunit,err=110,end=120) (paramr4(idime),idime=1,ndime)
       do idime=1, ndime
          param(idime+1) = REAL(paramr4(idime),rp)
       end do

       ! REAL*8 format
    else
       read(nunit,err=110,end=120) (paramr8(idime),idime=1,ndime)
       do idime=1, ndime
          param(idime+1) = REAL(paramr8(idime),rp)
       end do
    end if


    ! Number of parameters including id poin
    nnpar = ndime + 1

    if (echo) write(*,*) 'ECOUTE_BIN: READ --->', 'NNPAR:', nnpar, 'PARAM', (param(idime), idime=1,nnpar)

    return
    !
    ! Errors:
    !
110 call runend('ECOUTE_BIN: ERROR DETECTED WHEN READING')
111 call runend('ECOUTE_BIN: WRONG IFORM VALUE (4:INTEGER*4,8:INTEGER*8)')
112 call runend('ECOUTE_BIN: WRONG RFORM VALUE (0:NO_REALS,2:INTEGER*2,4:REAL*4,8:REAL*8)')
113 call runend('ECOUTE_BIN: INTEGER PRECISION IN BINARY IS HIGHER THAN IP')
114 call runend('ECOUTE_BIN: REAL PRECISION IN BINARY IS HIGHER THAN RP')
115 call runend('ECOUTE_BIN: WRONG IDPRE VALUE (0,1)')
116 call runend('ECOUTE_BIN: WRONG NPOIN VALUE (>0)')
117 call runend('ECOUTE_BIN: WRONG NDIME VALUE (>0)')
118 call runend('ECOUTE_BIN: WRONG LVALUE AND HVALUE RANGES (LVALUE<HVALUE)')

    !
    ! End of file:
    !
120 continue
    ! Signal closing
    header = .false.
    kfl_binin = 3
    if (echo) write(*,*) 'ECOUTE_BIN: CLOSING --->', 'HEADER:', header, 'KFL_BININ:', kfl_binin,&
         'nunit:', nunit, 'lisin:', lisin, 'lisi1:', lisi1

  end subroutine ecoute_bin

!#if  defined (__GFORTRAN__)  || defined (_MERCURIUM) || defined (_CRAYFTN) || defined(__ibmxl__)
!  function ZEXT(ushort)
!    implicit none
!    integer(2) :: ushort
!    integer(4) :: ZEXT
!
!    ZEXT = ushort
!
!    ! Extend to unsigned format
!    if (ZEXT < 0) then
!       ZEXT = ZEXT + 65536
!    endif
!  end function ZEXT
!#endif


  subroutine opincl()
    !
    ! Include file cannot be opened: look the directory of the data file
    !
    integer(ip)   :: i,istat
    integer(4)    :: istat4,nunit4
    character(20) :: wstat

    !inquire(unit=nunit,opened=lopen,exist=lexis,form=wform)
    !print*,'popo=',trim(ccard)
    !print*,'caca=',nunit,lopen,lexis
    !print*,'pipi=',trim(wform)
    nunit4 = int(nunit,4)
    if (kfl_binin==1) then
       call runend('OPINCL: OBSOLETE OPTION')
!       ! Open INCLUDE file in binary mode (data is always stored in LITTLE_ENDIAN format
!       open(unit=nunit4,file=adjustl(trim(ccard)),err=7,form='UNFORMATTED',&
!#ifdef __INTEL_COMPILER
!            access='STREAM',buffered='YES',convert='LITTLE_ENDIAN',status='OLD',IOSTAT=istat4)
!#else
!       access='STREAM',convert='LITTLE_ENDIAN',status='OLD',IOSTAT=istat4)
!#endif
!       kfl_binin = 2
    else
       open(unit=nunit4,file=adjustl(trim(ccard)),err=7,form='FORMATTED',&
#ifdef __INTEL_COMPILER
            buffered='YES',status='OLD',IOSTAT=istat4)
#else
       status='OLD',IOSTAT=istat4)
#endif
    end if
    return

7   i=len(trim(fil_pdata))
    do while(i>1)
       i=i-1
       if(fil_pdata(i:i)=='/'.or.fil_pdata(i:i)==achar(92)) i=-i
    end do
    if(i<0) then
       i=-i
       if (kfl_binin==1) then
          call runend('OPINCL: OBSOLETE OPTION')
!          ! Open INCLUDE file in binary mode (data is always stored in LITTLE_ENDIAN format
!          open(unit=nunit,file=trim(fil_pdata(1:i))//trim(ccard),err=8,&
!#ifdef __INTEL_COMPILER
!               form='unformatted',access='stream',buffered='yes',convert='LITTLE_ENDIAN',status='old')
!#else
!          form='unformatted',access='stream',convert='LITTLE_ENDIAN',status='old')
!#endif
!          kfl_binin = 2
       else
          open(unit=nunit,file=trim(fil_pdata(1:i))//trim(ccard),err=8,&
#ifdef __INTEL_COMPILER
               form='formatted',buffered='YES',status='old')
#else
          form='formatted',status='old')
#endif
       end if
       return
    else
       istat = int(istat4,ip)
       wstat = intost(istat)
       call runend('ECOUTE: COULD NOT OPEN FILE: '//trim(ccard)//'. ERROR '//trim(wstat))
    end if

8   continue
    istat = int(istat4,ip)
    wstat = intost(istat)
    call runend('ECOUTE: COULD NOT OPEN FILE: '//adjustl(trim(fil_pdata(1:i))&
         //trim(ccard))//'. ERROR '//trim(wstat))

  end subroutine opincl

end module mod_ecoute
