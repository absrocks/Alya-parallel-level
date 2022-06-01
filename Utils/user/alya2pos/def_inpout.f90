module def_inpout

  !-----------------------------------------------------------------------
  !    
  !    This common block contains the parameters needed for some input-
  !    output operations
  ! 
  !-----------------------------------------------------------------------
  use def_kintyp
  !
  ! Listen files
  !
  integer(ip)            :: lisda,lisre,lispa,nunit
  integer(ip), parameter :: lisin = 10, lisi1=9
  !
  ! Listen parameters:
  !
  integer(ip), parameter :: maxwp=60
  !
  ! Listen:
  !
  character(251)         :: ccard,wname             ! it must be 150+1
  integer(ip)            :: nwopa,nnwor,nnpar,endst
  real(rp)               :: param(maxwp)
  character(5)           :: words(maxwp)

contains

  function getint(fword,defau,textvari)    
    !-----------------------------------------------------------------------
    !     
    !     This routine gets the integer value associated with fword
    !
    !      - texts(1:1).eq.'!' => compulsory parameter
    !      - texts(1:1).eq.'*' => compulsory parameter
    !
    !-----------------------------------------------------------------------
    implicit         none    
    character(5)  :: fword
    integer(ip)   :: defau
    character(len = *)                :: textvari
    character(len = len(textvari))    :: texts
    integer(ip)   :: value
    integer(ip)   :: iword
    logical(lg)   :: found
    character(1)  :: markc
    integer(ip)   :: getint

    value=defau 
    found=.false.
    iword=0
    texts=textvari
    markc=texts(1:1)

    do while((iword<nwopa).and.(.not.found))
       iword=iword+1
       if(words(iword)==fword) then
          found=.true.
          value=int(param(iword))
       end if
    end do

    if(((markc=='!').or.(markc=='*')).and.(.not.found)) then
       write(lisre,1) fword,texts(2:35)
       call runend('GETINT: COMPULSORY PARAM. NOT FOUND')
    end if

    getint=value

1   format(/,4x,'*** ERROR: ',a5,' = ',a34,/,&
         & 15x,'IS A COMPULSORY PARAMETER. SPECIFY IT !')

  end function getint

  function getrea(fword,defau,textvari)
    !-----------------------------------------------------------------------
    !     
    !     This routine gets the real value associated with fword
    !
    !      - texts(1:1).eq.'!' => compulsory parameter
    !      - texts(1:1).eq.'*' => compulsory parameter
    !
    !-----------------------------------------------------------------------
    implicit         none
    character(5)  :: fword
    real(rp)  :: defau
    character(len = *)                :: textvari
    character(len = len(textvari))    :: texts
    integer(ip)   :: iword
    logical(lg)   :: found
    real(rp)      :: value
    character(1)  :: markc
    real(rp)      :: getrea

    value=defau 
    found=.false.
    iword=0
    texts=textvari
    markc=texts(1:1)

    do while((iword<nwopa).and.(.not.found))
       iword=iword+1
       if(words(iword)==fword) then
          found=.true.
          value=param(iword)
       end if
    end do

    if(((markc=='!').or.(markc=='*')).and.(.not.found)) then
       write(lisre,1) fword,texts(2:35)
       call runend('GETREA: COMPULSORY PARAM. NOT FOUND')
    end if

    getrea=value

1   format(/,4x,'*** ERROR: ',a5,' = ',a34,/,&
         & 15x,'IS A COMPULSORY PARAMETER. SPECIFY IT !')

  end function getrea

  function getcha(fword,defau,textvari)    
    !-----------------------------------------------------------------------
    !     
    !     This routine gets the character value associated with fword
    !
    !      - texts(1:1).eq.'!' => compulsory parameter
    !      - texts(1:1).eq.'*' => compulsory parameter
    !
    !-----------------------------------------------------------------------
    implicit         none    
    character(5)  :: fword
    character(5)  :: defau
    character(len = *)                :: textvari
    character(len = len(textvari))    :: texts
    character(5)  :: value
    integer(ip)   :: iword
    logical(lg)   :: found
    character(1)  :: markc
    character(5)  :: getcha

    value=defau 
    found=.false.
    iword=0
    texts=textvari
    markc=texts(1:1)

    do while((iword<nwopa).and.(.not.found))
       iword=iword+1
       if(words(iword)==fword) then
          found=.true.
          value=words(iword+1)
       end if
    end do

    if(((markc=='!').or.(markc=='*')).and.(.not.found)) then
       write(lisre,1) fword,texts(2:35)
       call runend('GETCHA: COMPULSORY PARAM. NOT FOUND')
    end if

    getcha=value

1   format(/,4x,'*** ERROR: ',a5,' = ',a34,/,&
         & 15x,'IS A COMPULSORY PARAMETER. SPECIFY IT !')

  end function getcha

  function exists(fword)
    !-----------------------------------------------------------------------
    !     
    !     This routine verifies if fword exists in words
    !
    !-----------------------------------------------------------------------
    implicit none
    character(5) :: fword
    integer(ip)  :: iword
    logical(lg)  :: exists

    exists=.false.  
    iword=0
    do while((iword<nwopa).and.(.not.exists))
       iword=iword+1
       if(words(iword)==fword) then
          exists=.true.
       end if
    end do

  end function exists

  subroutine upcase(word)
    !-----------------------------------------------------------------------
    !
    !     This routine converts wopos to upper case 
    !
    !-----------------------------------------------------------------------
    implicit none
    character(5), intent(inout) :: word
    integer(ip)                 :: iposi,ioctv

    do iposi=1,5                                   ! process all positions
       ioctv=ichar(word(iposi:iposi))               ! octal value
       if(o'141'<=ioctv.and.o'172'>=ioctv) then ! it is a lower case
          ioctv=ioctv-o'40'                          ! equivalent upper case
          word(iposi:iposi)=char(ioctv)              ! convert it to upcase
       end if
    end do ! iposi=1,5

  end subroutine upcase

  subroutine decod1(nstri,strin,lflag,digit)
    !-----------------------------------------------------------------------
    !
    !     This routine decodified a string
    !
    !-----------------------------------------------------------------------
    implicit none
    integer(ip)    :: nstri,lflag,istri,decim
    real(rp)       :: digit
    character(*)   :: strin                           ! (nstri+1)
    character(251) :: stri1

    lflag=0                                           ! It is a parameter
    istri=1                                           
    do while(istri<=nstri)                          
       decim=ichar(strin(istri:istri))                ! decimal value
       if (decim<48.or.decim>57) then                 ! It is not a num.
          if (&
               decim/= 43.and.decim/= 45.and. &       ! discard + -
               decim/= 68.and.decim/= 69.and. &       ! discard D E
               decim/=100.and.decim/=101.and. &       ! discard d e
               decim/= 46) then                       ! discard .
             lflag=1
             istri=nstri
          end if
       end if
       istri=istri+1
    end do

    if (lflag==0) then                                ! It's a number
       do istri=1,nstri
          stri1(istri:istri)=strin(istri:istri)
       end do
       stri1(nstri+1:nstri+1)=' '
       read(stri1(1:nstri),*)digit                             ! DIGIT<-'STRIN'
    end if

  end subroutine decod1

  function getpos(fword)    
    !-----------------------------------------------------------------------
    !     
    !     This routine gets the integer value associated with fword
    !
    !      - texts(1:1).eq.'!' => compulsory parameter
    !      - texts(1:1).eq.'*' => compulsory parameter
    !
    !-----------------------------------------------------------------------
    implicit         none    
    character(5)  :: fword
    integer(ip)   :: iword
    logical(lg)   :: found
    integer(ip)   :: getpos

    found=.false.
    iword=0

    do while((iword<nwopa).and.(.not.found))
       iword=iword+1
       if(words(iword)==fword) then         
          found=.true.
       end if
    end do

    if( found ) then
       getpos=iword
    else
       call runend('GETPOS: CANNOT FIND WORD FILTER')
    end if

  end function getpos
 
  subroutine livinf(itask,messa,inume)
    implicit none 
    integer(ip)  :: itask,inume
    character(*) :: messa

    if( itask == 0 ) then
       write(6,99)
    else if( itask == 1 ) then
       write(6,1) trim(messa)
    else if( itask == 2 ) then
       write(6,2) trim(messa)
    else if( itask == 3 ) then
       write(6,3) trim(messa)
    end if

99   format('--|')
1   format('--| ALYA2POS: ',a)
2   format('--| ALYA2POS: START READING ',a)
3   format('--| ALYA2POS: END   READING ',a)

  end subroutine livinf

end module def_inpout

